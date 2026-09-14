"""Versioned source connectors for the target evidence hub.

Connectors return :class:`EvidenceRecord` objects and never convert a failed
query into a negative association. Network errors are surfaced so the hub can
record them without fabricating zero scores.
"""

from __future__ import annotations

import json
import logging
import math
import re
import urllib.parse
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping

import numpy as np
import pandas as pd

try:
    from ..common.http import HttpError, http_get
except ImportError:  # pragma: no cover - direct src import in tests/scripts
    from common.http import HttpError, http_get
from .context import EvidenceContext
from .models import EvidenceRecord, EvidenceTier

LOG = logging.getLogger("evidence.connectors")


class ConnectorError(RuntimeError):
    """Raised when a source response cannot be interpreted safely."""


class BaseConnector(ABC):
    name = "base"
    source_version = "unknown"
    source_group = "base"
    tier = EvidenceTier.CURATED
    requires_network = True

    @abstractmethod
    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        raise NotImplementedError


def _json_get(
    url: str,
    *,
    params: Mapping[str, Any] | None = None,
    timeout: int = 120,
    retries: int = 3,
) -> Any:
    if params:
        query = urllib.parse.urlencode(
            {key: value for key, value in params.items() if value is not None},
            doseq=True,
        )
        separator = "&" if "?" in url else "?"
        url = f"{url}{separator}{query}"
    try:
        body = http_get(
            url,
            timeout=timeout,
            retries=retries,
            max_bytes=64 * 1024 * 1024,
        )
    except HttpError as exc:
        raise ConnectorError(str(exc)) from exc
    try:
        return json.loads(body.decode("utf-8", "replace"))
    except json.JSONDecodeError as exc:
        raise ConnectorError(f"invalid JSON from {url}: {exc}") from exc


def _post_json(
    url: str,
    payload: Mapping[str, Any],
    *,
    timeout: int = 120,
    retries: int = 3,
) -> Any:
    import urllib.request

    last_error: Exception | None = None
    for attempt in range(1, max(1, retries) + 1):
        try:
            request = urllib.request.Request(
                url,
                data=json.dumps(payload).encode("utf-8"),
                headers={
                    "Content-Type": "application/json",
                    "User-Agent": "liver-cancer-pipeline/1.7",
                },
            )
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8", "replace"))
        except Exception as exc:  # noqa: BLE001
            last_error = exc
    raise ConnectorError(f"POST {url} failed: {last_error}")


def _normalise_symbol(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9-]{1,14}", text):
        return ""
    if not (text.isupper() or re.search(r"\d", text)):
        return ""
    return text.upper()


def _float_or_none(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _score_from_pchembl(value: Any) -> float | None:
    number = _float_or_none(value)
    if number is None:
        return None
    return float(np.clip(number / 14.0, 0.0, 1.0))


def _score_from_pvalue(value: Any) -> float | None:
    number = _float_or_none(value)
    if number is None or number <= 0:
        return None
    return float(np.clip(-math.log10(number) / 50.0, 0.0, 1.0))


def _recursive_dicts(value: Any) -> Iterator[dict[str, Any]]:
    if isinstance(value, dict):
        yield value
        for nested in value.values():
            yield from _recursive_dicts(nested)
    elif isinstance(value, list):
        for nested in value:
            yield from _recursive_dicts(nested)


def _first_key(record: Mapping[str, Any], aliases: Iterable[str]) -> Any:
    normalised = {
        re.sub(r"[^a-z0-9]+", "", str(key).lower()): value
        for key, value in record.items()
    }
    for alias in aliases:
        key = re.sub(r"[^a-z0-9]+", "", alias.lower())
        if key in normalised and normalised[key] not in (None, ""):
            return normalised[key]
    return None


class OpenTargetsConnector(BaseConnector):
    name = "OpenTargets"
    source_version = "platform-api-v4"
    source_group = "opentargets"
    tier = EvidenceTier.CURATED

    ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        disease_name = str(context.disease.get("name") or "").strip()
        if not disease_name:
            return []
        search_query = """
        query searchDisease($q: String!) {
          search(queryString: $q, entityNames: ["disease"]) {
            hits { id name entity }
          }
        }
        """
        payload = _post_json(
            self.ENDPOINT,
            {"query": search_query, "variables": {"q": disease_name}},
            timeout=context.timeout_seconds,
        )
        hits = (
            (((payload.get("data") or {}).get("search") or {}).get("hits") or [])
        )
        disease = next(
            (row for row in hits if row.get("entity") == "disease"),
            None,
        )
        if not disease:
            raise ConnectorError(f"disease not found in Open Targets: {disease_name}")
        disease_id = str(disease.get("id") or "")
        records: list[EvidenceRecord] = []
        page_size = min(500, max(1, int(context.max_records_per_source)))
        for page in range(
            math.ceil(context.max_records_per_source / page_size)
        ):
            query = """
            query diseaseTargets($id: String!, $index: Int!, $size: Int!) {
              disease(efoId: $id) {
                id
                name
                associatedTargets(page: {index: $index, size: $size}) {
                  count
                  rows {
                    score
                    target { id approvedSymbol }
                  }
                }
              }
            }
            """
            response = _post_json(
                self.ENDPOINT,
                {
                    "query": query,
                    "variables": {
                        "id": disease_id,
                        "index": page,
                        "size": page_size,
                    },
                },
                timeout=context.timeout_seconds,
            )
            disease_data = (response.get("data") or {}).get("disease") or {}
            rows = (
                (disease_data.get("associatedTargets") or {}).get("rows") or []
            )
            if not rows:
                break
            for row in rows:
                target = row.get("target") or {}
                symbol = _normalise_symbol(target.get("approvedSymbol"))
                target_id = str(target.get("id") or "")
                if not symbol or not target_id:
                    continue
                records.append(
                    EvidenceRecord(
                        source=self.name,
                        source_record_id=f"{disease_id}:{target_id}",
                        evidence_type="curated_association",
                        subject_type="disease",
                        subject_id=disease_id,
                        relation="associated_with",
                        object_type="target",
                        object_id=target_id,
                        target_symbol=symbol,
                        tier=self.tier,
                        source_version=self.source_version,
                        source_group=self.source_group,
                        score=_float_or_none(row.get("score")),
                        species="Homo sapiens",
                        url=self.ENDPOINT,
                        license="Open Targets Platform terms",
                        payload={"disease_name": disease.get("name")},
                    )
                )
                if len(records) >= context.max_records_per_source:
                    break
            if len(records) >= context.max_records_per_source:
                break
        return records


class ChEMBLConnector(BaseConnector):
    name = "ChEMBL"
    source_version = "ChEMBL REST API"
    source_group = "chembl"
    tier = EvidenceTier.EXPERIMENTAL

    BASE = "https://www.ebi.ac.uk/chembl/api/data"

    def _resolve_molecule_id(self, context: EvidenceContext) -> str:
        explicit = str(context.compound.get("chembl_id") or "").strip()
        if explicit:
            return explicit
        query = str(
            context.compound.get("inchi_key")
            or context.compound.get("canonical_smiles")
            or context.compound.get("isomeric_smiles")
            or ""
        ).strip()
        if not query:
            return ""
        try:
            payload = _json_get(
                f"{self.BASE}/molecule/search.json",
                params={"q": query, "limit": 20},
                timeout=context.timeout_seconds,
            )
        except ConnectorError:
            smiles = str(
                context.compound.get("isomeric_smiles")
                or context.compound.get("canonical_smiles")
                or ""
            ).strip()
            if not smiles:
                raise
            encoded = urllib.parse.quote(smiles, safe="")
            payload = _json_get(
                f"{self.BASE}/similarity/{encoded}/40.json",
                params={"limit": 20},
                timeout=context.timeout_seconds,
            )
        molecules = payload.get("molecules") or []
        if not molecules:
            return ""
        return str(molecules[0].get("molecule_chembl_id") or "")

    def _target_symbol(
        self,
        target_id: str,
        context: EvidenceContext,
        cache: dict[str, str],
    ) -> str:
        if target_id in cache:
            return cache[target_id]
        payload = _json_get(
            f"{self.BASE}/target/{target_id}.json",
            timeout=context.timeout_seconds,
        )
        if str(payload.get("organism") or "") != "Homo sapiens":
            cache[target_id] = ""
            return ""
        symbol = ""
        for component in payload.get("target_components") or []:
            for synonym in component.get("target_component_synonyms") or []:
                if str(synonym.get("syn_type") or "").upper().startswith(
                    "GENE_SYMBOL"
                ):
                    symbol = _normalise_symbol(synonym.get("component_synonym"))
                    if symbol:
                        break
            if symbol:
                break
        cache[target_id] = symbol
        return symbol

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        molecule_id = self._resolve_molecule_id(context)
        if not molecule_id:
            return []
        payload = _json_get(
            f"{self.BASE}/activity.json",
            params={
                "molecule_chembl_id": molecule_id,
                "limit": min(1000, context.max_records_per_source),
            },
            timeout=context.timeout_seconds,
        )
        activities = payload.get("activities") or []
        records: list[EvidenceRecord] = []
        target_cache: dict[str, str] = {}
        for activity in activities:
            target_id = str(activity.get("target_chembl_id") or "")
            activity_id = str(activity.get("activity_id") or "")
            if not target_id or not activity_id:
                continue
            symbol = self._target_symbol(target_id, context, target_cache)
            if not symbol:
                continue
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=activity_id,
                    evidence_type="direct_bioactivity",
                    subject_type="compound",
                    subject_id=molecule_id,
                    relation="targets",
                    object_type="target",
                    object_id=target_id,
                    target_symbol=symbol,
                    tier=self.tier,
                    source_version=self.source_version,
                    source_group=self.source_group,
                    score=_score_from_pchembl(activity.get("pchembl_value")),
                    effect_size=_float_or_none(activity.get("standard_value")),
                    assay=str(activity.get("assay_description") or ""),
                    direction=str(activity.get("action_type") or "unknown"),
                    species="Homo sapiens",
                    url=f"{self.BASE}/activity/{activity_id}",
                    license="CC BY-SA 3.0",
                    payload={
                        "standard_type": activity.get("standard_type"),
                        "standard_units": activity.get("standard_units"),
                        "target_pref_name": activity.get("target_pref_name"),
                    },
                )
            )
            if len(records) >= context.max_records_per_source:
                break
        return records


class BindingDBConnector(BaseConnector):
    name = "BindingDB"
    source_version = "BindingDB REST API"
    source_group = "bindingdb"
    tier = EvidenceTier.EXPERIMENTAL

    BASE = "https://bindingdb.org"

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        smiles = str(
            context.compound.get("isomeric_smiles")
            or context.compound.get("canonical_smiles")
            or ""
        ).strip()
        if not smiles:
            return []
        payload = _json_get(
            f"{self.BASE}/rest/getTargetsByCompound",
            params={
                "smiles": smiles,
                "cutoff": 0.85,
                "identity": 85,
                "response": "application/json",
            },
            timeout=context.timeout_seconds,
        )
        records: list[EvidenceRecord] = []
        for index, row in enumerate(_recursive_dicts(payload)):
            target_name = _first_key(
                row,
                [
                    "Target Name",
                    "target_name",
                    "protein_name",
                    "UniProt (SwissProt) Recommended Name",
                ],
            )
            gene_name = _first_key(
                row,
                ["Gene Name", "gene_name", "gene", "symbol"],
            )
            uniprot = _first_key(
                row,
                [
                    "UniProt (SwissProt) Primary ID",
                    "uniprot",
                    "uniprot_id",
                    "Target Accession",
                ],
            )
            symbol = _normalise_symbol(gene_name) or _normalise_symbol(target_name)
            if not symbol:
                continue
            affinity = _first_key(
                row,
                [
                    "Ki (nM)",
                    "IC50 (nM)",
                    "Kd (nM)",
                    "EC50 (nM)",
                    "affinity",
                    "value",
                ],
            )
            affinity_nm = _float_or_none(affinity)
            score = None
            if affinity_nm is not None and affinity_nm > 0:
                score = float(np.clip((9.0 - math.log10(affinity_nm)) / 5.0, 0.0, 1.0))
            object_id = str(uniprot or target_name or symbol)
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=f"{context.compound.get('pubchem_cid', 'compound')}:{index}:{symbol}",
                    evidence_type="direct_binding",
                    subject_type="compound",
                    subject_id=str(
                        context.compound.get("pubchem_cid")
                        or context.compound.get("inchi_key")
                        or "compound"
                    ),
                    relation="binds",
                    object_type="target",
                    object_id=object_id,
                    target_symbol=symbol,
                    tier=self.tier,
                    source_version=self.source_version,
                    source_group=self.source_group,
                    score=score,
                    effect_size=affinity_nm,
                    assay=str(
                        _first_key(row, ["Assay Description", "assay", "assay_type"])
                        or ""
                    ),
                    species=str(_first_key(row, ["Species", "organism"]) or ""),
                    url=self.BASE,
                    license="BindingDB terms",
                    payload={"raw_affinity": affinity, "target_name": target_name},
                )
            )
            if len(records) >= context.max_records_per_source:
                break
        return records


class PubChemBioAssayConnector(BaseConnector):
    name = "PubChemBioAssay"
    source_version = "PubChem PUG REST"
    source_group = "pubchem"
    tier = EvidenceTier.EXPERIMENTAL

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        cid = str(context.compound.get("pubchem_cid") or "").strip()
        if not cid:
            return []
        payload = _json_get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON",
            timeout=context.timeout_seconds,
        )
        table = payload.get("Table") or {}
        columns = [str(value) for value in table.get("Columns") or []]
        rows = table.get("Row") or table.get("Rows") or []
        if not columns or not rows:
            return []
        index = {name: number for number, name in enumerate(columns)}

        def value(row: list[Any], *aliases: str) -> Any:
            for alias in aliases:
                for name, number in index.items():
                    if alias.lower() in name.lower():
                        return row[number] if number < len(row) else None
            return None

        records: list[EvidenceRecord] = []
        for row_number, row in enumerate(rows):
            aid = value(row, "AID")
            outcome = value(row, "Activity Outcome", "ActivityOutcome")
            target_name = value(row, "Target Name", "TargetName")
            target_accession = value(row, "Target Accession", "TargetAccession")
            symbol = _normalise_symbol(target_name)
            if not symbol:
                continue
            outcome_text = str(outcome or "").lower()
            score = 1.0 if "active" in outcome_text else 0.5
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=f"{cid}:{aid or row_number}",
                    evidence_type="bioassay_activity",
                    subject_type="compound",
                    subject_id=cid,
                    relation="targets",
                    object_type="target",
                    object_id=str(target_accession or target_name or symbol),
                    target_symbol=symbol,
                    tier=self.tier,
                    source_version=self.source_version,
                    source_group=self.source_group,
                    score=score,
                    assay=str(value(row, "Assay Name", "AssayName") or ""),
                    species="Homo sapiens",
                    url=f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{aid}" if aid else "",
                    license="PubChem terms",
                    payload={"activity_outcome": outcome},
                )
            )
            if len(records) >= context.max_records_per_source:
                break
        return records


class GWASCatalogConnector(BaseConnector):
    name = "GWASCatalog"
    source_version = "GWAS Catalog REST API v2"
    source_group = "gwas_catalog"
    tier = EvidenceTier.GENETIC

    BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"

    def _gene_symbols(self, record: Mapping[str, Any]) -> set[str]:
        symbols: set[str] = set()
        for nested in _recursive_dicts(record):
            for key, value in nested.items():
                key_text = re.sub(r"[^a-z0-9]+", "", str(key).lower())
                if "gene" not in key_text and "symbol" not in key_text:
                    continue
                values = value if isinstance(value, list) else [value]
                for item in values:
                    if isinstance(item, dict):
                        candidate = (
                            item.get("geneName")
                            or item.get("gene_name")
                            or item.get("symbol")
                            or item.get("name")
                        )
                    else:
                        candidate = item
                    symbol = _normalise_symbol(candidate)
                    if symbol:
                        symbols.add(symbol)
        return symbols

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        disease_name = str(context.disease.get("name") or "").strip()
        if not disease_name:
            return []
        studies = _json_get(
            f"{self.BASE}/studies",
            params={"efo_trait": disease_name, "size": min(100, context.max_records_per_source)},
            timeout=context.timeout_seconds,
        ).get("_embedded", {}).get("studies", [])
        records: list[EvidenceRecord] = []
        for study in studies:
            accession = str(study.get("accessionId") or study.get("accession") or "")
            if not accession:
                continue
            payload = _json_get(
                f"{self.BASE}/studies/{urllib.parse.quote(accession)}/associations",
                timeout=context.timeout_seconds,
            )
            associations = (
                payload.get("_embedded", {}).get("associations")
                or payload.get("associations")
                or []
            )
            for association in associations:
                symbols = self._gene_symbols(association)
                if not symbols:
                    continue
                p_value = _first_key(
                    association,
                    ["pvalue", "p_value", "pValue"],
                )
                if p_value is None:
                    for nested in _recursive_dicts(association):
                        p_value = _first_key(
                            nested,
                            ["pvalue", "p_value", "pValue"],
                        )
                        if p_value is not None:
                            break
                association_id = str(
                    association.get("associationId")
                    or association.get("accessionId")
                    or f"{accession}:{len(records)}"
                )
                for symbol in sorted(symbols):
                    records.append(
                        EvidenceRecord(
                            source=self.name,
                            source_record_id=f"{association_id}:{symbol}",
                            evidence_type="genetic_association",
                            subject_type="disease",
                            subject_id=disease_name,
                            relation="associated_with",
                            object_type="target",
                            object_id=symbol,
                            target_symbol=symbol,
                            tier=self.tier,
                            source_version=self.source_version,
                            source_group=self.source_group,
                            score=_score_from_pvalue(p_value),
                            p_value=_float_or_none(p_value),
                            species="Homo sapiens",
                            url=f"https://www.ebi.ac.uk/gwas/studies/{accession}",
                            license="GWAS Catalog terms",
                            payload={"study_accession": accession},
                        )
                    )
                    if len(records) >= context.max_records_per_source:
                        return records
        return records


class ClinVarConnector(BaseConnector):
    name = "ClinVar"
    source_version = "NCBI E-utilities/ClinVar"
    source_group = "clinvar"
    tier = EvidenceTier.GENETIC

    BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        options = dict(context.source_options.get(self.name) or {})
        limit = int(options.get("max_targets", 100) or 100)
        records: list[EvidenceRecord] = []
        for symbol in context.target_symbols[:limit]:
            gene = _normalise_symbol(symbol)
            if not gene:
                continue
            payload = _json_get(
                self.BASE,
                params={
                    "db": "clinvar",
                    "term": f"{gene}[gene] AND human[orgn]",
                    "retmode": "json",
                    "retmax": 1,
                },
                timeout=context.timeout_seconds,
            )
            count = int(
                ((payload.get("esearchresult") or {}).get("count") or 0)
            )
            if count <= 0:
                continue
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=f"{gene}:clinvar_count",
                    evidence_type="clinical_variant_association",
                    subject_type="disease",
                    subject_id=str(
                        context.disease.get("id")
                        or context.disease.get("name")
                        or "disease"
                    ),
                    relation="associated_with",
                    object_type="target",
                    object_id=gene,
                    target_symbol=gene,
                    tier=self.tier,
                    source_version=self.source_version,
                    source_group=self.source_group,
                    score=float(
                        np.clip(
                            math.log1p(count) / math.log1p(100.0),
                            0.0,
                            1.0,
                        )
                    ),
                    sample_size=count,
                    species="Homo sapiens",
                    url=(
                        "https://www.ncbi.nlm.nih.gov/clinvar/?term="
                        + urllib.parse.quote(f"{gene}[gene]")
                    ),
                    license="NCBI public data terms",
                    payload={"clinvar_record_count": count},
                )
            )
            if len(records) >= context.max_records_per_source:
                break
        return records


class ClinicalTrialsGovConnector(BaseConnector):
    name = "ClinicalTrialsGov"
    source_version = "ClinicalTrials.gov API v2"
    source_group = "clinicaltrials"
    tier = EvidenceTier.CURATED

    BASE = "https://clinicaltrials.gov/api/v2/studies"

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        disease = str(
            context.disease.get("name")
            or context.disease.get("id")
            or ""
        ).strip()
        if not disease:
            return []
        options = dict(context.source_options.get(self.name) or {})
        limit = int(options.get("max_targets", 50) or 50)
        records: list[EvidenceRecord] = []
        for symbol in context.target_symbols[:limit]:
            gene = _normalise_symbol(symbol)
            if not gene:
                continue
            payload = _json_get(
                self.BASE,
                params={
                    "query.term": f"{gene} {disease}",
                    "pageSize": 1,
                    "countTotal": "true",
                },
                timeout=context.timeout_seconds,
            )
            count = int(payload.get("totalCount") or 0)
            if count <= 0:
                continue
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=f"{gene}:{disease}:clinical_trials",
                    evidence_type="clinical_precedent",
                    subject_type="disease",
                    subject_id=disease,
                    relation="studied_in",
                    object_type="target",
                    object_id=gene,
                    target_symbol=gene,
                    tier=self.tier,
                    source_version=self.source_version,
                    source_group=self.source_group,
                    score=float(
                        np.clip(
                            math.log1p(count) / math.log1p(50.0),
                            0.0,
                            1.0,
                        )
                    ),
                    sample_size=count,
                    species="Homo sapiens",
                    url=(
                        "https://clinicaltrials.gov/search?term="
                        + urllib.parse.quote(f"{gene} {disease}")
                    ),
                    license="ClinicalTrials.gov public data terms",
                    payload={"clinical_trial_count": count},
                )
            )
            if len(records) >= context.max_records_per_source:
                break
        return records


class GTExConnector(BaseConnector):
    name = "GTEx"
    source_version = "GTEx API v2"
    source_group = "gtex"
    tier = EvidenceTier.CONTEXT

    BASE = "https://gtexportal.org/api/v2"

    def _ensembl_id(self, symbol: str, context: EvidenceContext) -> str:
        known = context.ensembl_ids.get(symbol.upper())
        if known:
            return known
        payload = _json_get(
            "https://mygene.info/v3/query",
            params={
                "q": symbol,
                "species": "human",
                "fields": "ensembl.gene",
                "size": 1,
            },
            timeout=context.timeout_seconds,
        )
        hits = payload if isinstance(payload, list) else [payload]
        for hit in hits:
            ensembl = (hit or {}).get("ensembl")
            if isinstance(ensembl, list):
                ensembl = ensembl[0] if ensembl else {}
            value = (ensembl or {}).get("gene")
            if value:
                return str(value)
        return ""

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        symbols = sorted(
            {
                _normalise_symbol(value)
                for value in context.target_symbols
                if _normalise_symbol(value)
            }
        )
        records: list[EvidenceRecord] = []
        for symbol in symbols:
            ensembl = self._ensembl_id(symbol, context)
            if not ensembl:
                continue
            payload = _json_get(
                f"{self.BASE}/expression/medianGeneExpression",
                params={
                    "gencodeId": ensembl,
                    "datasetId": "gtex_v8",
                    "itemsPerPage": 1000,
                },
                timeout=context.timeout_seconds,
            )
            rows = payload.get("medianGeneExpression") or payload.get("data") or []
            for row in rows:
                tissue = str(
                    row.get("tissueSiteDetailId")
                    or row.get("tissueSiteDetail")
                    or row.get("tissue")
                    or ""
                )
                if "liver" not in tissue.lower():
                    continue
                value = _float_or_none(
                    row.get("median")
                    or row.get("medianExpression")
                    or row.get("value")
                    or row.get("tpm")
                )
                if value is None:
                    continue
                records.append(
                    EvidenceRecord(
                        source=self.name,
                        source_record_id=f"{ensembl}:{tissue}",
                        evidence_type="expression_context",
                        subject_type="target",
                        subject_id=ensembl,
                        relation="expressed_in",
                        object_type="tissue",
                        object_id="liver",
                        target_symbol=symbol,
                        tier=self.tier,
                        source_version=self.source_version,
                        source_group=self.source_group,
                        score=float(
                            np.clip(math.log1p(max(value, 0.0)) / math.log1p(1000.0), 0.0, 1.0)
                        ),
                        effect_size=value,
                        tissue="liver",
                        species="Homo sapiens",
                        url="https://gtexportal.org/home/tissue/liver",
                        license="GTEx terms",
                        payload={"tissue_id": tissue},
                    )
                )
                break
            if len(records) >= context.max_records_per_source:
                break
        return records


class HPAConnector(BaseConnector):
    name = "HumanProteinAtlas"
    source_version = "HPA API"
    source_group = "hpa"
    tier = EvidenceTier.CONTEXT

    BASE = "https://www.proteinatlas.org"

    def _liver_expression(self, payload: Mapping[str, Any]) -> float | None:
        candidates: list[tuple[str, float]] = []
        for record in _recursive_dicts(payload):
            tissue = " ".join(
                str(value)
                for key, value in record.items()
                if "tissue" in str(key).lower() or str(key).lower() in {"name", "label"}
            )
            if "liver" not in tissue.lower():
                continue
            for key, value in record.items():
                key_text = str(key).lower()
                if not any(token in key_text for token in ("ntpm", "tpm", "value", "expression")):
                    continue
                number = _float_or_none(value)
                if number is not None:
                    candidates.append((tissue, number))
        if not candidates:
            return None
        return max(value for _, value in candidates)

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        if not context.allow_network:
            raise ConnectorError("network disabled")
        records: list[EvidenceRecord] = []
        for symbol, ensembl in sorted(context.ensembl_ids.items()):
            if not ensembl:
                continue
            payload = _json_get(
                f"{self.BASE}/{urllib.parse.quote(str(ensembl))}.json",
                timeout=context.timeout_seconds,
            )
            value = self._liver_expression(payload)
            if value is None:
                continue
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=f"{ensembl}:liver",
                    evidence_type="protein_expression_context",
                    subject_type="target",
                    subject_id=str(ensembl),
                    relation="expressed_in",
                    object_type="tissue",
                    object_id="liver",
                    target_symbol=_normalise_symbol(symbol),
                    tier=self.tier,
                    source_version=self.source_version,
                    source_group=self.source_group,
                    score=float(
                        np.clip(math.log1p(max(value, 0.0)) / math.log1p(1000.0), 0.0, 1.0)
                    ),
                    effect_size=value,
                    tissue="liver",
                    species="Homo sapiens",
                    url=f"{self.BASE}/{ensembl}-LIVER",
                    license="HPA terms",
                    payload={},
                )
            )
            if len(records) >= context.max_records_per_source:
                break
        return records


@dataclass(slots=True)
class DepMapLocalConnector(BaseConnector):
    path: str
    name: str = "DepMap"
    source_version: str = "local snapshot"
    source_group: str = "depmap"
    tier: EvidenceTier = EvidenceTier.CONTEXT
    disease_filter: str = "liver"
    requires_network = False

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        path = Path(self.path).expanduser()
        if not path.is_absolute():
            base = context.cache_dir or Path.cwd()
            path = (base / path).resolve()
        if not path.exists():
            raise ConnectorError(f"DepMap file not found: {path}")
        frame = pd.read_csv(path, low_memory=False)
        normalised = {
            re.sub(r"[^a-z0-9]+", "", str(column).lower()): str(column)
            for column in frame.columns
        }
        valid_symbols = {
            _normalise_symbol(value)
            for value in context.target_symbols
            if _normalise_symbol(value)
        }
        matrix_columns: dict[str, str] = {}
        for symbol in valid_symbols:
            exact = next(
                (
                    column
                    for column in frame.columns
                    if str(column).split(" ")[0].upper() == symbol
                ),
                None,
            )
            if exact:
                matrix_columns[symbol] = exact
        if matrix_columns:
            records: list[EvidenceRecord] = []
            for symbol, column in sorted(matrix_columns.items()):
                values = pd.to_numeric(frame[column], errors="coerce").dropna()
                if values.empty:
                    continue
                effect = float(values.median())
                records.append(
                    EvidenceRecord(
                        source=self.name,
                        source_record_id=f"{symbol}:{self.disease_filter}",
                        evidence_type="dependency",
                        subject_type="target",
                        subject_id=symbol,
                        relation="dependent_on",
                        object_type="disease_context",
                        object_id=self.disease_filter,
                        target_symbol=symbol,
                        tier=self.tier,
                        source_version=self.source_version,
                        source_group=self.source_group,
                        score=float(np.clip(-effect / 2.0, 0.0, 1.0)),
                        effect_size=effect,
                        species="Homo sapiens",
                        tissue=self.disease_filter,
                        url="https://depmap.org/portal/",
                        license="CC BY 4.0",
                        payload={
                            "aggregation": "median",
                            "file": str(path),
                            "column": str(column),
                        },
                    )
                )
                if len(records) >= context.max_records_per_source:
                    break
            return records
        gene_column = next(
            (
                normalised[key]
                for key in ("genename", "gene", "symbol", "genesyMBOL")
                if key in normalised
            ),
            None,
        )
        effect_column = next(
            (
                normalised[key]
                for key in ("geneeffect", "chronos", "effect", "dependency")
                if key in normalised
            ),
            None,
        )
        lineage_column = next(
            (
                normalised.get(key)
                for key in ("lineage", "primarydisease", "disease")
                if normalised.get(key)
            ),
            None,
        )
        cell_line_column = normalised.get("depmapid") or normalised.get("cellline")
        records: list[EvidenceRecord] = []
        if gene_column and effect_column:
            keep_columns = [gene_column, effect_column]
            if lineage_column:
                keep_columns.append(lineage_column)
            work = frame[keep_columns].copy()
            if self.disease_filter and lineage_column:
                work = work[
                    work[lineage_column].astype(str).str.contains(
                        self.disease_filter,
                        case=False,
                        na=False,
                    )
                ]
            work[effect_column] = pd.to_numeric(work[effect_column], errors="coerce")
            aggregated = (
                work.dropna()
                .groupby(gene_column, as_index=False)[effect_column]
                .median()
            )
            for row in aggregated.itertuples(index=False):
                symbol = _normalise_symbol(getattr(row, gene_column))
                effect = _float_or_none(getattr(row, effect_column))
                if not symbol or effect is None:
                    continue
                records.append(
                    EvidenceRecord(
                        source=self.name,
                        source_record_id=f"{symbol}:{self.disease_filter}",
                        evidence_type="dependency",
                        subject_type="target",
                        subject_id=symbol,
                        relation="dependent_on",
                        object_type="disease_context",
                        object_id=self.disease_filter,
                        target_symbol=symbol,
                        tier=self.tier,
                        source_version=self.source_version,
                        source_group=self.source_group,
                        score=float(np.clip(-effect / 2.0, 0.0, 1.0)),
                        effect_size=effect,
                        species="Homo sapiens",
                        tissue=self.disease_filter,
                        url="https://depmap.org/portal/",
                        license="CC BY 4.0",
                        payload={"aggregation": "median", "file": str(path)},
                    )
                )
                if len(records) >= context.max_records_per_source:
                    break
        elif cell_line_column:
            work = frame
            if self.disease_filter and lineage_column:
                mask = work[lineage_column].astype(str).str.contains(
                    self.disease_filter,
                    case=False,
                    na=False,
                )
                work = work[mask]
            if gene_column and effect_column:
                for row in work.itertuples(index=False):
                    symbol = _normalise_symbol(getattr(row, gene_column))
                    effect = _float_or_none(getattr(row, effect_column))
                    if not symbol or effect is None:
                        continue
                    cell_line = str(getattr(row, cell_line_column))
                    records.append(
                        EvidenceRecord(
                            source=self.name,
                            source_record_id=f"{symbol}:{cell_line}",
                            evidence_type="dependency",
                            subject_type="target",
                            subject_id=symbol,
                            relation="dependent_on",
                            object_type="cell_line",
                            object_id=cell_line,
                            target_symbol=symbol,
                            tier=self.tier,
                            source_version=self.source_version,
                            source_group=self.source_group,
                            score=float(np.clip(-effect / 2.0, 0.0, 1.0)),
                            effect_size=effect,
                            species="Homo sapiens",
                            tissue=self.disease_filter,
                            url="https://depmap.org/portal/",
                            license="CC BY 4.0",
                            payload={"file": str(path)},
                        )
                    )
                    if len(records) >= context.max_records_per_source:
                        break
        return records


class LocalTableConnector(BaseConnector):
    """Import a user-supplied table (CTD, Tox21, LINCS, DisGeNET, etc.)."""

    def __init__(
        self,
        *,
        name: str,
        path: str,
        source_version: str,
        evidence_type: str,
        relation: str,
        tier: EvidenceTier,
        subject_type: str,
        subject_column: str | None = None,
        subject_value: str | None = None,
        object_type: str = "target",
        object_column: str | None = None,
        target_column: str = "gene",
        score_column: str | None = None,
        effect_column: str | None = None,
        p_value_column: str | None = None,
        direction_column: str | None = None,
        source_group: str | None = None,
    ) -> None:
        self.name = name
        self.path = path
        self.source_version = source_version
        self.evidence_type = evidence_type
        self.relation = relation
        self.tier = tier
        self.subject_type = subject_type
        self.subject_column = subject_column
        self.subject_value = subject_value
        self.object_type = object_type
        self.object_column = object_column
        self.target_column = target_column
        self.score_column = score_column
        self.effect_column = effect_column
        self.p_value_column = p_value_column
        self.direction_column = direction_column
        self.source_group = source_group or name.lower()
        self.requires_network = False

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        path = Path(self.path).expanduser()
        if not path.is_absolute():
            base = context.cache_dir or Path.cwd()
            path = (base / path).resolve()
        if not path.exists():
            raise ConnectorError(f"{self.name} table not found: {path}")
        separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        frame = pd.read_csv(path, sep=separator, dtype=str, low_memory=False)
        if self.target_column not in frame.columns:
            raise ConnectorError(
                f"{self.name} table lacks target column {self.target_column!r}: "
                f"{list(frame.columns)}"
            )
        records: list[EvidenceRecord] = []
        for index, row in frame.iterrows():
            symbol = _normalise_symbol(row.get(self.target_column))
            if not symbol:
                continue
            subject_id = (
                str(row.get(self.subject_column) or "")
                if self.subject_column
                else str(self.subject_value or self.subject_type)
            )
            object_id = (
                str(row.get(self.object_column) or symbol)
                if self.object_column
                else symbol
            )
            score = _float_or_none(row.get(self.score_column)) if self.score_column else None
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=str(
                        row.get("id")
                        or row.get("record_id")
                        or f"{index}:{symbol}"
                    ),
                    evidence_type=self.evidence_type,
                    subject_type=self.subject_type,
                    subject_id=subject_id,
                    relation=self.relation,
                    object_type=self.object_type,
                    object_id=object_id,
                    target_symbol=symbol,
                    tier=self.tier,
                    source_version=self.source_version,
                    source_group=self.source_group,
                    score=score,
                    effect_size=(
                        _float_or_none(row.get(self.effect_column))
                        if self.effect_column
                        else None
                    ),
                    p_value=(
                        _float_or_none(row.get(self.p_value_column))
                        if self.p_value_column
                        else None
                    ),
                    direction=(
                        str(row.get(self.direction_column) or "unknown")
                        if self.direction_column
                        else "unknown"
                    ),
                    url="",
                    license=str(row.get("license") or ""),
                    payload={"file": str(path), "row": int(index)},
                )
            )
            if len(records) >= context.max_records_per_source:
                break
        return records


def connector_from_config(
    name: str,
    options: Mapping[str, Any],
) -> BaseConnector:
    """Instantiate one configured connector without importing hidden state."""
    key = name.strip().lower()
    local_names = {
        "ctd",
        "tox21",
        "toxcast",
        "comptox",
        "lincs",
        "disgenet",
        "drugbank",
        "iuphar",
        "genecards",
        "omim",
        "ttd",
        "local",
    }
    if key == "opentargets":
        return OpenTargetsConnector()
    if key == "chembl":
        return ChEMBLConnector()
    if key == "bindingdb":
        return BindingDBConnector()
    if key in {"pubchem", "pubchembioassay"}:
        return PubChemBioAssayConnector()
    if key in {"gwas", "gwascatalog"}:
        return GWASCatalogConnector()
    if key == "clinvar":
        return ClinVarConnector()
    if key in {"clinicaltrials", "clinicaltrialsgov"}:
        return ClinicalTrialsGovConnector()
    if key == "gtex":
        return GTExConnector()
    if key in {"hpa", "humanproteinatlas"}:
        return HPAConnector()
    if key == "depmap":
        path = str(options.get("path") or "").strip()
        if not path:
            raise ConnectorError("DepMap connector requires a local path")
        return DepMapLocalConnector(
            path=path,
            disease_filter=str(options.get("disease_filter") or "liver"),
        )
    if key in local_names or str(options.get("connector") or "").lower() == "local":
        required = ["path", "evidence_type", "relation", "tier"]
        missing = [field for field in required if not options.get(field)]
        if missing:
            raise ConnectorError(
                f"local connector missing options: {', '.join(missing)}"
            )
        return LocalTableConnector(
            name=str(options.get("name") or name),
            path=str(options["path"]),
            source_version=str(options.get("source_version") or "local snapshot"),
            evidence_type=str(options["evidence_type"]),
            relation=str(options["relation"]),
            tier=EvidenceTier.parse(options["tier"]),
            subject_type=str(options.get("subject_type") or "compound"),
            subject_column=options.get("subject_column"),
            subject_value=options.get("subject_value"),
            object_type=str(options.get("object_type") or "target"),
            object_column=options.get("object_column"),
            target_column=str(options.get("target_column") or "gene"),
            score_column=options.get("score_column"),
            effect_column=options.get("effect_column"),
            p_value_column=options.get("p_value_column"),
            direction_column=options.get("direction_column"),
            source_group=options.get("source_group"),
        )
    raise ConnectorError(f"unknown evidence connector: {name}")
