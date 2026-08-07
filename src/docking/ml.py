#!/usr/bin/env python3
"""Machine learning and deep learning rescoring for docking results."""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import (
    GradientBoostingClassifier,
    GradientBoostingRegressor,
    RandomForestClassifier,
    RandomForestRegressor,
)
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    mean_squared_error,
    r2_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.preprocessing import LabelEncoder, StandardScaler

from .config import ResolvedConfig
from .utils import DockingError, write_json

DESCRIPTOR_NAMES = [
    "MW",
    "LogP",
    "HBD",
    "HBA",
    "TPSA",
    "RotBonds",
    "AromaticRings",
    "HeavyAtoms",
]


def _features(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    desc = [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcTPSA(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
        mol.GetNumHeavyAtoms(),
    ]
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    return desc + list(fp)


def feature_names() -> list[str]:
    return DESCRIPTOR_NAMES + [f"morgan_{i}" for i in range(2048)]


def build_feature_matrix(df: pd.DataFrame) -> tuple[np.ndarray, pd.Index]:
    valid_idx: list[int] = []
    rows: list[list[float]] = []
    for i, smiles in enumerate(df["smiles"].astype(str)):
        feats = _features(smiles)
        if feats is None:
            continue
        valid_idx.append(i)
        rows.append(feats)
    if not rows:
        raise DockingError("no valid SMILES found for ML feature generation")
    return np.asarray(rows, dtype="float64"), pd.Index(valid_idx)


def train_ml(
    cfg: ResolvedConfig,
    log,
    model_type: str | None = None,
    label_column: str | None = None,
    training_csv: str | None = None,
) -> dict:
    model_type = (model_type or cfg.get("ml", "model", "rf")).lower()
    df, y, task = _load_training(cfg, label_column, training_csv)
    X, idx = build_feature_matrix(df)
    y = y.iloc[idx].to_numpy()
    if task == "classification":
        encoder = LabelEncoder()
        y = encoder.fit_transform(y)

    test_size = float(cfg.get("ml", "test_size", 0.2))
    stratify = y if task == "classification" and len(set(y)) > 1 else None
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=test_size,
        random_state=int(cfg.get("ml", "random_state", 42)),
        stratify=stratify,
    )

    model = _build_model(
        model_type,
        task,
        int(cfg.get("ml", "hidden_size", 128)),
        int(cfg.get("ml", "epochs", 80)),
        int(cfg.get("ml", "random_state", 42)),
    )
    model.fit(X_train, y_train)

    reports = cfg.reports_dir()
    reports.mkdir(parents=True, exist_ok=True)
    is_torch = model_type == "torch"
    model_file = reports / ("ml_model.pt" if is_torch else "ml_model.joblib")
    if is_torch:
        import torch

        torch.save(model, model_file)
    else:
        joblib.dump(model, model_file)

    metrics = _evaluate(model, X_test, y_test, task, encoder if task == "classification" else None)
    info = {
        "model_type": model_type,
        "task": task,
        "model_file": str(model_file),
        "feature_count": X.shape[1],
        "train_samples": int(len(y_train)),
        "test_samples": int(len(y_test)),
        "metrics": metrics,
    }
    if task == "classification":
        info["classes"] = encoder.classes_.tolist()
    write_json(reports / "ml_model_info.json", info)

    if hasattr(model, "feature_importances_"):
        _save_importance(model, X.shape[1], reports)
    if task == "classification" and hasattr(model, "predict_proba"):
        _save_roc(model, X_test, y_test, reports)

    log.info(
        "ML training complete: %s/%s, %s",
        model_type,
        task,
        metrics,
    )
    return info


def predict_ml(cfg: ResolvedConfig, log) -> dict:
    reports = cfg.reports_dir()
    info_path = reports / "ml_model_info.json"
    if not info_path.exists():
        raise DockingError("no trained ML model found; run ml-train first")
    info = _read_json(info_path)
    model_file = Path(info["model_file"])
    if not model_file.exists():
        raise DockingError(f"model file not found: {model_file}")
    if info["model_type"] == "torch":
        import torch

        model = torch.load(model_file, weights_only=False)
    else:
        model = joblib.load(model_file)
    task = info["task"]

    results_csv = reports / "ranked_results.csv"
    if not results_csv.exists():
        results_csv = cfg.results_path()
    if not results_csv.exists():
        raise DockingError("docking results not found; run dock/analyze first")
    df = pd.read_csv(results_csv, dtype={"id": str})
    if "smiles" not in df.columns:
        raise DockingError("results CSV has no smiles column")
    X, idx = build_feature_matrix(df)
    pred = model.predict(X) if not hasattr(model, "predict_proba") else model.predict_proba(X)
    if task == "classification" and pred.ndim == 2:
        ml_score = pred[:, 1]
    else:
        ml_score = np.asarray(pred, dtype="float64").ravel()

    out = df.iloc[idx].copy()
    out["ml_score"] = ml_score
    ascending = task == "regression"
    out["ml_rank"] = out["ml_score"].rank(ascending=ascending).astype(int)
    if "affinity" in out.columns and out["affinity"].notna().any():
        affinity_rank = out["affinity"].rank(ascending=True).astype(int)
        out["combined_rank"] = ((affinity_rank + out["ml_rank"]) / 2).astype(float)
    else:
        out["combined_rank"] = out["ml_rank"]
    out = out.sort_values("ml_rank")
    out_path = reports / "ml_ranked_results.csv"
    out.to_csv(out_path, index=False)
    summary = {
        "scored": int(len(out)),
        "task": task,
        "output": str(out_path),
    }
    write_json(reports / "ml_predict_summary.json", summary)
    log.info("ML prediction complete: %s ligands scored -> %s", len(out), out_path)
    return summary


def _load_training(cfg, label_column, training_csv):
    path = Path(training_csv) if training_csv else cfg.ml_training_csv()
    if not path.exists():
        raise DockingError(f"training CSV not found: {path}")
    df = pd.read_csv(path)
    if "smiles" not in df.columns:
        raise DockingError("training CSV must contain a smiles column")
    label = label_column or cfg.get("ml", "label_column", "active")
    task = str(cfg.get("ml", "task", "auto")).lower()
    if label not in df.columns and "affinity" in df.columns:
        label = "affinity"
    if label not in df.columns:
        raise DockingError(
            f"label column '{label}' not found in {path.name}"
        )
    y = df[label]
    if task == "auto":
        if pd.api.types.is_numeric_dtype(y):
            task = "classification" if y.nunique() <= 2 else "regression"
        else:
            task = "classification"
    if task not in ("classification", "regression"):
        raise DockingError(f"unsupported ML task: {task}")
    return df, y, task


def _build_model(model_type, task, hidden, epochs, random_state):
    if model_type == "rf":
        if task == "classification":
            return RandomForestClassifier(n_estimators=300, random_state=random_state, n_jobs=-1)
        return RandomForestRegressor(n_estimators=300, random_state=random_state, n_jobs=-1)
    if model_type == "gbm":
        if task == "classification":
            return GradientBoostingClassifier(random_state=random_state)
        return GradientBoostingRegressor(random_state=random_state)
    if model_type == "torch":
        try:
            import torch  # noqa: F401

            return _TorchMLP(
                hidden=hidden,
                task=task,
                epochs=epochs,
                random_state=random_state,
            )
        except Exception:
            model_type = "mlp"
    if task == "classification":
        return MLPClassifier(
            hidden_layer_sizes=(hidden, hidden),
            max_iter=epochs,
            early_stopping=True,
            random_state=random_state,
        )
    return MLPRegressor(
        hidden_layer_sizes=(hidden, hidden),
        max_iter=epochs,
        early_stopping=True,
        random_state=random_state,
    )


class _TorchMLP:
    def __init__(self, hidden=128, task="classification", epochs=80, random_state=42):
        import torch
        from torch import nn

        self.torch = torch
        self.task = task
        self.epochs = epochs
        self.random_state = random_state
        self.scaler = StandardScaler()
        self.net = None
        self.n_features = None

    def fit(self, X, y):
        import torch
        from torch import nn
        from torch.utils.data import DataLoader, TensorDataset

        torch.manual_seed(self.random_state)
        self.n_features = X.shape[1]
        out_units = 2 if self.task == "classification" else 1
        self.net = nn.Sequential(
            nn.Linear(self.n_features, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, out_units),
        )
        Xs = self.scaler.fit_transform(X).astype("float32")
        if self.task == "classification":
            ys = y.astype("int64")
        else:
            ys = y.astype("float32").reshape(-1, 1)
        dataset = TensorDataset(torch.tensor(Xs), torch.tensor(ys))
        loader = DataLoader(
            dataset,
            batch_size=min(32, len(dataset)),
            shuffle=True,
        )
        opt = torch.optim.Adam(self.net.parameters(), lr=1e-3)
        loss_fn = (
            nn.CrossEntropyLoss()
            if self.task == "classification"
            else nn.MSELoss()
        )
        self.net.train()
        for _ in range(self.epochs):
            for xb, yb in loader:
                opt.zero_grad()
                loss = loss_fn(self.net(xb), yb)
                loss.backward()
                opt.step()
        self.net.eval()
        return self

    def predict_proba(self, X):
        Xs = self.scaler.transform(X).astype("float32")
        with self.torch.no_grad():
            out = self.net(self.torch.tensor(Xs))
            if self.task == "classification":
                return self.torch.softmax(out, dim=1).numpy()
            return out.numpy().ravel()

    def predict(self, X):
        proba = self.predict_proba(X)
        if self.task == "classification":
            return proba.argmax(axis=1)
        return proba


def _evaluate(model, X_test, y_test, task, encoder=None):
    if task == "classification":
        pred = model.predict(X_test)
        acc = accuracy_score(y_test, pred)
        metrics = {"accuracy": round(float(acc), 4)}
        if len(set(y_test)) == 2 and hasattr(model, "predict_proba"):
            proba = model.predict_proba(X_test)[:, 1]
            metrics["roc_auc"] = round(float(roc_auc_score(y_test, proba)), 4)
            metrics["average_precision"] = round(
                float(average_precision_score(y_test, proba)), 4
            )
        return metrics
    pred = model.predict(X_test)
    rmse = float(np.sqrt(mean_squared_error(y_test, pred)))
    return {"r2": round(float(r2_score(y_test, pred)), 4), "rmse": round(rmse, 4)}


def _save_importance(model, n_features, reports):
    importance = np.asarray(model.feature_importances_)
    names = feature_names()[:n_features]
    frame = pd.DataFrame({"feature": names, "importance": importance})
    frame = frame.sort_values("importance", ascending=False)
    frame.to_csv(reports / "ml_feature_importance.csv", index=False)
    top = frame.head(20).iloc[::-1]
    if not top.empty:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.barh(top["feature"], top["importance"], color="#1665c0")
        ax.set_title("ML feature importance")
        fig.tight_layout()
        fig.savefig(reports / "ml_feature_importance.png", dpi=150)
        plt.close(fig)


def _save_roc(model, X_test, y_test, reports):
    proba = model.predict_proba(X_test)[:, 1]
    fpr, tpr, _ = roc_curve(y_test, proba)
    auc = roc_auc_score(y_test, proba)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.plot(fpr, tpr, label=f"AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#9aa5b1")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ML/DL rescoring ROC")
    ax.legend()
    fig.tight_layout()
    fig.savefig(reports / "ml_roc.png", dpi=150)
    plt.close(fig)


def _read_json(path: Path) -> dict:
    import json

    return json.loads(path.read_text(encoding="utf-8"))
