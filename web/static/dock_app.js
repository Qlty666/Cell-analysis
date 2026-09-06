function esc(value) {
  return String(value == null ? '' : value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}
let currentJob = null;
let pollTimer = null;
let jobAlerted = false;
const DOCK_JOB_KEY = 'liver_ui_dock_job';

function saveJobRecord() {
  try { sessionStorage.setItem(DOCK_JOB_KEY, currentJob); } catch (e) {}
}

function clearJobRecord() {
  try { sessionStorage.removeItem(DOCK_JOB_KEY); } catch (e) {}
}

function restoreJobRecord() {
  let saved = null;
  let target = null;
  try { saved = sessionStorage.getItem(DOCK_JOB_KEY); } catch (e) {}
  try { target = new URLSearchParams(window.location.search).get('job'); } catch (e) {}
  const job = target || saved;
  if (!job) return;
  currentJob = job;
  saveJobRecord();
  if (target) {
    try {
      const url = new URL(window.location.href);
      url.searchParams.delete('job');
      window.history.replaceState(null, '', url.toString());
    } catch (e) {}
  }
  msg('已恢复正在运行的虚拟筛选任务，日志继续刷新。', 'ok');
  document.getElementById('pauseBtn').disabled = false;
  document.getElementById('resumeBtn').disabled = true;
  if (pollTimer) clearInterval(pollTimer);
  pollTimer = setInterval(poll, 1000);
}

function showJobToast(kind, title, body) {
  const toast = document.getElementById('jobToast');
  if (!toast) return;
  toast.className = 'job-toast ' + (kind || 'info');
  document.getElementById('jobToastTitle').textContent = title;
  document.getElementById('jobToastBody').textContent = body || '';
  toast.hidden = false;
  jobAlerted = true;
}

function closeJobToast() {
  const toast = document.getElementById('jobToast');
  if (toast) toast.hidden = true;
}

function msg(text, cls) {
  const el = document.getElementById('message');
  el.className = cls || '';
  el.textContent = text;
}

function formValue(name) {
  const form = document.getElementById('form');
  if (!form) return '';
  const fd = new FormData(form);
  return fd.get(name) || '';
}

async function checkEnv(btn) {
  const box = document.getElementById('dockEnvOutput');
  btn.disabled = true;
  box.textContent = '正在检查环境...';
  try {
    await loadDockEnv();
    box.textContent = '环境检查完成。';
  } finally {
    btn.disabled = false;
  }
}

async function detectBox() {
  const workdir = formValue('workdir');
  const receptor = formValue('receptor');
  if (!workdir || !receptor) {
    msg('请先填写工作目录和受体文件路径', 'error');
    return;
  }
  const form = new URLSearchParams({workdir: workdir, receptor: receptor});
  try {
    const resp = await fetch('/dock/detect-box', {method: 'POST', body: form});
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.error || '检测失败');
    const setVal = (id, value) => {
      const el = document.getElementById(id);
      if (el) el.value = value;
    };
    setVal('center_x', data.center[0]);
    setVal('center_y', data.center[1]);
    setVal('center_z', data.center[2]);
    setVal('size_x', data.size[0]);
    setVal('size_y', data.size[1]);
    setVal('size_z', data.size[2]);
    msg('对接盒已自动检测并写入配置。', 'ok');
  } catch (e) {
    msg(String(e.message || e), 'error');
  }
}

async function launchDock(data, btn) {
  btn.disabled = true;
  jobAlerted = false;
  closeJobToast();
  msg('正在启动...');
  try {
    const resp = await fetch('/dock/start', {method: 'POST', body: data});
    const result = await resp.json();
    if (!resp.ok) throw new Error(result.error || '启动失败');
    currentJob = result.job;
    saveJobRecord();
    document.getElementById('log').textContent = '';
    document.getElementById('resultTable').innerHTML = '';
    document.getElementById('gallery').innerHTML = '';
    document.getElementById('summaryRow').textContent = '';
    document.getElementById('fileLinks').textContent = '';
    const mdMessage = document.getElementById('mdMessage');
    if (mdMessage) mdMessage.textContent = '';
    msg('任务已启动，日志实时刷新。', 'ok');
    document.getElementById('pauseBtn').disabled = false;
    document.getElementById('resumeBtn').disabled = true;
    if (pollTimer) clearInterval(pollTimer);
    pollTimer = setInterval(poll, 1000);
  } catch (e) {
    msg(String(e.message || e), 'error');
  } finally {
    btn.disabled = false;
  }
}

function startRun() {
  const form = document.getElementById('form');
  const data = new URLSearchParams(new FormData(form));
  launchDock(data, document.getElementById('startBtn'));
}

async function startMdSimulation(btn) {
  const form = document.getElementById('mdForm');
  const data = new URLSearchParams(new FormData(form));
  data.set('stage', 'md-simulation');
  if (!data.get('md_mode')) data.set('md_mode', 'prepare');
  await launchDock(data, btn || document.getElementById('mdStartBtn'));
}

async function poll() {
  if (!currentJob) return;
  try {
    const logResp = await fetch('/dock/log?job=' + currentJob);
    const text = await logResp.text();
    const box = document.getElementById('log');
    box.textContent = text;
    box.scrollTop = box.scrollHeight;
    const statusResp = await fetch('/dock/status?job=' + currentJob);
    const status = await statusResp.json();
    const pauseBtn = document.getElementById('pauseBtn');
    const resumeBtn = document.getElementById('resumeBtn');
    if (status.queued) {
      pauseBtn.disabled = true;
      resumeBtn.disabled = true;
      msg('任务已进入队列，等待前面的任务完成...', 'ok');
      return;
    }
    if (status.running) {
      pauseBtn.disabled = false;
      resumeBtn.disabled = true;
    } else if (status.paused) {
      pauseBtn.disabled = true;
      resumeBtn.disabled = false;
    } else {
      pauseBtn.disabled = true;
      resumeBtn.disabled = true;
    }
    if (!status.running) {
      clearInterval(pollTimer);
      const m = document.getElementById('message');
      if (status.paused) {
        m.className = 'ok';
        m.textContent = '任务已暂停，可点击继续任务。';
        if (!jobAlerted) {
          showJobToast('paused', '任务已暂停', '当前阶段：' + (status.stage || '未知'));
        }
      } else {
        m.className = status.ok ? 'ok' : 'error';
        if (status.ok) {
          m.textContent = '任务已完成';
          if (!jobAlerted) showJobToast('ok', '任务已完成', '虚拟筛选任务已全部完成。');
        } else {
          m.textContent = '任务失败，请查看日志';
          if (!jobAlerted) {
            showJobToast('error', '任务中断', '运行到阶段：' + (status.stage || '未知') + '\n原因：' + (status.error || '请查看日志'));
          }
        }
      }
      if (status.ok) await loadResults(currentJob);
      if (!status.paused) {
        currentJob = null;
        clearJobRecord();
      }
    }
  } catch (e) {
    // keep polling on transient errors
  }
}

async function pauseDockJob() {
  if (!currentJob) return;
  const form = new URLSearchParams({job: currentJob});
  try {
    const resp = await fetch('/dock/pause', {method: 'POST', body: form});
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.error || '暂停失败');
    msg('暂停请求已发送。', 'ok');
  } catch (e) {
    msg(String(e.message || e), 'error');
  }
}

async function resumeDockJob() {
  if (!currentJob) return;
  const form = new URLSearchParams({job: currentJob});
  try {
    const resp = await fetch('/dock/resume', {method: 'POST', body: form});
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.error || '恢复失败');
    msg('任务已恢复。', 'ok');
    jobAlerted = false;
    closeJobToast();
    document.getElementById('pauseBtn').disabled = false;
    document.getElementById('resumeBtn').disabled = true;
    if (pollTimer) clearInterval(pollTimer);
    pollTimer = setInterval(poll, 1000);
  } catch (e) {
    msg(String(e.message || e), 'error');
  }
}

async function loadResults(job) {
  try {
    const resp = await fetch('/dock/results?job=' + job);
    const data = await resp.json();
    if (data.stage === 'md-simulation') {
      renderMdResults(job, data);
      return;
    }
    const summary = data.summary || {};
    const summaryRow = document.getElementById('summaryRow');
    summaryRow.innerHTML = '<p class="muted">成功对接 ' + esc(summary.total_docked || 0) +
      ' 个配体，命中 ' + esc(summary.hits || 0) +
      ' 个（强 ' + esc(summary.strong_hits || 0) +
      ' / 中 ' + esc(summary.moderate_hits || 0) +
      ' / 弱 ' + esc(summary.weak_hits || 0) +
      '），最佳亲和力 ' + esc(summary.best_affinity === null || summary.best_affinity === undefined ? '-' : summary.best_affinity) +
      ' kcal/mol。输出目录：' + esc(data.output_dir) + '</p>';

    const links = document.getElementById('fileLinks');
    links.innerHTML = (data.files || []).map(name =>
      '<a class="file-link" href="/dock/file?job=' + job + '&name=' + encodeURIComponent(name) + '" download>' + esc(name) + '</a>'
    ).join('');

    const table = document.getElementById('resultTable');
    if (data.rows && data.rows.length) {
      table.innerHTML = '<thead><tr><th>排名</th><th>ID</th><th>亲和力</th><th>分级</th><th>模式</th><th>SMILES</th></tr></thead><tbody>' +
        data.rows.map(r =>
          '<tr><td>' + esc(r.rank || '') + '</td><td>' + esc(r.id) + '</td><td>' +
          esc(r.affinity) + '</td><td>' + esc(r.affinity_class || '') + '</td><td>' + esc(r.mode || '') + '</td><td style="max-width:280px;word-break:break-all;">' +
          esc(r.smiles || '') + '</td></tr>'
        ).join('') + '</tbody>';
    } else {
      table.innerHTML = '<tbody><tr><td>暂无结果</td></tr></tbody>';
    }

    const gallery = document.getElementById('gallery');
    gallery.innerHTML = (data.figures || []).map(name =>
      '<figure><img src="/dock/file?job=' + job + '&name=' + encodeURIComponent(name) + '"><figcaption>' + esc(name) + '</figcaption></figure>'
    ).join('');
  } catch (e) {
    document.getElementById('summaryRow').innerHTML = '<div class="error">结果加载失败：' + esc(e.message || e) + '</div>';
  }
}

function renderMdResults(job, data) {
  const md = data.md || {};
  const summary = md.summary || {};
  const summaryRow = document.getElementById('summaryRow');
  summaryRow.innerHTML = '<p class="muted">模式 ' + esc(summary.mode || '') +
    '，请求 ' + esc(summary.requested || 0) +
    '，完成 ' + esc(summary.completed || 0) +
    '，准备 ' + esc(summary.prepared || 0) +
    '，失败 ' + esc(summary.failed || 0) +
    '。输出目录：' + esc(md.output_dir || '') + '</p>';

  const links = document.getElementById('fileLinks');
  links.innerHTML = (md.files || []).map(name =>
    '<a class="file-link" href="/dock/file?job=' + job + '&name=' + encodeURIComponent(name) + '" download>' + esc(name) + '</a>'
  ).join('');

  const table = document.getElementById('resultTable');
  const headers = ['ID', '状态', '时间 (ns)', '蛋白 RMSD 均值 (nm)', '配体 RMSD 均值 (nm)', '配体 RMSF 均值 (nm)', '结合口袋 RMSF 均值 (nm)', 'Rg 均值 (nm)', 'SASA 均值 (nm²)', '氢键均值', '稳定性', '错误'];
  const keys = ['id', 'status', 'time_ns', 'rmsd_protein_mean_nm', 'rmsd_ligand_mean_nm', 'rmsf_ligand_mean_nm', 'rmsf_contact_residue_mean_nm', 'rg_protein_mean_nm', 'sasa_protein_mean_nm2', 'hbonds_protein_ligand_mean', 'stability_label', 'error'];
  if (md.rows && md.rows.length) {
    table.innerHTML = '<thead><tr>' + headers.map(h => '<th>' + h + '</th>').join('') + '</tr></thead><tbody>' +
      md.rows.map(row => '<tr>' + keys.map(key => {
        let value = row[key];
        const numeric = ['rmsd_protein_mean_nm', 'rmsd_ligand_mean_nm', 'rmsf_ligand_mean_nm', 'rmsf_contact_residue_mean_nm', 'rg_protein_mean_nm', 'sasa_protein_mean_nm2', 'hbonds_protein_ligand_mean'];
        if (numeric.includes(key) &&
            value !== '' && value !== undefined && value !== null) {
          value = Number(value).toFixed(3);
        }
        if (value === '' || value === undefined || value === null) value = '-';
        return '<td>' + esc(value) + '</td>';
      }).join('') + '</tr>').join('') + '</tbody>';
  } else {
    table.innerHTML = '<tbody><tr><td>暂无结果</td></tr></tbody>';
  }

  const gallery = document.getElementById('gallery');
  gallery.innerHTML = (md.figures || []).map(name =>
    '<figure><img src="/dock/file?job=' + job + '&name=' + encodeURIComponent(name) + '"><figcaption>' + esc(name) + '</figcaption></figure>'
  ).join('');
}

async function loadDockHistory() {
  const box = document.getElementById('dockHistory');
  try {
    const resp = await fetch('/dock/history');
    const items = await resp.json();
    if (!items.length) {
      box.textContent = '暂无分析历史。';
      return;
    }
    const rows = items.map(item =>
      '<tr><td>' + esc(item.stage || '') + '</td><td>' + esc(item.status) + '</td><td>' +
      new Date(item.finished * 1000).toLocaleString() + '</td><td>' +
      esc(item.output) + '</td></tr>'
    ).join('');
    box.innerHTML = '<div class="tablebox"><table><thead><tr><th>阶段</th><th>状态</th><th>完成时间</th><th>输出路径</th></tr></thead><tbody>' + rows + '</tbody></table></div>';
  } catch (e) {
    box.textContent = '历史记录加载失败：' + String(e.message || e);
  }
}

async function loadDockEnv() {
  const box = document.getElementById('dockEnv');
  try {
    const resp = await fetch('/dock/check-env', {method: 'POST'});
    const data = await resp.json();
    if (!data.checks || !data.checks.length) {
      box.textContent = '无法获取环境信息。';
      return;
    }
    const badge = c => c.ok ? '<span class="badge ok">OK</span>' : '<span class="badge bad">FAIL</span>';
    const software = data.checks.filter(c => c.kind === 'software');
    const packages = data.checks.filter(c => c.kind !== 'software');
    const softwareRows = software.map(c =>
      '<tr><td>' + esc(c.name) + '</td><td>' + esc(c.version) + '</td><td>' + badge(c) +
      '</td><td><a href="' + esc(c.url) + '" target="_blank" rel="noopener">' + esc(c.url) + '</a></td></tr>'
    ).join('');
    const packageRows = packages.map(c =>
      '<tr><td>' + esc(c.name) + '</td><td>' + esc(c.version) + '</td><td>' + badge(c) +
      '</td><td><code>' + esc(c.install || '') + '</code></td></tr>'
    ).join('');
    let html = '<h3>软件</h3>';
    html += '<div class="tablebox"><table><thead><tr><th>软件</th><th>版本/路径</th><th>状态</th><th>下载地址</th></tr></thead><tbody>' + softwareRows + '</tbody></table></div>';
    html += '<h3>Python 库 / 插件</h3>';
    html += '<div class="tablebox"><table><thead><tr><th>插件</th><th>版本</th><th>状态</th><th>安装命令</th></tr></thead><tbody>' + packageRows + '</tbody></table></div>';
    box.innerHTML = html;
  } catch (e) {
    box.textContent = '环境信息加载失败：' + String(e.message || e);
  }
}

async function installEnv(btn) {
  const box = document.getElementById('dockEnvOutput');
  const target = document.getElementById('dockInstallPath').value.trim();
  btn.disabled = true;
  box.innerHTML = '<p class="muted">正在安装缺失依赖，请稍候...</p>';
  try {
    const form = new URLSearchParams({project: 'dock'});
    if (target) form.set('target', target);
    await fetch('/install', {method: 'POST', body: form});
    const timer = setInterval(async () => {
      try {
        const resp = await fetch('/install-status');
        const data = await resp.json();
        if (!data.running) {
          clearInterval(timer);
          const ok = data.ok;
          box.innerHTML = '<div class="' + (ok ? 'ok' : 'error') + '">' +
            (ok ? '环境补全完成。' : '环境补全失败，请查看日志。') +
            '</div>' + (data.log ? '<pre style="height:220px;">' + esc(data.log) + '</pre>' : '');
          await loadDockEnv();
          btn.disabled = false;
        } else {
          box.innerHTML = '<pre style="height:220px;">' + esc(data.log || '正在安装...') + '</pre>';
        }
      } catch (e) {
        clearInterval(timer);
        box.innerHTML = '<div class="error">状态读取失败：' + esc(e.message || e) + '</div>';
        btn.disabled = false;
      }
    }, 1000);
  } catch (e) {
    box.innerHTML = '<div class="error">环境补全启动失败：' + esc(e.message || e) + '</div>';
    btn.disabled = false;
  }
}

function koMsg(text, cls) {
  const el = document.getElementById('koMessage');
  el.className = cls || '';
  el.textContent = text;
}

async function runKnockout(btn) {
  const form = document.getElementById('koForm');
  const data = new URLSearchParams(new FormData(form));
  btn.disabled = true;
  koMsg('正在运行虚拟敲除...');
  try {
    const resp = await fetch('/dock/knockout', {method: 'POST', body: data});
    const result = await resp.json();
    if (!resp.ok) throw new Error(result.error || '虚拟敲除失败');
    renderKnockout(result);
    koMsg('虚拟敲除完成。', 'ok');
  } catch (e) {
    koMsg(String(e.message || e), 'error');
  } finally {
    btn.disabled = false;
  }
}

async function exportValidation(btn) {
  const form = document.getElementById('koForm');
  const data = new URLSearchParams(new FormData(form));
  if (!data.get('validation_top_n')) data.set('validation_top_n', '10');
  btn.disabled = true;
  const box = document.getElementById('koValidation');
  box.textContent = '正在导出湿实验验证方案...';
  try {
    const resp = await fetch('/dock/knockout/validate', {method: 'POST', body: data});
    const result = await resp.json();
    if (!resp.ok) throw new Error(result.error || '验证方案导出失败');
    const linkBase = '/dock/knockout/file?workdir=' +
      encodeURIComponent(result.workdir) + '&name=';
    box.innerHTML = '<p class="muted">验证方案已导出：' + esc(result.output_dir) + '</p>' +
      (result.files || []).map(name =>
        '<a class="file-link" href="' + linkBase + encodeURIComponent(name) +
        '" download>' + esc(name) + '</a>'
      ).join('');
  } catch (e) {
    box.innerHTML = '<p class="error">' + esc(e.message || e) + '</p>';
  } finally {
    btn.disabled = false;
  }
}

async function loadValidationReport() {
  const box = document.getElementById('validationReport');
  box.textContent = '正在加载验证报告...';
  try {
    const resp = await fetch('/dock/validation-report');
    const data = await resp.json();
    box.textContent = data.report || '暂无报告';
  } catch (e) {
    box.textContent = '加载失败：' + String(e.message || e);
  }
}

async function runValidationJob(btn) {
  const box = document.getElementById('validationReport');
  btn.disabled = true;
  box.textContent = '正在启动随机真实 GSE 验证（10 个数据集，耗时较长）...';
  try {
    const resp = await fetch('/dock/validation/run', {method: 'POST'});
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.error || '启动失败');
    const timer = setInterval(async () => {
      try {
        const sresp = await fetch('/dock/validation-status');
        const status = await sresp.json();
        if (!status.running) {
          clearInterval(timer);
          btn.disabled = false;
          box.textContent = (status.ok ? '验证完成。' : '验证运行失败。') +
            '\n' + (status.log || '');
          if (status.ok) await loadValidationReport();
        } else {
          box.textContent = '验证运行中...\n' + (status.log || '');
        }
      } catch (e) {
        clearInterval(timer);
        btn.disabled = false;
        box.textContent = '状态读取失败：' + String(e.message || e);
      }
    }, 2000);
  } catch (e) {
    btn.disabled = false;
    box.textContent = '启动失败：' + String(e.message || e);
  }
}

function renderKnockout(data) {
  const summary = data.summary || {};
  const summaryEl = document.getElementById('koSummary');
  const classes = summary.target_class_counts || {};
  const classText = Object.keys(classes).map(k =>
    k + ' ' + classes[k]
  ).join('；');
  const insilico = summary.in_silico_knockout || {};
  summaryEl.innerHTML = '<p class="muted">评分基因 ' + esc(summary.genes_scored || 0) +
    ' 个，病例样本 ' + esc(summary.case_samples || 0) +
    '，正常样本 ' + esc(summary.normal_samples || 0) +
    '，DepMap 细胞系 ' + esc(summary.depmap_lines || 0) +
    '，多维评分 ' + (summary.multidimensional_scoring ? '是' : '否') +
    '。分类：' + esc(classText || '无') +
    (insilico.status ? '。单细胞虚拟敲除：' + esc(insilico.ko_gene || '') +
      ' ' + esc(insilico.status || '') : '') +
    '。输出目录：' + esc(data.output_dir) + '</p>';

  const linkBase = '/dock/knockout/file?workdir=' +
    encodeURIComponent(data.workdir) + '&name=';
  document.getElementById('koFiles').innerHTML = (data.files || []).map(name =>
    '<a class="file-link" href="' + linkBase + encodeURIComponent(name) +
    '" download>' + esc(name) + '</a>'
  ).join('');

  const headers = ['排名', '基因', '靶点评分', '分类', '敲除分', '逆转分', '通路分',
    '特异性分', '预后分', '成药性分', 'PPI hub', '安全标记', '优先级'];
  const keys = ['rank', 'gene', 'target_score', 'target_class', 'knockout_score',
    'reversal_score', 'pathway_score', 'specificity_score', 'prognosis_score',
    'druggability_score', 'ppi_hub_score', 'safety_concern', 'target_priority'];
  const scoreKeys = ['target_score', 'knockout_score', 'reversal_score',
    'pathway_score', 'specificity_score', 'prognosis_score', 'druggability_score',
    'ppi_hub_score', 'safety_concern'];
  const table = document.getElementById('koTable');
  if (data.rows && data.rows.length) {
    table.innerHTML = '<thead><tr>' + headers.map(h =>
      '<th>' + h + '</th>'
    ).join('') + '</tr></thead><tbody>' + data.rows.map(row => {
      return '<tr>' + keys.map(key => {
        let value = row[key];
        if (scoreKeys.indexOf(key) >= 0 && value !== '' &&
            value !== undefined && value !== null) {
          value = Number(value).toFixed(3);
        }
        if (value === '' || value === undefined || value === null) value = '-';
        return '<td>' + esc(value) + '</td>';
      }).join('') + '</tr>';
    }).join('') + '</tbody>';
  } else {
    table.innerHTML = '<tbody><tr><td>暂无结果</td></tr></tbody>';
  }

  document.getElementById('koGallery').innerHTML = (data.figures || []).map(name =>
    '<figure><img src="' + linkBase + encodeURIComponent(name) +
    '"><figcaption>' + esc(name) + '</figcaption></figure>'
  ).join('');
}

function netMsg(text, cls) {
  const el = document.getElementById('netMessage');
  el.className = cls || '';
  el.textContent = text;
}

async function runNetwork(btn) {
  const form = document.getElementById('netForm');
  const data = new URLSearchParams(new FormData(form));
  btn.disabled = true;
  netMsg('正在运行网络毒理学分析...');
  try {
    const resp = await fetch('/dock/network', {method: 'POST', body: data});
    const result = await resp.json();
    if (!resp.ok) throw new Error(result.error || '网络毒理学分析失败');
    renderNetwork(result);
    netMsg('网络毒理学分析完成。', 'ok');
  } catch (e) {
    netMsg(String(e.message || e), 'error');
  } finally {
    btn.disabled = false;
  }
}

function renderNetwork(data) {
  const summary = data.summary || {};
  document.getElementById('netSummary').innerHTML =
    '<p class="muted">化合物靶点 ' + esc(summary.compound_targets || 0) +
    ' 个，疾病基因 ' + esc(summary.disease_genes || 0) +
    '，交集基因 ' + esc(summary.overlap_genes || 0) +
    '，PPI hub 评分 ' + (summary.ppi_hub_scored ? '已启用' : '未启用') +
    '。输出目录：' + esc(data.output_dir) + '</p>';

  const linkBase = '/dock/network/file?workdir=' +
    encodeURIComponent(data.workdir) + '&name=';
  document.getElementById('netFiles').innerHTML = (data.files || []).map(name =>
    '<a class="file-link" href="' + linkBase + encodeURIComponent(name) +
    '" download>' + esc(name) + '</a>'
  ).join('');

  const headers = ['基因', '来源数', '来源', 'PPI degree', 'PPI betweenness',
    'PPI clustering', 'PPI hub'];
  const keys = ['gene', 'n_sources', 'sources', 'ppi_degree',
    'ppi_betweenness', 'ppi_clustering', 'ppi_hub_score'];
  const table = document.getElementById('netTable');
  if (data.rows && data.rows.length) {
    table.innerHTML = '<thead><tr>' + headers.map(h =>
      '<th>' + h + '</th>'
    ).join('') + '</tr></thead><tbody>' + data.rows.map(row => {
      return '<tr>' + keys.map(key => {
        let value = row[key];
        if (value === '' || value === undefined || value === null) value = '-';
        return '<td>' + esc(value) + '</td>';
      }).join('') + '</tr>';
    }).join('') + '</tbody>';
  } else {
    table.innerHTML = '<tbody><tr><td>暂无交集基因</td></tr></tbody>';
  }
  document.getElementById('netGallery').innerHTML = (data.figures || []).map(name =>
    '<figure><img src="' + linkBase + encodeURIComponent(name) +
    '"><figcaption>' + esc(name) + '</figcaption></figure>'
  ).join('');
}

function faersMsg(text, cls) {
  const el = document.getElementById('faersMessage');
  el.className = cls || '';
  el.textContent = text;
}

async function runFaers(btn) {
  const form = document.getElementById('faersForm');
  const data = new URLSearchParams(new FormData(form));
  btn.disabled = true;
  faersMsg('正在运行 FAERS 信号检测...');
  try {
    const resp = await fetch('/dock/faers', {method: 'POST', body: data});
    const result = await resp.json();
    if (!resp.ok) throw new Error(result.error || 'FAERS 信号检测失败');
    renderFaers(result);
    faersMsg('FAERS 信号检测完成。', 'ok');
  } catch (e) {
    faersMsg(String(e.message || e), 'error');
  } finally {
    btn.disabled = false;
  }
}

function renderFaers(data) {
  const summary = data.summary || {};
  document.getElementById('faersSummary').innerHTML =
    '<p class="muted">药物-事件组合 ' + esc(summary.pairs || 0) +
    ' 个，信号 ' + esc(summary.signals || 0) +
    ' 个，最小计数 ' + esc(summary.min_count || 3) +
    '。输出目录：' + esc(data.output_dir) + '</p>';

  const linkBase = '/dock/faers/file?workdir=' +
    encodeURIComponent(data.workdir) + '&name=';
  document.getElementById('faersFiles').innerHTML = (data.files || []).map(name =>
    '<a class="file-link" href="' + linkBase + encodeURIComponent(name) +
    '" download>' + esc(name) + '</a>'
  ).join('');

  const headers = ['药物', '事件', '计数', 'ROR', 'ROR 下限', 'PRR',
    'IC', 'EBGM', '信号'];
  const keys = ['drug', 'event', 'a', 'ror', 'ror_lower', 'prr',
    'ic', 'ebgm', 'signal'];
  const table = document.getElementById('faersTable');
  if (data.rows && data.rows.length) {
    table.innerHTML = '<thead><tr>' + headers.map(h =>
      '<th>' + h + '</th>'
    ).join('') + '</tr></thead><tbody>' + data.rows.map(row => {
      return '<tr>' + keys.map(key => {
        let value = row[key];
        if (key === 'ror' || key === 'ror_lower' || key === 'prr' ||
            key === 'ic' || key === 'ebgm') {
          if (value !== '' && value !== undefined && value !== null) {
            value = Number(value).toFixed(3);
          }
        }
        if (value === '' || value === undefined || value === null) value = '-';
        return '<td>' + esc(value) + '</td>';
      }).join('') + '</tr>';
    }).join('') + '</tbody>';
  } else {
    table.innerHTML = '<tbody><tr><td>暂无信号结果</td></tr></tbody>';
  }
}

function initDockPage() {
  loadDockHistory();
  loadDockEnv();
  restoreJobRecord();
}

function initValidationPage() {
  loadValidationReport();
}
