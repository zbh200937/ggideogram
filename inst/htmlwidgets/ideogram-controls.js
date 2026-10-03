function(el, x, data) {
  const svg = el.querySelector('svg');
  if (!svg) return;
  svg.style.width = data.width + 'mm';
  svg.style.height = data.height + 'mm';
  svg.parentElement.style.width = data.width + 'mm';
  svg.parentElement.style.height = 'auto';
  el.style.width = data.width + 'mm';
  el.style.maxWidth = '100%';
  el.style.height = 'auto';
  el.style.overflowX = 'auto';
  el.style.font = '11px/1.5 Arial, sans-serif';
  const records = Array.isArray(data.records) ? data.records : [];
  if (!records.length) return;
  const byId = new Map(records.map(record => [String(record.id), record]));
  const selected = new Set();
  const marks = Array.from(svg.querySelectorAll('[data-id]'));
  const controls = document.createElement('div');
  controls.className = 'ideogram-controls';
  controls.style.cssText = 'display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:8px 0;color:#202020;';
  const status = document.createElement('span');
  status.setAttribute('role', 'status');
  const table = document.createElement('table');
  table.className = 'ideogram-selected-records';
  table.style.cssText = 'border-collapse:collapse;margin:6px 0;width:100%;text-align:left;';
  const head = table.createTHead().insertRow();
  ['ID', 'Chr', 'Start', 'End'].forEach(label => {
    const cell = document.createElement('th');
    cell.textContent = label;
    cell.style.cssText = 'padding:3px 10px;border-bottom:1px solid #d7dde2;';
    head.appendChild(cell);
  });
  const body = table.createTBody();
  const style = document.createElement('style');
  style.textContent = '.ideogram-controls button,.ideogram-controls select,.ideogram-controls a{font:inherit;padding:3px 7px;background:white;border:1px solid #bac4cd;border-radius:3px;color:#202020;text-decoration:none}.ideogram-controls button:disabled,.ideogram-controls a[aria-disabled="true"]{color:#8c969f}.ideogram-selected:not(text){stroke:#202020!important;stroke-width:1.6px!important}';
  el.appendChild(style);
  el.prepend(controls);
  el.appendChild(table);
  let filter = '';
  const refresh = () => {
    marks.forEach(mark => {
      const id = mark.getAttribute('data-id');
      const record = byId.get(id);
      if (!record) return;
      mark.style.display = !filter || String(record.category) === filter ? '' : 'none';
      mark.classList.toggle('ideogram-selected', selected.has(id));
    });
    body.replaceChildren();
    records.filter(record => selected.has(String(record.id))).forEach(record => {
      const row = body.insertRow();
      [record.id, record.chr, record.start, record.end].forEach((value, i) => {
        const cell = row.insertCell();
        cell.style.cssText = 'padding:3px 10px;';
        if (i === 0 && record.url) {
          const link = document.createElement('a');
          link.textContent = value;
          link.href = record.url;
          link.target = '_blank';
          link.rel = 'noopener';
          cell.appendChild(link);
        } else cell.textContent = value;
      });
    });
    table.hidden = !selected.size;
    status.textContent = '已选 ' + selected.size + ' 个区间';
    download.setAttribute('aria-disabled', String(!selected.size));
    if (selected.size) {
      const rows = records.filter(record => selected.has(String(record.id)));
      const content = [data.fields.map(tsv).join('\t')].concat(
        rows.map(record => record.values.map(tsv).join('\t'))).join('\n') + '\n';
      download.href = 'data:text/tab-separated-values;charset=utf-8,' + encodeURIComponent(content);
    } else download.removeAttribute('href');
    clear.disabled = !selected.size;
  };
  if (data.filter) {
    const label = document.createElement('label');
    label.textContent = '显示类别 ';
    const select = document.createElement('select');
    select.setAttribute('aria-label', '显示类别');
    const all = new Option('全部类别', '');
    select.appendChild(all);
    Array.from(new Set(records.map(record => String(record.category)))).forEach(category => {
      select.appendChild(new Option(category, category));
    });
    select.addEventListener('change', () => {
      filter = select.value;
      selected.clear();
      refresh();
    });
    label.appendChild(select);
    controls.appendChild(label);
  }
  const clear = document.createElement('button');
  clear.type = 'button';
  clear.textContent = '清空选择';
  clear.addEventListener('click', () => { selected.clear(); refresh(); });
  const download = document.createElement('a');
  download.setAttribute('role', 'button');
  download.download = 'ideogram-selected.tsv';
  download.textContent = '导出所选区间';
  const tsv = value => {
    if (value === null || value === undefined) return '';
    const text = typeof value === 'object' ? JSON.stringify(value) : String(value);
    return /[\t\r\n"]/.test(text) ? '"' + text.replace(/"/g, '""') + '"' : text;
  };
  controls.append(clear, download, status);
  svg.addEventListener('click', event => {
    const mark = event.target.closest('[data-id]');
    if (!mark || !svg.contains(mark)) return;
    const id = mark.getAttribute('data-id');
    if (!byId.has(id)) return;
    if (selected.has(id)) selected.delete(id); else selected.add(id);
    refresh();
  });
  refresh();
}
