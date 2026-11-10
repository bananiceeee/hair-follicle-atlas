import {decodeCells, decodeGene, findGene} from './data.mjs';

const select = (s) => document.querySelector(s);
const all = (s) => Array.from(document.querySelectorAll(s));
const geneInput = select('#gene');
const status = select('#status');
const plots = select('#plots');
const legend = select('#legend');
const summary = select('#expression-summary');
const query = new URLSearchParams(location.search);
let atlas = query.get('atlas') === 'follicle' ? 'follicle' : 'skin';
let mode = ['split','annotation'].includes(query.get('view')) ? query.get('view') : 'expression';
let revision = 0;
const atlasCache = new Map();
let expressionCache;
let geneSets;
const names = {skin:'Whole skin', follicle:'Hair follicle subset'};

async function load(path, binary = false) {
  const response = await fetch(path);
  if (!response.ok) throw Error('Could not load atlas data');
  return binary ? response.arrayBuffer() : response.text();
}
async function getAtlas(kind) {
  if (!atlasCache.has(kind)) {
    const base = './data/' + kind + '/';
    const promise = Promise.all(['format.tsv','genes.txt','cell-types.tsv','cells.bin']
      .map((file) => load(base + file, file.endsWith('.bin'))))
      .then(([format,genes,types,cells]) => {
        const [size,chunk,scale] = format.trim().split('\t').map(Number);
        const labels = types.trim().split(/\r?\n/).slice(1).map((line) => line.split('\t'));
        return { ...decodeCells(cells,size,labels), chunk, scale, base, genes:genes.trim().split(/\r?\n/) };
      }).catch((e) => { atlasCache.delete(kind); throw e; });
    atlasCache.set(kind,promise);
  }
  return atlasCache.get(kind);
}
async function getExpression(data,index) {
  if (expressionCache?.data === data && expressionCache.index === index) return expressionCache.result;
  const block = Math.floor(index / data.chunk);
  const width = Math.min(data.chunk, data.genes.length - block * data.chunk);
  const result = await decodeGene(await load(data.base + String(block).padStart(4,'0') + '.bin', true),
    index % data.chunk,width,data.cells,data.scale);
  expressionCache = {data,index,result};
  return result;
}

const colorStops = [[226,229,229],[69,117,180],[232,178,77],[191,46,38]];
const colors = Array.from({length:101}, (_,i) => {
  const f = i / 100 * (colorStops.length-1), j = Math.min(2,Math.floor(f));
  return 'rgb(' + colorStops[j].map((v,k) => Math.round(v+(f-j)*(colorStops[j+1][k]-v))).join(',') + ')';
});
function node(tag, text, css) {
  const item = document.createElement(tag);
  if (text) item.textContent = text;
  if (css) item.className = css;
  return item;
}
function umap(data,result,condition,title) {
  const figure = node('figure');
  const canvas = node('canvas');
  canvas.width = 720;
  canvas.height = 620;
  canvas.setAttribute('aria-label',title);
  canvas.setAttribute('role','img');
  figure.append(canvas,node('figcaption',title));
  const c = canvas.getContext('2d');
  const [xmin,xmax,ymin,ymax] = data.bounds;
  const xRange = xmax - xmin, yRange = ymax - ymin;
  const scale = Math.min((canvas.width-68)/xRange,(canvas.height-64)/yRange);
  const x0 = (canvas.width-xRange*scale)/2, y0 = (canvas.height-yRange*scale)/2;
  const radius = data.cells > 10000 ? 1.65 : 2.05;
  const value = result?.expression;
  function draw(positive) {
    for (let i = 0; i < data.cells; ++i) {
      if (condition !== null && data.condition[i] !== condition) continue;
      if (value && (value[i] > 0) !== positive) continue;
      if (!value && positive) continue;
      c.fillStyle = value ? (positive ? colors[Math.max(1,Math.min(100,Math.round(value[i]/result.cap*100)))] : '#e2e5e5')
        : data.types[data.type[i]][1];
      const x = x0+(data.x[i]-xmin)*scale, y = y0+(ymax-data.y[i])*scale;
      c.beginPath(); c.arc(x,y,radius,0,Math.PI*2); c.fill();
    }
  }
  draw(false);
  if (value) draw(true);
  c.strokeStyle = '#344844'; c.lineWidth = 1.4; c.fillStyle = '#344844'; c.font = '14px system-ui';
  c.beginPath(); c.moveTo(23,572); c.lineTo(23,539); c.moveTo(23,572); c.lineTo(62,572); c.stroke();
  c.fillText('UMAP 1',23,591);
  c.save(); c.translate(14,571); c.rotate(-Math.PI/2); c.fillText('UMAP 2',0,0); c.restore();
  plots.append(figure);
}
function drawLegend(data,result) {
  if (!result) {
    for (const [name,color] of data.types) {
      const key = node('span'); const dot = node('span','','swatch'); dot.style.background = color;
      key.append(dot,document.createTextNode(name)); legend.append(key);
    }
  } else if (result.detected) {
    legend.append(node('span','0'),node('span','','gradient'),node('span',result.cap.toFixed(2)),
      node('span','Log-normalized expression'));
  } else legend.append(node('span','Grey · no detected expression'));
}
async function drawMembership(symbol,id) {
  const message = select('#membership-status');
  const grid = select('#membership-grid');
  grid.replaceChildren();
  message.textContent = 'Checking ' + symbol + '…';
  try {
    geneSets ??= load('./data/goat-top50.json').then(JSON.parse).catch((error) => {
      geneSets = undefined;
      throw error;
    });
    const sets = await geneSets;
    if (id !== revision) return;
    const key = symbol.toUpperCase();
    const count = sets.filter(({members}) => Object.hasOwn(members,key)).length;
    message.textContent = symbol + ' · ' + count + ' of ' + sets.length + ' top-50 sets by symbol';
    for (const {label,members} of sets) {
      const member = Object.hasOwn(members,key);
      const state = member ? (members[key] ? 'Member · matched in goat skin' : 'Member · no goat-skin match') : 'Outside top 50';
      const item = node('div','','membership-item');
      item.append(node('strong',label),node('span',state,member ? 'membership-hit' : ''));
      grid.append(item);
    }
  } catch {
    if (id === revision) message.textContent = 'Gene-set membership is temporarily unavailable.';
  }
}

function downloadSummary(gene, kind, types, groups) {
  const rows = [['gene','atlas','cell_type','condition','total_cells','detected_cells','percent_detected']];
  for (const [condition,counts] of groups) {
    counts.forEach(([total,detected],i) => {
      rows.push([gene,kind,types[i][0],condition,total,detected,
        total ? (100 * detected / total).toFixed(6) : '']);
    });
  }
  const escape = (value) => '"' + String(value).replaceAll('"','""') + '"';
  const csv = rows.map((row) => row.map(escape).join(',')).join('\r\n') + '\r\n';
  const url = URL.createObjectURL(new Blob([csv],{type:'text/csv;charset=utf-8'}));
  const link = node('a');
  link.href = url;
  link.download = gene.replace(/[^a-z0-9._-]/gi,'_') + '_' + kind + '_expression_detection.csv';
  link.click();
  setTimeout(() => URL.revokeObjectURL(url),1000);
}

function drawSummary(data, result, split, gene, kind) {
  const counts = data.types.map(() => [[0,0],[0,0]]);
  for (let i = 0; i < data.cells; ++i) {
    const entry = counts[data.type[i]][data.condition[i]];
    ++entry[0];
    if (result.expression[i] > 0) ++entry[1];
  }
  const allCounts = counts.map(([low,high]) => low.map((value,j) => value + high[j]));
  const splitCounts = [['L',counts.map(([low]) => low)],['H',counts.map(([,high]) => high)]];
  const groups = split ? splitCounts : [['All',allCounts]];
  const grid = node('div','','summary-grid');
  for (const [label,groupCounts] of groups) {
    const column = node('div');
    if (split) column.append(node('h4',label));
    data.types.forEach(([name,color], i) => {
      const [total,detected] = groupCounts[i];
      const percent = total ? 100 * detected / total : 0;
      const row = node('div','','summary-row');
      const track = node('span','','summary-track');
      const bar = node('span','','summary-bar');
      bar.style.width = percent + '%';
      bar.style.backgroundColor = color;
      track.append(bar);
      row.title = name + ' · ' + detected + ' / ' + total + ' cells';
      row.append(node('span',name),track,node('span',total ? percent.toFixed(1) + '%' : '–'));
      column.append(row);
    });
    grid.append(column);
  }
  const heading = node('div','','summary-heading');
  const download = node('button','Download CSV');
  download.type = 'button';
  download.addEventListener('click',() =>
    downloadSummary(gene,kind,data.types,[['All',allCounts],...splitCounts]));
  heading.append(node('h3','Cells with detected expression by cell type'),download);
  summary.append(heading,grid);
  summary.hidden = false;
}

async function render() {
  const id = ++revision, kind = atlas, view = mode;
  for (const b of all('[data-atlas]')) b.setAttribute('aria-pressed',b.dataset.atlas === kind);
  for (const b of all('[data-view]')) b.setAttribute('aria-pressed',b.dataset.view === view);
  status.textContent = 'Loading…';
  drawMembership(geneInput.value.trim(),id);
  plots.replaceChildren(); legend.replaceChildren();
  summary.replaceChildren(); summary.hidden = true;
  try {
    const data = await getAtlas(kind);
    if (id !== revision) return;
    const index = findGene(data.genes,geneInput.value);
    if (index === -1 && view !== 'annotation') {
      select('#plot-title').textContent = 'Expression UMAP';
      status.textContent = 'Gene not found in this atlas. Try an official mouse gene symbol.';
      return;
    }
    const result = view === 'annotation' ? null : await getExpression(data,index);
    if (id !== revision) return;
    const canonical = data.genes[index];
    if (canonical) geneInput.value = canonical;
    const url = new URL(location.href);
    if (canonical) url.searchParams.set('gene',canonical);
    url.searchParams.set('atlas',kind); url.searchParams.set('view',view);
    history.replaceState(null,'',url);
    select('#plot-title').textContent = result ? canonical + ' expression' : 'Cell types';
    status.textContent = names[kind] + ' · ' + data.cells.toLocaleString('en-US') + ' cells' +
      (result ? ' · ' + result.detected.toLocaleString('en-US') + ' with detected expression (' +
        (result.detected/data.cells*100).toFixed(1) + '%)' : '');
    if (view === 'split') {
      umap(data,result,0,'L · Low density');
      umap(data,result,1,'H · High density');
    } else umap(data,result,null,names[kind]);
    drawLegend(data,result);
    if (result) drawSummary(data,result,view === 'split',canonical,kind);
  } catch (error) {
    if (id === revision) status.textContent = error.message + '. Please reload or try again.';
  }
}
async function suggest() {
  const kind = atlas, q = geneInput.value.trim().toLowerCase(), list = select('#genes');
  if (!q) { list.replaceChildren(); return; }
  try {
    const data = await getAtlas(kind);
    if (kind !== atlas || q !== geneInput.value.trim().toLowerCase()) return;
    list.replaceChildren(...data.genes.filter((gene) => gene.toLowerCase().startsWith(q)).slice(0,35)
      .map((gene) => { const option = node('option'); option.value = gene; return option; }));
  } catch { /* The main panel reports data errors. */ }
}
select('#search').addEventListener('submit',(event) => {
  event.preventDefault(); if (mode === 'annotation') mode = 'expression'; render();
});
geneInput.addEventListener('input',suggest);
for (const b of all('[data-gene]')) b.addEventListener('click',() => {
  geneInput.value = b.dataset.gene; if (mode === 'annotation') mode = 'expression'; render();
});
for (const b of all('[data-atlas]')) b.addEventListener('click',() => { atlas = b.dataset.atlas; render(); });
for (const b of all('[data-view]')) b.addEventListener('click',() => { mode = b.dataset.view; render(); });
geneInput.value = query.get('gene') || 'Krt79';
render();
