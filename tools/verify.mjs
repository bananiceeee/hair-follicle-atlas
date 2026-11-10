import assert from 'node:assert/strict';
import {existsSync, readdirSync, readFileSync, statSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {decodeCells, decodeGene, findGene} from '../docs/data.mjs';

const examples = ['Krt79', 'Il11ra1', 'Nr4a1'];
const toArrayBuffer = (b) => b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);
const metadata = (dir) => {
  const read = (name) => readFileSync(dir + name, 'utf8').trim();
  const [cells, chunk, scale] = read('format.tsv').split('\t').map(Number);
  const genes = read('genes.txt').split(/\r?\n/);
  const types = read('cell-types.tsv').split(/\r?\n/).slice(1).map((line) => line.split('\t'));
  return {cells, chunk, scale, genes, types, dir};
};
const shard = (index) => String(index).padStart(4, '0') + '.bin';
async function expression(data, gene) {
  const index = findGene(data.genes, gene);
  assert(index >= 0, gene);
  const block = Math.floor(index / data.chunk);
  const width = Math.min(data.chunk, data.genes.length - block * data.chunk);
  const file = readFileSync(data.dir + shard(block));
  return [await decodeGene(toArrayBuffer(file), index % data.chunk, width, data.cells, data.scale), file.length];
}

for (const [kind, expectedCells] of [['skin', 27051], ['follicle', 6026]]) {
  const dir = 'docs/data/' + kind + '/';
  const data = metadata(dir);
  assert.equal(data.cells, expectedCells);
  assert(data.genes.length > 10000);
  assert.equal(new Set(data.genes).size, data.genes.length);
  assert.equal(findGene(data.genes, '__unknown_gene__'), -1);
  const cells = decodeCells(toArrayBuffer(readFileSync(dir + 'cells.bin')), data.cells, data.types);
  assert.equal(new Set(cells.condition).size, 2);
  assert.equal(new Set(cells.type).size, data.types.length);
  const files = readdirSync(dir);
  const blocks = Math.ceil(data.genes.length / data.chunk);
  assert.equal(files.filter((file) => /^\d+\.bin$/.test(file)).length, blocks);
  let largest = 0;
  for (let i = 0; i < blocks; i++) {
    const zipped = readFileSync(dir + shard(i));
    largest = Math.max(largest, zipped.length);
    const bytes = gunzipSync(zipped);
    const width = Math.min(data.chunk, data.genes.length - i * data.chunk);
    let offset = width * 4;
    for (let j = 0; j < width; j++) {
      const count = bytes.readUInt32LE(j * 4);
      assert(count <= data.cells && offset + count * 4 <= bytes.length, shard(i));
      let previous = -1;
      for (let k = 0; k < count; k++) {
        const cell = bytes.readUInt16LE(offset + k * 2);
        const value = bytes.readUInt16LE(offset + (count + k) * 2);
        assert(cell > previous && cell < data.cells && value <= Math.round(Math.log1p(10000) * data.scale), shard(i));
        previous = cell;
      }
      offset += count * 4;
    }
    assert.equal(offset, bytes.length, shard(i));
  }
  const demoPath = 'preview/data/' + kind + '/';
  const demo = existsSync(demoPath + 'format.tsv') ? metadata(demoPath) : null;
  const checked = {};
  for (const gene of examples) {
    const [full, size] = await expression(data, gene.toUpperCase());
    if (demo) {
      const [sample] = await expression(demo, gene);
      assert.deepEqual(full, sample);
    }
    checked[gene] = {detected: full.detected, shardBytes: size};
  }
  const bytes = files.reduce((sum, file) => sum + statSync(dir + file).size, 0);
  console.log(JSON.stringify({kind, cells: data.cells, genes: data.genes.length, types: data.types.length,
    low: cells.condition.filter((x) => x === 0).length, blocks, bytes, largest, checked}));
}
