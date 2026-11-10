const integer = (view, offset, width) => width === 2
  ? view.getUint16(offset, true) : view.getUint32(offset, true);

export function decodeCells(buffer, cells, types) {
  if (buffer.byteLength !== cells * 10) throw Error('Invalid cell coordinates');
  const view = new DataView(buffer);
  const x = new Float32Array(cells), y = new Float32Array(cells);
  const type = new Uint8Array(buffer, cells * 8, cells);
  const condition = new Uint8Array(buffer, cells * 9, cells);
  for (let i = 0; i < cells; ++i) {
    x[i] = view.getFloat32(i * 4, true);
    y[i] = view.getFloat32((cells + i) * 4, true);
    if (!Number.isFinite(x[i]) || !Number.isFinite(y[i]) || type[i] >= types.length || condition[i] > 1)
      throw Error('Invalid cell metadata');
  }
  return { x, y, type, condition, types, cells,
    bounds: [Math.min(...x), Math.max(...x), Math.min(...y), Math.max(...y)] };
}

export async function decodeGene(buffer, index, width, cells, scale) {
  const stream = new Blob([buffer]).stream().pipeThrough(new DecompressionStream('gzip'));
  const view = new DataView(await new Response(stream).arrayBuffer());
  const count = integer(view, index * 4, 4);
  let offset = width * 4;
  for (let j = 0; j < index; ++j) offset += integer(view, j * 4, 4) * 4;
  if (count > cells || offset + count * 4 > view.byteLength) throw Error('Invalid expression block');
  const expression = new Float32Array(cells);
  const positive = [];
  let previous = -1;
  for (let j = 0; j < count; ++j) {
    const cell = integer(view, offset + j * 2, 2);
    if (cell >= cells || cell <= previous) throw Error('Invalid expression index');
    previous = cell;
    const value = integer(view, offset + (count + j) * 2, 2) / scale;
    expression[cell] = value;
    if (value > 0) positive.push(value);
  }
  positive.sort((a, b) => a - b);
  const q = Math.max(0, (positive.length - 1) * .99);
  const cap = positive.length ? positive[Math.floor(q)] +
    (q % 1) * (positive[Math.ceil(q)] - positive[Math.floor(q)]) : 0;
  return { expression, detected: positive.length, cap };
}

export function findGene(genes, input) {
  const query = input.trim();
  const exact = genes.indexOf(query);
  if (exact >= 0) return exact;
  const matching = genes.reduce((acc, gene, i) => {
    if (gene.toLowerCase() === query.toLowerCase()) acc.push(i);
    return acc;
  }, []);
  return matching.length === 1 ? matching[0] : -1;
}
