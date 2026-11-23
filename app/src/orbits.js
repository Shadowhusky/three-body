const SVG = 'http://www.w3.org/2000/svg';

function glyph(thumb, colors) {
  const svg = document.createElementNS(SVG, 'svg');
  svg.setAttribute('viewBox', '-1.08 -1.08 2.16 2.16');
  svg.setAttribute('aria-hidden', 'true');
  thumb.forEach((pts, b) => {
    const path = document.createElementNS(SVG, 'polyline');
    path.setAttribute('points', pts.join(' '));
    path.setAttribute('fill', 'none');
    path.setAttribute('stroke', colors[b]);
    path.setAttribute('stroke-width', '0.03');
    path.setAttribute('stroke-opacity', '0.9');
    svg.append(path);
  });
  return svg;
}

const GROUPS = [
  ['flat', 'Planar orbits', 'Zero angular momentum. All three bodies move in one plane.'],
  [
    'family',
    'One family',
    'A continuous family of 3D orbits, each periodic in a rotating frame. Every entry is the same pattern rotating at a different rate; any fraction can be computed on request.',
  ],
  ['spatial', '3D orbits', 'Systems with angular momentum whose paths leave the plane.'],
  ['orderly', 'Quasi-periodic orbits', 'They never repeat, yet stay regular. A formula fitted to a short span predicts them far ahead.'],
];

let rows = new Map();

export function buildOrbitList(list, orbits, colors, onPick) {
  const onlyStable = document.getElementById('only-stable');
  const count = document.getElementById('orbit-count');

  const render = () => {
    list.replaceChildren();
    rows = new Map();
    let shown = 0;
    for (const [kind, title, note] of GROUPS) {
      const group = orbits.filter((o) => o.kind === kind && (!onlyStable.checked || o.stable));
      if (!group.length) continue;
      const head = document.createElement('li');
      head.innerHTML = `<h2>${title}</h2><p class="group-note">${note}</p>`;
      list.append(head);
      for (const o of group) {
        const li = document.createElement('li');
        const row = document.createElement('button');
        row.className = 'orbit';
        row.append(glyph(o.thumb, colors));
        const text = document.createElement('span');
        text.innerHTML = `<span class="name"></span><br><span class="meta"></span>`;
        text.querySelector('.name').textContent = o.name;
        text.querySelector('.meta').textContent = [o.period ? (o.stable ? 'stable' : 'unstable') : '', o.tag].filter(Boolean).join(', ');
        const period = document.createElement('span');
        period.className = 'period';
        period.textContent = o.period ? `period ${o.period.toFixed(2)}` : '';
        row.append(text, period);
        row.addEventListener('click', () => onPick(o));
        li.append(row);
        list.append(li);
        rows.set(o.id, row);
        shown++;
      }
    }
    count.textContent = orbits.length ? `${shown} of ${orbits.length}` : 'No orbits loaded';
  };

  onlyStable.addEventListener('change', render);
  render();
}

export function markCurrentOrbit(id) {
  rows.forEach((row, key) => row.setAttribute('aria-current', String(key === id)));
}
