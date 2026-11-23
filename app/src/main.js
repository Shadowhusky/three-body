import { Stage, BODY_COLORS } from './scene.js';
import { ThreeBody, energy, minSeparation, signedArea, middleBody, escaper } from './physics.js';
import { buildOrbitList, markCurrentOrbit } from './orbits.js';
import { MapPlate } from './map.js';
import findingsHtml from './findings.html?raw';

const $ = (id) => document.getElementById(id);
const stage = new Stage($('stage'));

const EQUAL = [1, 1, 1];
const TAPE_LENGTH = 34;

const app = {
  system: null,
  sim: null,
  playing: true,
  speed: 1,
  trail: 0.5,
  baseSpan: 12,
  maxGap: 0.01,
  e0: 0,
  normal: [0, 0, 1],
  areaSign: 0,
  gone: -1,
  showFormula: false,
  frame: 0,
};

export function collinearStart(p1, p2) {
  return [-1, 0, 0, 1, 0, 0, 0, 0, 0, p1, p2, 0, p1, p2, 0, -2 * p1, -2 * p2, 0];
}

function spinAxis(s, m) {
  const L = [0, 0, 0];
  for (let b = 0; b < 3; b++) {
    const [x, y, z] = [s[3 * b], s[3 * b + 1], s[3 * b + 2]];
    const [u, v, w] = [s[9 + 3 * b], s[10 + 3 * b], s[11 + 3 * b]];
    L[0] += m[b] * (y * w - z * v);
    L[1] += m[b] * (z * u - x * w);
    L[2] += m[b] * (x * v - y * u);
  }
  const n = Math.hypot(...L);
  return n > 1e-9 ? L.map((c) => c / n) : [0, 0, 1];
}

// Run the system ahead once to learn how much room it needs and, for a
// repeating orbit, the full closed path.
function preview(system) {
  const sim = new ThreeBody(system.masses);
  sim.setState(system.state);
  const horizon = system.period || 25;
  const samples = [];
  const reach = [];
  sim.advanceTo(horizon, horizon / 1500, (t, s) => {
    let r = 0;
    for (let b = 0; b < 3; b++) r = Math.max(r, Math.hypot(s[3 * b], s[3 * b + 1], s[3 * b + 2]));
    reach.push(r);
    if (system.period) samples.push(Float64Array.from(s.subarray(0, 9)));
  }, 400000);
  reach.sort((a, b) => a - b);
  // a repeating orbit is framed whole; a wild one by where it spends most of its time
  const extent = system.period ? reach[reach.length - 1] : Math.min(Math.max(2 * reach[reach.length >> 1], 1.3), 2.6);
  return { extent: Math.max(extent, 0.5), samples: system.period ? samples : null };
}

export function loadSystem(system) {
  system.masses ||= EQUAL;
  app.system = system;
  app.sim = new ThreeBody(system.masses);
  app.sim.setState(system.state);
  app.e0 = energy(system.state, system.masses);
  app.normal = spinAxis(system.state, system.masses);
  app.areaSign = 0;
  app.gone = -1;

  const { extent, samples } = preview(system);
  app.baseSpan = system.period ? system.period * 0.92 : 12;
  app.maxGap = (system.period || 25) / 2500;
  stage.frame(extent, system.masses);
  stage.setPath(samples);
  stage.clearTrails();
  stage.pushTrail(0, system.state);
  applyTrail();

  $('system-name').textContent = system.name;
  $('system-note').textContent = system.note || '';
  $('r-fate').textContent = '';
  $('tape-beads').replaceChildren();
  $('formula-toggle').hidden = !system.formula;
  $('r-formula-row').hidden = !(system.formula && app.showFormula);
  stage.setGhosts(null);
  markCurrentOrbit(system.id);
  mapPlate?.setCurrent(system.mapPoint || null);
  setPlaying(true);
}

function applyTrail() {
  stage.span = app.baseSpan * Math.pow(10, 2 * app.trail - 1) * (stage.braidTarget ? 2.5 : 1);
}

function onSample(t, s) {
  stage.pushTrail(t, s);
  const area = signedArea(s, app.normal);
  const sign = Math.sign(area);
  if (sign && app.areaSign && sign !== app.areaSign) {
    const bead = document.createElement('i');
    bead.className = `b${middleBody(s, app.normal) + 1}`;
    const tape = $('tape-beads');
    tape.append(bead);
    if (tape.childElementCount > TAPE_LENGTH) tape.firstElementChild.remove();
  }
  if (sign) app.areaSign = sign;
}

const ghostBuf = new Float64Array(9);

// Positions straight from the fitted formula: a sum of circular motions.
function formulaPositions(formula, t) {
  ghostBuf.fill(0);
  for (let k = 0; k < formula.freqs.length; k++) {
    const c = Math.cos(formula.freqs[k] * t), s = Math.sin(formula.freqs[k] * t);
    for (let b = 0; b < 3; b++) {
      const [re, im] = formula.amps[b][k];
      ghostBuf[3 * b] += re * c - im * s;
      ghostBuf[3 * b + 1] += re * s + im * c;
    }
  }
  return ghostBuf;
}

function setPlaying(on) {
  app.playing = on;
  $('play').textContent = on ? 'Pause' : 'Play';
  $('play').setAttribute('aria-label', on ? 'Pause' : 'Play');
}

const SUPER = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
function sci(x) {
  if (x === 0) return '0';
  const [mantissa, exp] = x.toExponential(1).split('e');
  return `${mantissa} × 10${[...String(+exp)].map((ch) => SUPER[ch]).join('')}`;
}

function updateReadout() {
  const { sim, system } = app;
  $('r-time').textContent = sim.t.toFixed(2);
  $('r-energy').textContent = sci(Math.abs((energy(sim.s, system.masses) - app.e0) / app.e0));
  $('r-sep').textContent = minSeparation(sim.s).toFixed(4);
  if (system.formula && app.showFormula) {
    const g = formulaPositions(system.formula, sim.t);
    let miss = 0;
    for (let b = 0; b < 3; b++) miss = Math.max(miss, Math.hypot(g[3 * b] - sim.s[3 * b], g[3 * b + 1] - sim.s[3 * b + 1]));
    $('r-formula').textContent = sci(miss);
  }
}

let last = performance.now();
function tick(nowMs) {
  const dt = Math.min((nowMs - last) / 1000, 0.05);
  last = nowMs;
  const { sim, system } = app;
  if (sim) {
    if (app.playing) {
      sim.advanceTo(sim.t + dt * app.speed, app.maxGap, onSample, 6000);
      if (app.gone < 0) {
        app.gone = escaper(sim.s, system.masses);
        if (app.gone >= 0) {
          $('r-fate').textContent = `Body ${app.gone + 1} is ejected at time ${sim.t.toFixed(1)}. The other two remain bound as a binary.`;
        }
      } else if (Math.hypot(sim.s[3 * app.gone], sim.s[3 * app.gone + 1], sim.s[3 * app.gone + 2]) > 60 * stage.scale) {
        setPlaying(false);
      }
    }
    stage.setBodies(sim.s);
    stage.setGhosts(system.formula && app.showFormula ? formulaPositions(system.formula, sim.t) : null);
    if (app.frame++ % 5 === 0) updateReadout();
    stage.render(sim.t);
  }
  requestAnimationFrame(tick);
}

// controls

$('play').addEventListener('click', () => setPlaying(!app.playing));
$('restart').addEventListener('click', () => loadSystem(app.system));
$('speed').addEventListener('input', (e) => {
  app.speed = Math.pow(10, +e.target.value);
  $('speed-out').textContent = `${app.speed < 10 ? app.speed.toFixed(1) : Math.round(app.speed)}×`;
});
$('trail').addEventListener('input', (e) => {
  app.trail = +e.target.value;
  applyTrail();
});
$('braid').addEventListener('change', (e) => {
  stage.braidTarget = e.target.checked ? 1 : 0;
  applyTrail();
});
$('formula').addEventListener('change', (e) => {
  app.showFormula = e.target.checked;
  $('r-formula-row').hidden = !(app.system.formula && app.showFormula);
});
$('spin').addEventListener('change', (e) => (stage.controls.autoRotate = e.target.checked));

window.addEventListener('keydown', (e) => {
  if (e.target.matches('input, button')) return;
  if (e.code === 'Space') {
    e.preventDefault();
    setPlaying(!app.playing);
  } else if (e.key === 'r') {
    loadSystem(app.system);
  } else if (e.key === 'b') {
    $('braid').click();
  }
});

// drawer

const tabs = [...document.querySelectorAll('.tabs button')];
function openTab(name) {
  const current = tabs.find((t) => t.getAttribute('aria-pressed') === 'true')?.dataset.tab;
  const next = current === name ? null : name;
  tabs.forEach((t) => t.setAttribute('aria-pressed', String(t.dataset.tab === next)));
  $('drawer').hidden = !next;
  document.body.classList.toggle('drawer-open', !!next);
  document.querySelectorAll('[data-panel]').forEach((p) => (p.hidden = p.dataset.panel !== next));
  stage.setInset(next && window.innerWidth > 1100 ? $('drawer').offsetWidth : 0);
}
tabs.forEach((t) => t.addEventListener('click', () => openTab(t.dataset.tab)));

// data

let mapPlate = null;
let orbits = [];

function orbitSystem(o) {
  let note;
  if (!o.period) {
    const f = o.formulaStats;
    note = `Quasi-periodic: it never repeats. A ${o.formula.freqs.length}-term formula fitted to the first ${f.fitSpan} time units stays within ${sci(f.missAt10x)} of the true motion out to time ${f.testedTo}. Turn on "Show formula" to compare.`;
  } else {
    const stability = o.stable
      ? 'Linearly stable: small disturbances stay small.'
      : `Unstable: a small disturbance grows ${o.growth < 1000 ? o.growth.toFixed(1) : o.growth.toExponential(0)}-fold per period.`;
    note = `Period ${o.period.toFixed(3)}. ${stability}`;
    if (o.turn) note += ` Family member: the whole pattern rotates by ${o.turn} of a half-turn every half cycle.`;
  }
  return {
    id: o.id,
    name: o.name,
    note,
    state: o.state,
    masses: o.masses || EQUAL,
    period: o.period,
    formula: o.formula || null,
    mapPoint: o.kind === 'flat' || o.kind === 'orderly' ? o.params : null,
  };
}

function describeMapPoint(p1, p2, cell) {
  if (!cell) return '';
  if (cell.escaper >= 0) return `In the survey, body ${cell.escaper + 1} is ejected around time ${cell.tEnd.toFixed(0)}.`;
  return cell.orderly
    ? 'The survey classifies this system as regular: bound and predictable.'
    : 'The survey classifies this system as chaotic: still bound when the survey ended.';
}

async function getJson(url) {
  try {
    const r = await fetch(url);
    return r.ok ? await r.json() : null;
  } catch {
    return null;
  }
}

async function init() {
  $('findings').innerHTML = findingsHtml;
  $('findings').addEventListener('click', (e) => {
    const name = e.target.closest('[data-orbit]')?.dataset.orbit;
    const o = name && orbits.find((x) => x.name === name);
    if (o) loadSystem(orbitSystem(o));
  });

  const [orbitData, mapMeta] = await Promise.all([getJson('/data/orbits.json'), getJson('/data/map.json')]);
  orbits = orbitData || [];
  buildOrbitList($('orbit-list'), orbits, BODY_COLORS, (o) => loadSystem(orbitSystem(o)));

  if (mapMeta) {
    mapPlate = new MapPlate($('map-canvas'), $('map-readout'), mapMeta, orbits.filter((o) => o.kind === 'flat'));
    await mapPlate.load();
    mapPlate.onPick = (p1, p2, cell) =>
      loadSystem({
        name: `Initial velocity (${p1.toFixed(3)}, ${p2.toFixed(3)})`,
        note: describeMapPoint(p1, p2, cell),
        state: collinearStart(p1, p2),
        mapPoint: [p1, p2],
      });
  }

  const first = orbits.find((o) => o.featured) || orbits[0];
  loadSystem(
    first
      ? orbitSystem(first)
      : { name: 'Initial velocity (0.347, 0.533)', note: 'No research data loaded yet.', state: collinearStart(0.347, 0.533) },
  );
  requestAnimationFrame(tick);
}

init();
