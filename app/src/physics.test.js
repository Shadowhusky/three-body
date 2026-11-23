import assert from 'node:assert/strict';
import { ThreeBody, energy } from './physics.js';

const m = [1, 1, 1];
const start = (p1, p2) => [-1, 0, 0, 1, 0, 0, 0, 0, 0, p1, p2, 0, p1, p2, 0, -2 * p1, -2 * p2, 0];

{
  const sim = new ThreeBody(m);
  const s0 = start(0.2, 0.31);
  sim.setState(s0);
  const e0 = energy(s0, m);
  let samples = 0;
  sim.advanceTo(40, 0.01, () => samples++, 1e7);
  assert.ok(samples >= 4000);
  assert.ok(Math.abs(energy(sim.s, m) - e0) / Math.abs(e0) < 1e-7, 'energy drifts on a violent run');
}

{
  // Reference state at t = 50 from the C integrator (order 30, tol 1e-19).
  const ref = [
    -1.0312301082774806, -0.3142942566332039, 0, 0.4753769092109488, -0.1956135773885795, 0,
    0.5558531990666115, 0.5099078340217735, 0, -0.1958895887258697, 0.4329335982153318, 0,
    1.1808435383207878, -0.19190615892204738, 0, -0.9849539495949153, -0.24102743929328496, 0,
  ];
  {
    const sim = new ThreeBody(m, 24, 1e-17);
    sim.setState(start(0.347, 0.533));
    sim.advanceTo(50, 1, null, 1e7);
    const err = Math.hypot(...ref.map((x, i) => x - sim.s[i]));
    assert.ok(err < 1e-10, `disagrees with the C integrator by ${err}`);
    console.log('matches C integrator to', err.toExponential(2));
  }
}
console.log('physics ok');
