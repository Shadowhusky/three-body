// Same method as research/core/tb.c: the solution is expanded as a power series
// in time, with coefficients generated from the force law by recurrence (G = 1).

const PI = [0, 0, 1];
const PJ = [1, 2, 2];

export class ThreeBody {
  constructor(masses = [1, 1, 1], order = 18, tol = 1e-13) {
    this.m = Float64Array.from(masses);
    this.N = order;
    this.tol = tol;
    this.t = 0;
    this.s = new Float64Array(18);
    const n = order + 2;
    this.X = new Float64Array(9 * n);
    this.V = new Float64Array(9 * n);
    this.D = new Float64Array(9 * n);
    this.S = new Float64Array(3 * n);
    this.U = new Float64Array(3 * n);
    this.A = new Float64Array(9);
    this.n = n;
    this.h = 0;
    this.fresh = false;
  }

  setState(s, t = 0) {
    this.s.set(s);
    this.t = t;
    this.fresh = false;
  }

  computeSeries() {
    const { X, V, D, S, U, A, m, N, n, s } = this;
    for (let c = 0; c < 9; c++) {
      X[c * n] = s[c];
      V[c * n] = s[9 + c];
    }
    for (let k = 0; k < N; k++) {
      A.fill(0);
      for (let p = 0; p < 3; p++) {
        const i = PI[p], j = PJ[p];
        for (let d = 0; d < 3; d++) D[(3 * p + d) * n + k] = X[(3 * j + d) * n + k] - X[(3 * i + d) * n + k];
        let sk = 0;
        for (let l = 0; l <= k; l++) {
          for (let d = 0; d < 3; d++) sk += D[(3 * p + d) * n + l] * D[(3 * p + d) * n + k - l];
        }
        S[p * n + k] = sk;
        if (k === 0) {
          U[p * n] = 1 / (sk * Math.sqrt(sk));
        } else {
          let acc = 0;
          for (let q = 1; q <= k; q++) acc += (-1.5 * q - (k - q)) * S[p * n + q] * U[p * n + k - q];
          U[p * n + k] = acc / (k * S[p * n]);
        }
        for (let d = 0; d < 3; d++) {
          let f = 0;
          for (let l = 0; l <= k; l++) f += D[(3 * p + d) * n + l] * U[p * n + k - l];
          A[3 * i + d] += m[j] * f;
          A[3 * j + d] -= m[i] * f;
        }
      }
      const inv = 1 / (k + 1);
      for (let c = 0; c < 9; c++) {
        X[c * n + k + 1] = V[c * n + k] * inv;
        V[c * n + k + 1] = A[c] * inv;
      }
    }
    let h = Infinity;
    for (let k = N - 1; k <= N; k++) {
      let nrm = 0;
      for (let c = 0; c < 9; c++) {
        nrm = Math.max(nrm, Math.abs(X[c * n + k]), Math.abs(V[c * n + k]));
      }
      if (nrm > 0) h = Math.min(h, Math.pow(this.tol / nrm, 1 / k));
    }
    this.h = 0.9 * h;
    this.fresh = true;
  }

  evalAt(tau, out) {
    const { X, V, N, n } = this;
    for (let c = 0; c < 9; c++) {
      let x = X[c * n + N], v = V[c * n + N];
      for (let k = N - 1; k >= 0; k--) {
        x = x * tau + X[c * n + k];
        v = v * tau + V[c * n + k];
      }
      out[c] = x;
      out[9 + c] = v;
    }
  }

  // Advance to tTarget. onSample(t, state) fires along the way, at least twice per
  // series step and no further apart than maxGap, so tight encounters stay smooth.
  advanceTo(tTarget, maxGap, onSample, maxSteps = 20000) {
    const buf = (this._buf ||= new Float64Array(18));
    let steps = 0;
    while (this.t < tTarget && steps < maxSteps) {
      if (!this.fresh) this.computeSeries();
      const h = Math.min(this.h, tTarget - this.t);
      if (onSample) {
        const parts = Math.max(2, Math.ceil(h / maxGap));
        for (let q = 1; q < parts; q++) {
          const tau = (h * q) / parts;
          this.evalAt(tau, buf);
          onSample(this.t + tau, buf);
        }
      }
      this.evalAt(h, buf);
      this.s.set(buf);
      this.t += h;
      this.fresh = false;
      steps++;
      if (onSample) onSample(this.t, this.s);
    }
    return steps < maxSteps;
  }
}

export function energy(s, m) {
  let ke = 0, pe = 0;
  for (let b = 0; b < 3; b++) {
    ke += 0.5 * m[b] * (s[9 + 3 * b] ** 2 + s[10 + 3 * b] ** 2 + s[11 + 3 * b] ** 2);
  }
  for (let p = 0; p < 3; p++) {
    const i = PI[p], j = PJ[p];
    pe -= (m[i] * m[j]) / Math.hypot(s[3 * j] - s[3 * i], s[3 * j + 1] - s[3 * i + 1], s[3 * j + 2] - s[3 * i + 2]);
  }
  return ke + pe;
}

export function minSeparation(s) {
  let d = Infinity;
  for (let p = 0; p < 3; p++) {
    const i = PI[p], j = PJ[p];
    d = Math.min(d, Math.hypot(s[3 * j] - s[3 * i], s[3 * j + 1] - s[3 * i + 1], s[3 * j + 2] - s[3 * i + 2]));
  }
  return d;
}

// Signed area of the triangle projected along `normal`; its sign flips when the bodies line up.
export function signedArea(s, normal) {
  const ax = s[3] - s[0], ay = s[4] - s[1], az = s[5] - s[2];
  const bx = s[6] - s[0], by = s[7] - s[1], bz = s[8] - s[2];
  return (ay * bz - az * by) * normal[0] + (az * bx - ax * bz) * normal[1] + (ax * by - ay * bx) * normal[2];
}

// Which body sits between the other two, seen along `normal`.
export function middleBody(s, normal) {
  const q = new Float64Array(9);
  for (let b = 0; b < 3; b++) {
    const along = s[3 * b] * normal[0] + s[3 * b + 1] * normal[1] + s[3 * b + 2] * normal[2];
    for (let d = 0; d < 3; d++) q[3 * b + d] = s[3 * b + d] - along * normal[d];
  }
  for (let k = 0; k < 3; k++) {
    const i = (k + 1) % 3, j = (k + 2) % 3;
    let dot = 0;
    for (let d = 0; d < 3; d++) dot += (q[3 * i + d] - q[3 * k + d]) * (q[3 * j + d] - q[3 * k + d]);
    if (dot < 0) return k;
  }
  return 0;
}

// Body index that has left for good, or -1 (same test as escape_check in tb.c).
export function escaper(s, m, rEsc = 10) {
  let best = 0, dmin = Infinity;
  for (let p = 0; p < 3; p++) {
    const i = PI[p], j = PJ[p];
    const q = (s[3 * j] - s[3 * i]) ** 2 + (s[3 * j + 1] - s[3 * i + 1]) ** 2 + (s[3 * j + 2] - s[3 * i + 2]) ** 2;
    if (q < dmin) {
      dmin = q;
      best = p;
    }
  }
  const i = PI[best], j = PJ[best], k = 3 - i - j;
  const mb = m[i] + m[j];
  let r2 = 0, v2 = 0, rv = 0, w2 = 0;
  for (let d = 0; d < 3; d++) {
    const R = (m[i] * s[3 * i + d] + m[j] * s[3 * j + d]) / mb;
    const W = (m[i] * s[9 + 3 * i + d] + m[j] * s[9 + 3 * j + d]) / mb;
    const r = s[3 * k + d] - R, v = s[9 + 3 * k + d] - W, w = s[9 + 3 * j + d] - s[9 + 3 * i + d];
    r2 += r * r;
    v2 += v * v;
    rv += r * v;
    w2 += w * w;
  }
  const r = Math.sqrt(r2), dij = Math.sqrt(dmin);
  if (r < rEsc || rv <= 0) return -1;
  if (0.5 * v2 - (mb + m[k]) / r <= 0) return -1;
  const ebin = 0.5 * w2 - mb / dij;
  if (ebin >= 0) return -1;
  if (r < 5 * (-mb / (2 * ebin))) return -1;
  return k;
}
