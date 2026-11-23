// Taylor-series integrator for the gravitational three-body problem (G = 1).
//
// At each step the solution is expanded as a power series in time. The
// coefficients follow from  x'' = sum m_j (x_j - x_i) / r^3  by term-by-term
// recurrence, so the order is arbitrary and no coefficient tables are needed.

#include <math.h>
#include <string.h>

#define MAXN 48
#define SYM_MAX 96

static const int PI[3] = {0, 0, 1}, PJ[3] = {1, 2, 2};

typedef struct {
    int ord;
    double X[3][3][MAXN + 2];
    double V[3][3][MAXN + 2];
} jet_t;

static void jet_compute(jet_t *J, const double *s, const double *m) {
    int N = J->ord;
    double D[3][3][MAXN + 2], S[3][MAXN + 2], U[3][MAXN + 2];
    for (int b = 0; b < 3; b++)
        for (int d = 0; d < 3; d++) {
            J->X[b][d][0] = s[3 * b + d];
            J->V[b][d][0] = s[9 + 3 * b + d];
        }
    for (int k = 0; k < N; k++) {
        double A[3][3] = {{0}};
        for (int p = 0; p < 3; p++) {
            int i = PI[p], j = PJ[p];
            for (int d = 0; d < 3; d++) D[p][d][k] = J->X[j][d][k] - J->X[i][d][k];
            double sk = 0;
            for (int l = 0; l <= k; l++)
                sk += D[p][0][l] * D[p][0][k - l] + D[p][1][l] * D[p][1][k - l] + D[p][2][l] * D[p][2][k - l];
            S[p][k] = sk;
            if (k == 0) {
                U[p][0] = 1.0 / (sk * sqrt(sk));
            } else {
                // u = s^(-3/2)  =>  s u' = -3/2 s' u
                double acc = 0;
                for (int q = 1; q <= k; q++) acc += (-1.5 * q - (k - q)) * S[p][q] * U[p][k - q];
                U[p][k] = acc / (k * S[p][0]);
            }
            for (int d = 0; d < 3; d++) {
                double f = 0;
                for (int l = 0; l <= k; l++) f += D[p][d][l] * U[p][k - l];
                A[i][d] += m[j] * f;
                A[j][d] -= m[i] * f;
            }
        }
        double inv = 1.0 / (k + 1);
        for (int b = 0; b < 3; b++)
            for (int d = 0; d < 3; d++) {
                J->X[b][d][k + 1] = J->V[b][d][k] * inv;
                J->V[b][d][k + 1] = A[b][d] * inv;
            }
    }
}

// Largest step for which the last two series terms stay below tol.
static double jet_step(const jet_t *J, double tol) {
    int N = J->ord;
    double h = 1e300;
    for (int k = N - 1; k <= N; k++) {
        double nrm = 0;
        for (int b = 0; b < 3; b++)
            for (int d = 0; d < 3; d++) {
                double a = fabs(J->X[b][d][k]), c = fabs(J->V[b][d][k]);
                if (a > nrm) nrm = a;
                if (c > nrm) nrm = c;
            }
        if (nrm > 0) {
            double hk = pow(tol / nrm, 1.0 / k);
            if (hk < h) h = hk;
        }
    }
    return 0.9 * h;
}

static void jet_eval(const jet_t *J, double tau, double *s) {
    int N = J->ord;
    for (int b = 0; b < 3; b++)
        for (int d = 0; d < 3; d++) {
            double x = J->X[b][d][N], v = J->V[b][d][N];
            for (int k = N - 1; k >= 0; k--) {
                x = x * tau + J->X[b][d][k];
                v = v * tau + J->V[b][d][k];
            }
            s[3 * b + d] = x;
            s[9 + 3 * b + d] = v;
        }
}

// Flow map: advance s0 by exactly t_end. Returns step count, or -1 if max_steps hit.
long tb_flow(const double *s0, const double *m, double t_end, double *s1, int ord, double tol, long max_steps) {
    jet_t J;
    J.ord = ord;
    double s[18], t = 0;
    memcpy(s, s0, sizeof s);
    double dir = t_end >= 0 ? 1.0 : -1.0, T = fabs(t_end);
    long n = 0;
    while (t < T) {
        if (n >= max_steps) return -1;
        jet_compute(&J, s, m);
        double h = jet_step(&J, tol);
        if (t + h > T) h = T - t;
        jet_eval(&J, dir * h, s);
        t += h;
        n++;
    }
    memcpy(s1, s, sizeof s);
    return n;
}

// Uniformly sampled trajectory: out[i] is the state at t = i * t_end / (n_out - 1).
long tb_sample(const double *s0, const double *m, double t_end, int n_out, double *out, int ord, double tol,
               long max_steps) {
    jet_t J;
    J.ord = ord;
    double s[18], t = 0, dt = t_end / (n_out - 1);
    memcpy(s, s0, sizeof s);
    memcpy(out, s, sizeof s);
    int next = 1;
    long n = 0;
    while (next < n_out) {
        if (n >= max_steps) return -1;
        jet_compute(&J, s, m);
        double h = jet_step(&J, tol);
        while (next < n_out && next * dt <= t + h) {
            jet_eval(&J, next * dt - t, out + 18 * next);
            next++;
        }
        jet_eval(&J, h, s);
        t += h;
        n++;
    }
    return n;
}

static double tri_area(const double *s) {
    return (s[3] - s[0]) * (s[7] - s[1]) - (s[4] - s[1]) * (s[6] - s[0]);
}

static int middle_body(const double *s) {
    for (int k = 0; k < 3; k++) {
        int i = (k + 1) % 3, j = (k + 2) % 3;
        double dot = 0;
        for (int d = 0; d < 3; d++) dot += (s[3 * i + d] - s[3 * k + d]) * (s[3 * j + d] - s[3 * k + d]);
        if (dot < 0) return k;
    }
    return 0;
}

static double state_dist(const double *a, const double *b) {
    double q = 0;
    for (int i = 0; i < 18; i++) q += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(q);
}

// Returns the escaping body's index if the system has split into a bound binary
// and a third body receding on a hyperbolic orbit, else -1.
static int escape_check(const double *s, const double *m, double r_esc) {
    int best = 0;
    double dmin = 1e300;
    for (int p = 0; p < 3; p++) {
        double q = 0;
        for (int d = 0; d < 3; d++) {
            double e = s[3 * PJ[p] + d] - s[3 * PI[p] + d];
            q += e * e;
        }
        if (q < dmin) {
            dmin = q;
            best = p;
        }
    }
    int i = PI[best], j = PJ[best], k = 3 - i - j;
    double mb = m[i] + m[j], r2 = 0, v2 = 0, rv = 0, w2 = 0;
    for (int d = 0; d < 3; d++) {
        double R = (m[i] * s[3 * i + d] + m[j] * s[3 * j + d]) / mb;
        double W = (m[i] * s[9 + 3 * i + d] + m[j] * s[9 + 3 * j + d]) / mb;
        double r = s[3 * k + d] - R, v = s[9 + 3 * k + d] - W, w = s[9 + 3 * j + d] - s[9 + 3 * i + d];
        r2 += r * r;
        v2 += v * v;
        rv += r * v;
        w2 += w * w;
    }
    double r = sqrt(r2), dij = sqrt(dmin);
    if (r < r_esc || rv <= 0) return -1;
    if (0.5 * v2 - (mb + m[k]) / r <= 0) return -1;
    double ebin = 0.5 * w2 - mb / dij;
    if (ebin >= 0) return -1;
    double a = -mb / (2 * ebin);
    if (r < 5 * a) return -1;
    return k;
}

// Diagnostics for one orbit.
//   out[0] status      0 survived to t_max, 1 broke up, 2 step budget exhausted
//   out[1] t_end       time reached
//   out[2] escaper     body index that left (or -1)
//   out[3] log_growth  ln of the total growth of a 1e-8 shadow-orbit separation
//   out[4] min_sep     closest approach of any pair
//   out[5] n_syz       number of syzygies (collinear configurations; planar motion only)
//   out[6] ret_min     closest phase-space return to the initial state after ret_tmin
//   out[7] ret_t       time of that closest return
//   out[8] steps
// sym: the middle body at each syzygy (first SYM_MAX events).
void tb_scan(const double *s0, const double *m, double t_max, int ord, double tol, double ret_tmin, double r_esc,
             long max_steps, double *out, unsigned char *sym) {
    const double d0 = 1e-8;
    const int NSUB = 6;
    jet_t J, K;
    J.ord = K.ord = ord;
    double s[18], sh[18], sub[18], prev[18], t = 0;
    memcpy(s, s0, sizeof s);
    memcpy(sh, s0, sizeof s);
    static const double dirv[12] = {0.31, -0.22, 0, -0.17, 0.41, 0, -0.14, -0.19, 0, 0.27, 0.13, 0};
    {
        double nrm = 0, pert[18] = {0};
        for (int b = 0; b < 2; b++)
            for (int d = 0; d < 2; d++) {
                pert[3 * b + d] = dirv[3 * b + d];
                pert[9 + 3 * b + d] = dirv[6 + 3 * b + d];
            }
        for (int i = 0; i < 18; i++) nrm += pert[i] * pert[i];
        nrm = sqrt(nrm);
        for (int i = 0; i < 18; i++) sh[i] += d0 * pert[i] / nrm;
    }
    double logsum = 0, min_sep = 1e300, ret_min = 1e300, ret_t = 0;
    int n_syz = 0, status = 0, escaper = -1;
    long n = 0;
    double area_prev = tri_area(s);
    memcpy(prev, s, sizeof s);
    while (t < t_max) {
        if (n >= max_steps) {
            status = 2;
            break;
        }
        jet_compute(&J, s, m);
        double h = jet_step(&J, tol);
        if (t + h > t_max) h = t_max - t;
        for (int q = 1; q <= NSUB; q++) {
            double tau = h * q / NSUB;
            jet_eval(&J, tau, sub);
            for (int p = 0; p < 3; p++) {
                double e2 = 0;
                for (int d = 0; d < 3; d++) {
                    double e = sub[3 * PJ[p] + d] - sub[3 * PI[p] + d];
                    e2 += e * e;
                }
                if (e2 < min_sep * min_sep) min_sep = sqrt(e2);
            }
            double area = tri_area(sub);
            if ((area > 0) != (area_prev > 0) && area != 0) {
                if (n_syz < SYM_MAX) {
                    double lo = h * (q - 1) / NSUB, hi = tau, mid[18];
                    double flo = area_prev;
                    for (int it = 0; it < 40; it++) {
                        double c = 0.5 * (lo + hi);
                        jet_eval(&J, c, mid);
                        double fc = tri_area(mid);
                        if ((fc > 0) == (flo > 0)) {
                            lo = c;
                            flo = fc;
                        } else
                            hi = c;
                    }
                    jet_eval(&J, 0.5 * (lo + hi), mid);
                    sym[n_syz] = (unsigned char)(middle_body(mid) + 1);
                }
                n_syz++;
            }
            if (area != 0) area_prev = area;
            if (t + tau > ret_tmin) {
                double dd = state_dist(sub, s0);
                if (dd < ret_min) {
                    ret_min = dd;
                    ret_t = t + tau;
                }
            }
        }
        jet_compute(&K, sh, m);
        jet_eval(&K, h, sh);
        jet_eval(&J, h, s);
        t += h;
        n++;
        double dd = state_dist(s, sh);
        if (dd > 1e-5 || dd < 1e-11) {
            logsum += log(dd / d0);
            for (int i = 0; i < 18; i++) sh[i] = s[i] + (sh[i] - s[i]) * d0 / dd;
        }
        escaper = escape_check(s, m, r_esc);
        if (escaper >= 0) {
            status = 1;
            break;
        }
    }
    double dd = state_dist(s, sh);
    out[0] = status;
    out[1] = t;
    out[2] = escaper;
    out[3] = logsum + log(dd / d0);
    out[4] = min_sep;
    out[5] = n_syz;
    out[6] = ret_min;
    out[7] = ret_t;
    out[8] = (double)n;
}

void tb_scan_batch(const double *s0s, int count, const double *m, double t_max, int ord, double tol, double ret_tmin,
                   double r_esc, long max_steps, double *outs, unsigned char *syms) {
    for (int i = 0; i < count; i++) {
        memset(syms + (size_t)i * SYM_MAX, 0, SYM_MAX);
        tb_scan(s0s + 18 * (size_t)i, m, t_max, ord, tol, ret_tmin, r_esc, max_steps, outs + 9 * (size_t)i,
                syms + (size_t)i * SYM_MAX);
    }
}

// Reversing symmetries of the equations of motion:
//   type 0: reflect through the origin, swap bodies 0 and 1, reverse velocities
//   type 1: half-turn about the y axis,  swap bodies 0 and 1, reverse velocities
// An orbit that starts from a state fixed by one of them and reaches another
// fixed state at time t is periodic with period 2t.
static void sym_apply(int type, const double *s, double *r) {
    static const int perm[3] = {1, 0, 2};
    double sx = -1, sy = type == 0 ? -1 : 1, sz = -1;
    for (int b = 0; b < 3; b++) {
        int c = perm[b];
        r[3 * b] = sx * s[3 * c];
        r[3 * b + 1] = sy * s[3 * c + 1];
        r[3 * b + 2] = sz * s[3 * c + 2];
        r[9 + 3 * b] = -sx * s[9 + 3 * c];
        r[9 + 3 * b + 1] = -sy * s[9 + 3 * c + 1];
        r[9 + 3 * b + 2] = -sz * s[9 + 3 * c + 2];
    }
}

static double sym_dist2(int type, const double *s) {
    double r[18], q = 0;
    sym_apply(type, s, r);
    for (int i = 0; i < 18; i++) q += (s[i] - r[i]) * (s[i] - r[i]);
    return q;
}

// Finds the times at which the orbit comes closest to a symmetric state.
// Keeps the kmax closest approaches below dthresh. Returns how many were kept.
int tb_symscan(const double *s0, const double *m, double t_max, int type, int ord, double tol, double t_min,
               double dthresh, double r_esc, long max_steps, int kmax, double *tmins, double *dmins) {
    const int NSUB = 6;
    jet_t J;
    J.ord = ord;
    double s[18], sub[18], t = 0;
    memcpy(s, s0, sizeof s);
    double ta = 0, tb_ = 0, qa = 0, qb = 0;
    int have = 0, kept = 0;
    long n = 0;
    while (t < t_max && n < max_steps) {
        jet_compute(&J, s, m);
        double h = jet_step(&J, tol);
        if (t + h > t_max) h = t_max - t;
        for (int q = 1; q <= NSUB; q++) {
            double tc = t + h * q / NSUB;
            jet_eval(&J, h * q / NSUB, sub);
            double qc = sym_dist2(type, sub);
            if (have >= 2 && qb < qa && qb <= qc && tb_ > t_min) {
                // vertex of the parabola through the three samples of squared distance
                double x1 = ta - tb_, x2 = tc - tb_, y1 = qa - qb, y2 = qc - qb;
                double den = x1 * x2 * (x1 - x2);
                double A = (x2 * y1 - x1 * y2) / den, B = (x1 * x1 * y2 - x2 * x2 * y1) / den;
                double tv = tb_, qv = qb;
                if (A > 0) {
                    tv = tb_ - B / (2 * A);
                    qv = qb - B * B / (4 * A);
                    if (qv < 0) qv = 0;
                }
                double dv = sqrt(qv);
                if (dv < dthresh) {
                    int pos = kept;
                    if (kept < kmax)
                        kept++;
                    else {
                        pos = -1;
                        double worst = dv;
                        for (int i = 0; i < kmax; i++)
                            if (dmins[i] > worst) {
                                worst = dmins[i];
                                pos = i;
                            }
                    }
                    if (pos >= 0) {
                        tmins[pos] = tv;
                        dmins[pos] = dv;
                    }
                }
            }
            ta = tb_;
            qa = qb;
            tb_ = tc;
            qb = qc;
            have++;
        }
        jet_eval(&J, h, s);
        t += h;
        n++;
        if (escape_check(s, m, r_esc) >= 0) break;
    }
    return kept;
}

void tb_symscan_batch(const double *s0s, int count, const double *m, double t_max, int type, int ord, double tol,
                      double t_min, double dthresh, double r_esc, long max_steps, int kmax, double *tmins,
                      double *dmins, int *kepts) {
    for (int i = 0; i < count; i++)
        kepts[i] = tb_symscan(s0s + 18 * (size_t)i, m, t_max, type, ord, tol, t_min, dthresh, r_esc, max_steps, kmax,
                              tmins + (size_t)i * kmax, dmins + (size_t)i * kmax);
}
