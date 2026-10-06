// DecayAnalysis.cpp
#include "DecayAnalysis.h"

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <numeric>
#include <utility>

namespace {

const double SQRT2 = 1.4142135623730951;
const double BIG   = 1e300;

// A step at t0, smoothed by a Gaussian of width sigma.
double stepTerm(double t, double t0, double sigma)
{
    return 0.5 * std::erfc(-(t - t0) / (sigma * SQRT2));
}

// exp(-(t - t0) / tau) starting at t0, smoothed by a Gaussian of width sigma.
double expTerm(double t, double t0, double sigma, double tau)
{
    double s = t - t0;
    double z = (sigma / tau - s / sigma) / SQRT2;
    if (z > 25.0) return 0.0;
    double e = sigma * sigma / (2.0 * tau * tau) - s / tau;
    if (e > 700.0) return 0.0;
    return 0.5 * std::exp(e) * std::erfc(z);
}

// q = [t0, sigma, c, A_1, tau_1, A_2, tau_2, ...]
double modelValue(double t, const std::vector<double>& q, int n)
{
    double v = q[2] * stepTerm(t, q[0], q[1]);
    for (int i = 0; i < n; i++) v += q[3 + 2 * i] * expTerm(t, q[0], q[1], q[4 + 2 * i]);
    return v;
}

// Gauss-Jordan inverse of an n x n matrix (row major). False if singular.
bool invert(std::vector<double> a, int n, std::vector<double>& inv)
{
    inv.assign(n * n, 0.0);
    for (int i = 0; i < n; i++) inv[i * n + i] = 1.0;
    for (int col = 0; col < n; col++) {
        int piv = col;
        double best = std::fabs(a[col * n + col]);
        for (int r = col + 1; r < n; r++) {
            if (std::fabs(a[r * n + col]) > best) { best = std::fabs(a[r * n + col]); piv = r; }
        }
        if (!(best > 1e-300)) return false;
        if (piv != col) {
            for (int c = 0; c < n; c++) {
                std::swap(a[piv * n + c], a[col * n + c]);
                std::swap(inv[piv * n + c], inv[col * n + c]);
            }
        }
        double d = a[col * n + col];
        for (int c = 0; c < n; c++) { a[col * n + c] /= d; inv[col * n + c] /= d; }
        for (int r = 0; r < n; r++) {
            if (r == col) continue;
            double f = a[r * n + col];
            if (f == 0.0) continue;
            for (int c = 0; c < n; c++) {
                a[r * n + c]   -= f * a[col * n + c];
                inv[r * n + c] -= f * inv[col * n + c];
            }
        }
    }
    return true;
}

// For fixed t0, sigma and decay times the amplitudes (and the offset) enter
// linearly, so they are found exactly by weighted least squares. theta =
// [t0, ln sigma, ln tau_1, ...] if fitRise, otherwise [ln tau_1, ...] (and t0 = 0,
// sigma = 1 ns are fixed). Returns chi2 (BIG if the fit is impossible).
// beta, if given, is [c (if withOffset), A_1, ..., A_n].
double chi2Projected(const std::vector<double>& theta, int n, bool withOffset, bool fitRise,
                     const std::vector<double>& t, const std::vector<double>& y,
                     const std::vector<double>& w, std::vector<double>* beta)
{
    double t0 = 0.0, sigma = 1.0;
    const int o = fitRise ? 2 : 0;
    if (fitRise) {
        t0 = theta[0];
        sigma = std::exp(theta[1]);
        if (!(sigma > 0.05 && sigma < 200.0)) return BIG;
    }
    std::vector<double> tau(n);
    for (int i = 0; i < n; i++) {
        tau[i] = std::exp(theta[o + i]);
        if (!(tau[i] > 0.5 && tau[i] < 1e6)) return BIG;
    }

    const int first = withOffset ? 1 : 0;
    const int m = n + first;
    std::vector<double> ata(m * m, 0.0), atb(m, 0.0), col(m);
    for (size_t k = 0; k < t.size(); k++) {
        if (withOffset) col[0] = stepTerm(t[k], t0, sigma);
        for (int i = 0; i < n; i++) col[first + i] = expTerm(t[k], t0, sigma, tau[i]);
        for (int a = 0; a < m; a++) {
            atb[a] += w[k] * col[a] * y[k];
            for (int b = 0; b <= a; b++) ata[a * m + b] += w[k] * col[a] * col[b];
        }
    }
    for (int a = 0; a < m; a++) {
        for (int b = a + 1; b < m; b++) ata[a * m + b] = ata[b * m + a];
    }

    // Scale the columns so that the matrix is well conditioned.
    std::vector<double> d(m);
    for (int a = 0; a < m; a++) {
        if (!(ata[a * m + a] > 0.0)) return BIG;
        d[a] = std::sqrt(ata[a * m + a]);
    }
    std::vector<double> scaled(m * m);
    for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) scaled[a * m + b] = ata[a * m + b] / (d[a] * d[b]);
    }
    std::vector<double> inv;
    if (!invert(scaled, m, inv)) return BIG;

    std::vector<double> sol(m, 0.0);
    for (int a = 0; a < m; a++) {
        double s = 0.0;
        for (int b = 0; b < m; b++) s += inv[a * m + b] * atb[b] / d[b];
        sol[a] = s / d[a];
    }

    double chi2 = 0.0;
    for (size_t k = 0; k < t.size(); k++) {
        double v = 0.0;
        if (withOffset) v += sol[0] * stepTerm(t[k], t0, sigma);
        for (int i = 0; i < n; i++) v += sol[first + i] * expTerm(t[k], t0, sigma, tau[i]);
        double e = y[k] - v;
        chi2 += w[k] * e * e;
    }
    if (beta) *beta = sol;
    return chi2;
}

// Nelder-Mead minimisation (no derivatives needed).
std::vector<double> nelderMead(const std::function<double(const std::vector<double>&)>& f,
                               const std::vector<double>& start, const std::vector<double>& step,
                               int maxEval, double tol, double& fbest)
{
    const int n = (int)start.size();
    std::vector<std::vector<double> > pts(n + 1, start);
    for (int i = 0; i < n; i++) pts[i + 1][i] += step[i];
    std::vector<double> fv(n + 1);
    int evals = 0;
    for (int i = 0; i <= n; i++) { fv[i] = f(pts[i]); evals++; }

    while (evals < maxEval) {
        std::vector<int> idx(n + 1);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return fv[a] < fv[b]; });
        std::vector<std::vector<double> > p2(n + 1);
        std::vector<double> f2(n + 1);
        for (int i = 0; i <= n; i++) { p2[i] = pts[idx[i]]; f2[i] = fv[idx[i]]; }
        pts.swap(p2);
        fv.swap(f2);

        if (std::fabs(fv[n] - fv[0]) <= tol * (std::fabs(fv[0]) + tol)) break;

        std::vector<double> centroid(n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) centroid[j] += pts[i][j] / n;
        }
        // The point centroid + coef * (worst - centroid): coef -1 reflects, -2 expands,
        // -0.5 contracts outside, +0.5 contracts inside.
        auto along = [&](double coef) {
            std::vector<double> p(n);
            for (int j = 0; j < n; j++) p[j] = centroid[j] + coef * (pts[n][j] - centroid[j]);
            return p;
        };

        std::vector<double> xr = along(-1.0);
        double fr = f(xr);
        evals++;
        if (fr < fv[0]) {
            std::vector<double> xe = along(-2.0);
            double fe = f(xe);
            evals++;
            if (fe < fr) { pts[n] = xe; fv[n] = fe; }
            else         { pts[n] = xr; fv[n] = fr; }
        } else if (fr < fv[n - 1]) {
            pts[n] = xr;
            fv[n] = fr;
        } else {
            std::vector<double> xc = (fr < fv[n]) ? along(-0.5) : along(0.5);
            double fc = f(xc);
            evals++;
            if (fc < std::min(fr, fv[n])) {
                pts[n] = xc;
                fv[n] = fc;
            } else {
                for (int i = 1; i <= n; i++) {
                    for (int j = 0; j < n; j++) pts[i][j] = pts[0][j] + 0.5 * (pts[i][j] - pts[0][j]);
                    fv[i] = f(pts[i]);
                    evals++;
                }
            }
        }
    }

    int best = 0;
    for (int i = 1; i <= n; i++) if (fv[i] < fv[best]) best = i;
    fbest = fv[best];
    return pts[best];
}

} // namespace

DecayAnalysis::DecayAnalysis(const Settings& settings)
    : cfg(settings), events(0) {}

int DecayAnalysis::preSamples() const  { return (int)std::lround(cfg.preNs / cfg.nsPerSample); }
int DecayAnalysis::postSamples() const { return (int)std::lround(cfg.postNs / cfg.nsPerSample); }

bool DecayAnalysis::add(const std::vector<double>& x, int onsetSample)
{
    const int n = (int)x.size();
    const int pre = preSamples();
    const int post = postSamples();
    const int areaN = (int)std::lround(cfg.areaNs / cfg.nsPerSample);
    if (sumY.empty()) {
        sumY.assign(pre + post, 0.0);
        sumYY.assign(pre + post, 0.0);
    }
    if (onsetSample < 1 || onsetSample >= n) return false;

    // The level before the pulse, and the height of the pulse.
    int lo = std::max(0, onsetSample - 25);
    int hi = std::max(lo + 1, onsetSample - 5);
    double before = 0.0;
    for (int j = lo; j < hi; j++) before += x[j];
    before /= (hi - lo);

    int end = std::min(n - 1, onsetSample + cfg.peakSamples);
    double peak = x[onsetSample];
    for (int j = onsetSample; j <= end; j++) peak = std::max(peak, x[j]);
    double height = peak - before;
    if (!(height > 0.0)) return false;
    double level = before + cfg.alignFraction * height;

    // The first sample at or above the alignment level.
    int a = -1;
    int jEnd = std::min(n - 1, onsetSample + cfg.searchAfter);
    for (int j = std::max(1, onsetSample - cfg.searchBefore); j <= jEnd; j++) {
        if (x[j] >= level) { a = j; break; }
    }
    if (a < 0) return false;
    if (a - pre < 0 || a + post > n || a + areaN > n) return false;

    double area = 0.0;
    for (int k = 0; k < areaN; k++) area += x[a + k];
    area *= cfg.nsPerSample;
    if (!(area > 0.0)) return false;

    for (int k = 0; k < pre + post; k++) {
        double y = x[a - pre + k] / area;
        sumY[k]  += y;
        sumYY[k] += y * y;
    }
    events++;
    return true;
}

DecayAnalysis::Average DecayAnalysis::average() const
{
    Average a;
    a.events = events;
    if (events < 2) return a;

    const int pre = preSamples();
    const double N = (double)events;
    for (size_t k = 0; k < sumY.size(); k++) {
        double mean = sumY[k] / N;
        double var = (sumYY[k] / N - mean * mean) * N / (N - 1.0);
        if (var < 0.0) var = 0.0;
        a.t.push_back(((int)k - pre) * cfg.nsPerSample);
        a.mean.push_back(mean);
        a.error.push_back(std::sqrt(var / N));
    }
    return a;
}

double DecayAnalysis::peakTimeNs() const
{
    Average a = average();
    double best = 0.0, top = -1e300;
    for (size_t k = 0; k < a.t.size(); k++) {
        if (a.mean[k] > top) { top = a.mean[k]; best = a.t[k]; }
    }
    return best;
}

DecayAnalysis::FitResult DecayAnalysis::fit(int nComponents, bool withOffset, double startNs, bool fitRise) const
{
    Average a = average();
    return fitData(a.t, a.mean, a.error, nComponents, withOffset, startNs, fitRise, events);
}

DecayAnalysis::FitResult DecayAnalysis::fitData(const std::vector<double>& t, const std::vector<double>& y,
                                                const std::vector<double>& err, int N, bool withOffset,
                                                double startNs, bool fitRise, long eventCount) const
{
    FitResult r;
    r.ok = false;
    r.fitRise = fitRise;
    r.fitStartNs = startNs;
    r.fitEndNs = cfg.fitEndNs;
    r.errorsOk = false;
    r.withOffset = withOffset;
    r.nComponents = N;
    r.events = eventCount;
    r.t0 = r.t0Err = r.sigma = r.sigmaErr = r.rise1090 = r.rise1090Err = 0.0;
    r.offset = r.offsetErr = r.offsetFraction = r.offsetFractionErr = r.offsetShare = r.offsetShareErr = 0.0;
    r.shareSum = 0.0;
    r.chi2 = 0.0;
    r.ndf = 0;
    r.errorScale = 1.0;
    r.allPositive = true;

    if (N < 1 || N > 3 || t.size() != y.size() || t.size() != err.size()) return r;

    // The points to fit. Each error is increased by a small systematic floor, so
    // that no single very precise point dominates the fit.
    double ymax = 0.0;
    for (double v : y) ymax = std::max(ymax, v);
    const double floorErr = cfg.systematicFloor * ymax;
    std::vector<double> ft, fy, fe, fw;
    for (size_t k = 0; k < t.size(); k++) {
        if (t[k] < startNs || t[k] > cfg.fitEndNs) continue;
        double e = std::sqrt(err[k] * err[k] + floorErr * floorErr);
        if (!(e > 0.0)) continue;
        ft.push_back(t[k]);
        fy.push_back(y[k]);
        fe.push_back(e);
        fw.push_back(1.0 / (e * e));
    }
    const int nFreeNeeded = (fitRise ? 2 : 0) + (withOffset ? 1 : 0) + 2 * N;
    if ((int)ft.size() < nFreeNeeded + 3) return r;

    // Every combination of starting values is tried.
    std::vector<const std::vector<double>*> lists;
    lists.push_back(&cfg.fastStarts);
    lists.push_back(&cfg.mediumStarts);
    lists.push_back(&cfg.slowStarts);
    std::vector<std::vector<double> > tauSets(1);
    for (int i = 0; i < N; i++) {
        std::vector<std::vector<double> > next;
        for (size_t s = 0; s < tauSets.size(); s++) {
            for (double v : *lists[i]) {
                std::vector<double> s2 = tauSets[s];
                s2.push_back(v);
                next.push_back(s2);
            }
        }
        tauSets.swap(next);
    }
    if (tauSets.empty() || cfg.sigmaStarts.empty()) return r;

    auto objective = [&](const std::vector<double>& th) {
        return chi2Projected(th, N, withOffset, fitRise, ft, fy, fw, nullptr);
    };

    double bestChi2 = BIG;
    std::vector<double> bestTheta;
    std::vector<double> sigmaList = cfg.sigmaStarts;
    if (!fitRise) sigmaList.assign(1, 1.0);
    for (double s0 : sigmaList) {
        for (size_t s = 0; s < tauSets.size(); s++) {
            const int o = fitRise ? 2 : 0;
            std::vector<double> start(o + N), step(o + N);
            if (fitRise) {
                start[0] = 0.0;           step[0] = 2.0;
                start[1] = std::log(s0);  step[1] = 0.3;
            }
            for (int i = 0; i < N; i++) {
                start[o + i] = std::log(tauSets[s][i]);
                step[o + i] = 0.3;
            }
            double fb = BIG;
            std::vector<double> th = nelderMead(objective, start, step, cfg.maxEvaluations, cfg.tolerance, fb);
            std::vector<double> step2 = step;
            for (size_t j = 0; j < step2.size(); j++) step2[j] *= 0.3;
            th = nelderMead(objective, th, step2, cfg.maxEvaluations, cfg.tolerance, fb);
            if (fb < bestChi2) { bestChi2 = fb; bestTheta = th; }
        }
    }
    if (bestTheta.empty() || !(bestChi2 < BIG)) return r;

    // Put the components in order of decay time and get the amplitudes.
    const int o = fitRise ? 2 : 0;
    std::vector<double> taus(N);
    for (int i = 0; i < N; i++) taus[i] = bestTheta[o + i];
    std::sort(taus.begin(), taus.end());
    std::vector<double> theta = bestTheta;
    for (int i = 0; i < N; i++) theta[o + i] = taus[i];
    std::vector<double> beta;
    const double chi2 = chi2Projected(theta, N, withOffset, fitRise, ft, fy, fw, &beta);
    if (!(chi2 < BIG)) return r;

    const int P = 3 + 2 * N;                       // t0, sigma, c, then A and tau for each component
    std::vector<double> q(P, 0.0);
    q[0] = fitRise ? theta[0] : 0.0;
    q[1] = fitRise ? std::exp(theta[1]) : 1.0;
    q[2] = withOffset ? beta[0] : 0.0;
    const int first = withOffset ? 1 : 0;
    for (int i = 0; i < N; i++) {
        q[3 + 2 * i] = beta[first + i];
        q[4 + 2 * i] = std::exp(theta[o + i]);
    }

    std::vector<int> freeIdx;
    for (int p = 0; p < P; p++) {
        if ((p == 0 || p == 1) && !fitRise) continue;
        if (p == 2 && !withOffset) continue;
        freeIdx.push_back(p);
    }
    const int nFree = (int)freeIdx.size();

    // The size of each parameter, used to scale the derivatives.
    const double ampScale = 1e-3 * ymax;
    auto scaleOf = [&](int p) {
        if (p == 0) return 1.0;
        if (p == 1) return q[1];
        if (p == 2 || (p >= 3 && (p - 3) % 2 == 0)) return std::max(std::fabs(q[p]), ampScale);
        return q[p];
    };

    // Covariance from the derivatives of the model with respect to the parameters.
    const int nPts = (int)ft.size();
    std::vector<std::vector<double> > J(nFree, std::vector<double>(nPts, 0.0));
    for (int a = 0; a < nFree; a++) {
        int p = freeIdx[a];
        double s = scaleOf(p);
        double h = 1e-4 * s;
        std::vector<double> qp = q, qm = q;
        qp[p] += h;
        qm[p] -= h;
        for (int k = 0; k < nPts; k++) {
            J[a][k] = (modelValue(ft[k], qp, N) - modelValue(ft[k], qm, N)) / (2.0 * h) * s;
        }
    }
    std::vector<double> fisher(nFree * nFree, 0.0);
    for (int a = 0; a < nFree; a++) {
        for (int b = 0; b < nFree; b++) {
            double s = 0.0;
            for (int k = 0; k < nPts; k++) s += fw[k] * J[a][k] * J[b][k];
            fisher[a * nFree + b] = s;
        }
    }
    std::vector<double> covScaled;
    bool covOk = invert(fisher, nFree, covScaled);
    std::vector<double> cov(P * P, 0.0);
    if (covOk) {
        for (int a = 0; a < nFree; a++) {
            for (int b = 0; b < nFree; b++) {
                cov[freeIdx[a] * P + freeIdx[b]] = covScaled[a * nFree + b] * scaleOf(freeIdx[a]) * scaleOf(freeIdx[b]);
            }
        }
        for (int a = 0; a < nFree; a++) {
            if (!(cov[freeIdx[a] * P + freeIdx[a]] >= 0.0)) covOk = false;
        }
    }

    const int ndf = nPts - nFree;
    const double errorScale = (ndf > 0 && chi2 / ndf > 1.0) ? std::sqrt(chi2 / ndf) : 1.0;

    auto paramErr = [&](int p) {
        double v = covOk ? cov[p * P + p] : 0.0;
        return std::sqrt(std::max(v, 0.0)) * errorScale;
    };

    // Quantities derived from the parameters: weights and shares of the light.
    const double T = cfg.areaNs;
    auto derived = [&](const std::vector<double>& qq) {
        std::vector<double> d(2 * N + 2, 0.0);
        double sumA = qq[2];
        for (int i = 0; i < N; i++) sumA += qq[3 + 2 * i];
        for (int i = 0; i < N; i++) {
            double A = qq[3 + 2 * i], tau = qq[4 + 2 * i];
            d[2 * i]     = (sumA != 0.0) ? A / sumA : 0.0;
            d[2 * i + 1] = A * tau * (1.0 - std::exp(-T / tau));
        }
        d[2 * N]     = (sumA != 0.0) ? qq[2] / sumA : 0.0;
        d[2 * N + 1] = qq[2] * T;
        return d;
    };
    std::vector<double> dval = derived(q);
    std::vector<double> derr(dval.size(), 0.0);
    if (covOk) {
        for (size_t di = 0; di < dval.size(); di++) {
            std::vector<double> grad(P, 0.0);
            for (int a = 0; a < nFree; a++) {
                int p = freeIdx[a];
                double h = 1e-4 * scaleOf(p);
                std::vector<double> qp = q, qm = q;
                qp[p] += h;
                qm[p] -= h;
                grad[p] = (derived(qp)[di] - derived(qm)[di]) / (2.0 * h);
            }
            double var = 0.0;
            for (int a = 0; a < nFree; a++) {
                for (int b = 0; b < nFree; b++) {
                    var += grad[freeIdx[a]] * grad[freeIdx[b]] * cov[freeIdx[a] * P + freeIdx[b]];
                }
            }
            derr[di] = std::sqrt(std::max(var, 0.0)) * errorScale;
        }
    }

    // Assemble the result.
    r.ok = true;
    r.errorsOk = covOk;
    if (fitRise) {
        r.t0 = q[0];            r.t0Err = paramErr(0);
        r.sigma = q[1];         r.sigmaErr = paramErr(1);
        r.rise1090 = 2.5631 * q[1];
        r.rise1090Err = 2.5631 * r.sigmaErr;
    }
    r.offset = q[2];        r.offsetErr = withOffset ? paramErr(2) : 0.0;
    r.offsetFraction = dval[2 * N];      r.offsetFractionErr = derr[2 * N];
    r.offsetShare = dval[2 * N + 1];     r.offsetShareErr = derr[2 * N + 1];
    r.shareSum = r.offsetShare;
    for (int i = 0; i < N; i++) {
        Component c;
        c.amplitude = q[3 + 2 * i];   c.amplitudeErr = paramErr(3 + 2 * i);
        c.tau = q[4 + 2 * i];         c.tauErr = paramErr(4 + 2 * i);
        c.fraction = dval[2 * i];     c.fractionErr = derr[2 * i];
        c.share = dval[2 * i + 1];    c.shareErr = derr[2 * i + 1];
        r.shareSum += c.share;
        if (!(c.amplitude > 0.0)) r.allPositive = false;
        r.components.push_back(c);
    }
    r.chi2 = chi2;
    r.ndf = ndf;
    r.errorScale = errorScale;
    for (int k = 0; k < nPts; k++) {
        double m = modelValue(ft[k], q, N);
        r.t.push_back(ft[k]);
        r.y.push_back(fy[k]);
        r.err.push_back(fe[k]);
        r.model.push_back(m);
        r.pull.push_back((fy[k] - m) / fe[k]);
    }
    return r;
}

void DecayAnalysis::print(const FitResult& r, const std::string& label)
{
    std::printf("\n%s: %ld waveforms, %d components\n", label.c_str(), r.events, r.nComponents);
    if (!r.ok) {
        std::printf("  the fit could not be done\n");
        return;
    }
    std::printf("  fit from %.0f to %.0f ns, chi2/ndf = %.2f\n", r.fitStartNs, r.fitEndNs,
                r.ndf > 0 ? r.chi2 / r.ndf : 0.0);
    for (size_t i = 0; i < r.components.size(); i++) {
        const Component& c = r.components[i];
        std::printf("  tau%d = %8.1f +/- %6.1f ns    weight = %.3f +/- %.3f\n",
                    (int)i + 1, c.tau, c.tauErr, c.fraction, c.fractionErr);
    }
    if (r.withOffset) {
        std::printf("  offset weight = %.3f +/- %.3f\n", r.offsetFraction, r.offsetFractionErr);
    }
    if (!r.allPositive) std::printf("  WARNING: a component has a negative amplitude, so this fit is not physical\n");
    if (!r.errorsOk)    std::printf("  WARNING: the errors could not be computed\n");
}

void DecayAnalysis::writePdf(const std::string& pdfPath, const FitResult& r) const
{
    Average a = average();
    if (!r.ok || a.t.empty()) return;

    double maxMean = 0.0;
    for (size_t k = 0; k < a.mean.size(); k++) maxMean = std::max(maxMean, a.mean[k]);
    if (!(maxMean > 0.0)) return;

    // Time is counted from the start of the averaged window.
    std::vector<double> ta(a.t.size()), ya(a.t.size());
    for (size_t k = 0; k < a.t.size(); k++) {
        ta[k] = a.t[k] + cfg.preNs;
        ya[k] = a.mean[k] / maxMean;
    }
    std::vector<double> tm(r.t.size()), ym(r.t.size());
    for (size_t k = 0; k < r.t.size(); k++) {
        tm[k] = r.t[k] + cfg.preNs;
        ym[k] = r.model[k] / maxMean;
    }

    TCanvas* c = new TCanvas("cDecay", "", 800, 800);

    TGraph* g = new TGraph((int)ta.size(), ta.data(), ya.data());
    g->SetTitle(";Time (ns);Normalized Pulse Amplitude");
    g->SetLineColor(kGray + 1);
    g->SetLineWidth(2);
    g->SetMinimum(0.0);
    g->SetMaximum(1.05);
    g->Draw("AL");

    TGraph* m = new TGraph((int)tm.size(), tm.data(), ym.data());
    m->SetLineColor(kBlack);
    m->SetLineStyle(2);
    m->SetLineWidth(3);
    m->Draw("L");

    TLegend* leg = new TLegend(0.58, 0.78, 0.88, 0.88);
    leg->AddEntry(g, "Average pulse", "l");
    leg->AddEntry(m, "Decay fit", "l");
    leg->Draw();

    c->Print(pdfPath.c_str());

    delete leg;
    delete m;
    delete g;
    delete c;
}
