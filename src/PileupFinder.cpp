// PileupFinder.cpp
#include "PileupFinder.h"

#include <algorithm>
#include <cmath>

double PileupFinder::noiseBefore(const std::vector<double>& x, int i) const
{
    int lo = std::max(0, i - cfg.noiseSamples);
    int m = i - lo;
    if (m < 8) return 0.0;

    // Least-squares straight line through the samples, then the scatter around it.
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int j = lo; j < i; j++) {
        double t = j - lo;
        sx += t;
        sy += x[j];
        sxx += t * t;
        sxy += t * x[j];
    }
    double den = m * sxx - sx * sx;
    double slope = (den != 0.0) ? (m * sxy - sx * sy) / den : 0.0;
    double offset = (sy - slope * sx) / m;

    double r2 = 0;
    for (int j = lo; j < i; j++) {
        double e = x[j] - (offset + slope * (j - lo));
        r2 += e * e;
    }
    return std::sqrt(r2 / m);
}

std::vector<PileupFinder::Pulse> PileupFinder::find(const std::vector<double>& x) const
{
    std::vector<Pulse> out;
    const int n = (int)x.size();
    const int w = cfg.windowSamples;
    if (w < 1 || n < 2 * w + 1) return out;

    const double eventPeak = *std::max_element(x.begin(), x.end());

    // Running sum, so the average over any stretch is cheap.
    std::vector<double> sum(n + 1, 0.0);
    for (int i = 0; i < n; i++) sum[i + 1] = sum[i] + x[i];

    bool inClimb = false;
    int lastPulse = -1000000;
    for (int i = w; i + w <= n; i++) {
        double before = (sum[i] - sum[i - w]) / w;
        double after  = (sum[i + w] - sum[i]) / w;

        if (after - before >= cfg.jump) {
            // First point of a new climb, not too soon after a pulse.
            if (!inClimb && i - lastPulse >= cfg.blankSamples) {
                int end = std::min(n - 1, i + cfg.peakSamples);
                double peak = x[i];
                for (int j = i; j <= end; j++) peak = std::max(peak, x[j]);
                double height = peak - before;

                double noise = noiseBefore(x, i);
                double needed = std::max(cfg.minFraction * eventPeak, cfg.noiseSigmas * noise);

                if (out.empty() || height >= needed) {
                    Pulse p;
                    p.onsetSample = i;
                    p.height = height;
                    p.noise = noise;
                    out.push_back(p);
                    lastPulse = i;
                    if ((int)out.size() >= cfg.maxPulses) break;
                }
            }
            inClimb = true;
        } else {
            inClimb = false;
        }
    }
    return out;
}
