// Waveform.cpp
#include "Waveform.h"
#include <cmath>

Waveform::Waveform(const int* wf, int nSamples, bool keepFullWaveform)
    : n(nSamples)
{
    double sigma_local = 0;
    baseline = computeBaselineSigmaCut(wf, n, SIGMA_CUT_A, SIGMA_CUT_MAX_ITER,
                                        SIGMA_CUT_TOL, &sigma_local);
    sigma = sigma_local;

    maxAmp = -1e18;

    if (keepFullWaveform) {
        subtracted.resize(n);
        for (int i = 0; i < n; i++) {
            subtracted[i] = baseline - wf[i];
            if (subtracted[i] > maxAmp) maxAmp = subtracted[i];
        }
    } else {
        for (int i = 0; i < n; i++) {
            double val = baseline - wf[i];
            if (val > maxAmp) maxAmp = val;
        }
    }
}

double Waveform::getBaseline() const { return baseline; }
double Waveform::getSigma()    const { return sigma; }

double Waveform::getMaxAmpADC()   const { return maxAmp; }
double Waveform::getMaxAmpVolts() const { return maxAmp * MV_PER_COUNT; }

const std::vector<double>& Waveform::getSubtractedADC() const { return subtracted; }

std::vector<double> Waveform::getSubtractedVolts() const {
    std::vector<double> mv(subtracted.size());
    for (size_t i = 0; i < subtracted.size(); i++) mv[i] = subtracted[i] * MV_PER_COUNT;
    return mv;
}

std::vector<double> Waveform::findPeaks(double threshold, int minSeparation)
{
    // STUB: not yet implemented.
    return std::vector<double>();
}

double Waveform::computeBaselineSigmaCut(const int* wf, int n,
                                          double A, int maxIter,
                                          double tol, double* sigmaOut)
{
    double sum = 0, sum2 = 0;
    for (int i = 0; i < n; i++) { sum += wf[i]; sum2 += (double)wf[i]*wf[i]; }
    double mean = sum/n, var = (sum2 - n*mean*mean)/(n-1);
    double sigma = (var > 0) ? std::sqrt(var) : 0;

    for (int iter = 0; iter < maxIter; iter++) {
        if (sigma <= 0) break;
        double lo = mean - A*sigma, hi = mean + A*sigma;
        sum = 0; sum2 = 0; long m = 0;
        for (int i = 0; i < n; i++)
            if (wf[i] >= lo && wf[i] <= hi) { sum += wf[i]; sum2 += (double)wf[i]*wf[i]; m++; }
        if (m < 2) break;
        double newMean = sum/m;
        var = (sum2 - m*newMean*newMean)/(m-1);
        double newSigma = (var > 0) ? std::sqrt(var) : 0;
        bool converged = std::fabs(newMean - mean) < tol;
        mean = newMean; sigma = newSigma;
        if (converged) break;
    }
    if (sigmaOut) *sigmaOut = sigma;
    return mean;
}