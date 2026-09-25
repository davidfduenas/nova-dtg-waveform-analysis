// Waveform.cpp
#include "Waveform.h"
#include <cmath>
#include <algorithm>
#include <vector>

#include <TH1D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>

// Default baseline window: 10% of the record.
double Waveform::pretriggerFraction = 0.10;

void   Waveform::setPretriggerFraction(double fraction) { pretriggerFraction = fraction; }
double Waveform::getPretriggerFraction() { return pretriggerFraction; }

Waveform::Waveform(const int* wf, int nSamples, bool keepFullWaveform)
    : n(nSamples)
{
    baseline = computePretriggerBaseline(wf, n, pretriggerFraction,
                                          PRETRIGGER_SIGMA_CUT_A,
                                          PRETRIGGER_SIGMA_CUT_ITERS,
                                          PRETRIGGER_SIGMA_CUT_TOL);

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

double Waveform::getMaxAmpADC()   const { return maxAmp; }
double Waveform::getMaxAmpVolts() const { return maxAmp * MV_PER_COUNT; }

const std::vector<double>& Waveform::getSubtractedADC() const { return subtracted; }

std::vector<double> Waveform::getSubtractedVolts() const {
    std::vector<double> mv(subtracted.size());
    for (size_t i = 0; i < subtracted.size(); i++) mv[i] = subtracted[i] * MV_PER_COUNT;
    return mv;
}

std::vector<double> Waveform::findPeaks(double amplitudeThreshold, double derivativeThreshold, int minSeparation)
{
    std::vector<double> peaks;
    int lastPeakSample = -minSeparation; // so the very first peak is never rejected

    for (int i = 1; i < n - 1; i++) {
        bool isLocalMax = (subtracted[i] >= subtracted[i-1]) && (subtracted[i] >= subtracted[i+1]);
        if (!isLocalMax) continue;

        if (subtracted[i] <= amplitudeThreshold) continue;

        double risingSlope = subtracted[i] - subtracted[i-1];
        if (risingSlope <= derivativeThreshold) continue;

        if ((i - lastPeakSample) < minSeparation) continue; // too close to the last accepted peak

        peaks.push_back(subtracted[i]);
        lastPeakSample = i;
    }

    return peaks;
}

bool Waveform::isNeutronLike(double amplitudeThreshold, int minStableSamples) const
{
    int run = 0;
    for (int i = 0; i < n; i++) {
        if (subtracted[i] > amplitudeThreshold) {
            run++;
            if (run >= minStableSamples) return true;
        } else {
            run = 0;
        }
    }
    return false;
}

bool Waveform::findNeutronLikeStretch(double amplitudeThreshold, int minStableSamples,
                                       int& outStartSample, double& outPeakAmplitude) const
{
    outStartSample = -1;
    outPeakAmplitude = 0;

    int run = 0, runStart = -1;
    double runMax = 0;
    for (int i = 0; i < n; i++) {
        if (subtracted[i] > amplitudeThreshold) {
            if (run == 0) { runStart = i; runMax = subtracted[i]; }
            else if (subtracted[i] > runMax) runMax = subtracted[i];
            run++;
            if (run >= minStableSamples) {
                outStartSample = runStart;
                outPeakAmplitude = runMax;
                return true;
            }
        } else {
            run = 0;
        }
    }
    return false;
}

bool Waveform::hasDecliningEnvelope(int startSample, int numChunks, int minDecliningPairs) const
{
    if (startSample < 0 || startSample >= n) return false;
    int rangeLen = n - startSample;
    if (rangeLen < numChunks) return false; // not enough samples to chunk meaningfully

    std::vector<double> chunkMedians(numChunks);
    int chunkSize = rangeLen / numChunks;

    for (int c = 0; c < numChunks; c++) {
        int chunkStart = startSample + c * chunkSize;
        int chunkEnd = (c == numChunks - 1) ? n : chunkStart + chunkSize; // last chunk absorbs remainder

        std::vector<double> vals(subtracted.begin() + chunkStart, subtracted.begin() + chunkEnd);
        std::sort(vals.begin(), vals.end());
        size_t m = vals.size();
        double median = (m % 2 == 0) ? 0.5 * (vals[m/2 - 1] + vals[m/2]) : vals[m/2];
        chunkMedians[c] = median;
    }

    int decliningPairs = 0;
    for (int c = 0; c < numChunks - 1; c++) {
        if (chunkMedians[c] > chunkMedians[c + 1]) decliningPairs++;
    }

    return decliningPairs >= minDecliningPairs;
}

std::vector<double> Waveform::countSeparatePulses(double amplitudeThreshold, int minQuietSamples) const
{
    std::vector<int> unused;
    return countSeparatePulses(amplitudeThreshold, minQuietSamples, unused);
}

std::vector<double> Waveform::countSeparatePulses(double amplitudeThreshold, int minQuietSamples,
                                                    std::vector<int>& outStartSamples) const
{
    std::vector<double> peaks;
    outStartSamples.clear();

    bool elevated = false;
    int quietRun = 0;
    double currentPeakMax = 0;
    int currentPeakStart = -1;

    for (int i = 0; i < n; i++) {
        double v = subtracted[i];

        if (v > amplitudeThreshold) {
            if (!elevated) {
                elevated = true;
                currentPeakMax = v;
                currentPeakStart = i;
            } else {
                if (v > currentPeakMax) currentPeakMax = v;
            }
            quietRun = 0;
        } else {
            if (elevated) {
                quietRun++;
                if (quietRun >= minQuietSamples) {
                    peaks.push_back(currentPeakMax);
                    outStartSamples.push_back(currentPeakStart);
                    elevated = false;
                    quietRun = 0;
                }
            }
        }
    }

    if (elevated) {
        peaks.push_back(currentPeakMax);
        outStartSamples.push_back(currentPeakStart);
    }

    return peaks;
}

double Waveform::computePretriggerBaseline(const int* wf, int n, double fraction,
                                            double A, int maxIter, double tol)
{
    int nPre = (int)std::round(fraction * n);
    if (nPre < 2) nPre = std::min(2, n);

    // Start from the plain mean/sigma of the pretrigger window.
    double sum = 0, sum2 = 0;
    for (int i = 0; i < nPre; i++) { sum += wf[i]; sum2 += (double)wf[i] * wf[i]; }
    double mean = sum / nPre;
    double var  = (sum2 - nPre * mean * mean) / (nPre - 1);
    double sigma = (var > 0) ? std::sqrt(var) : 0;

    // Iteratively clip samples more than A*sigma from the mean, restricted
    // to just this short window -- unlike the removed whole-record
    // version, a small number of spike samples here can actually be
    // excluded while leaving enough genuinely quiet samples to converge.
    for (int iter = 0; iter < maxIter; iter++) {
        if (sigma <= 0) break;
        double lo = mean - A * sigma, hi = mean + A * sigma;
        double s = 0, s2 = 0; int m = 0;
        for (int i = 0; i < nPre; i++) {
            if (wf[i] >= lo && wf[i] <= hi) { s += wf[i]; s2 += (double)wf[i] * wf[i]; m++; }
        }
        if (m < 2) break; // nothing left to converge on -- keep prior mean
        double newMean = s / m;
        var = (s2 - m * newMean * newMean) / (m - 1);
        double newSigma = (var > 0) ? std::sqrt(var) : 0;
        bool converged = std::fabs(newMean - mean) < tol;
        mean = newMean;
        sigma = newSigma;
        if (converged) break;
    }

    return mean;
}

void Waveform::drawWithThreshold(double amplitudeThreshold, int minStableSamples, long eventIndex,
                                  double nsPerSample, TCanvas* c, const std::string& pdfName,
                                  const std::string& extraLabel) const
{
    std::string hname = "hThresh" + std::to_string(eventIndex);
    TH1D* h = new TH1D(hname.c_str(),
                        ("Event " + std::to_string(eventIndex) + ";Time (ns);Baseline - Signal (ADC)").c_str(),
                        n, 0, n * nsPerSample);
    for (int s = 0; s < n; s++) h->SetBinContent(s + 1, subtracted[s]);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h->Draw("HIST");

    double yMin = h->GetMinimum();
    double yMax = h->GetMaximum();

    TLine* threshLine = new TLine(0, amplitudeThreshold, n * nsPerSample, amplitudeThreshold);
    threshLine->SetLineColor(kMagenta);
    threshLine->SetLineStyle(2);
    threshLine->SetLineWidth(1);
    threshLine->Draw();

    // Find the first stretch of >= minStableSamples consecutive samples
    // above threshold, same logic as isNeutronLike() but also tracking
    // where it starts, so we can mark the actual qualifying window.
    int runStart = -1, foundStart = -1;
    int run = 0;
    for (int i = 0; i < n; i++) {
        if (subtracted[i] > amplitudeThreshold) {
            if (run == 0) runStart = i;
            run++;
            if (run >= minStableSamples) { foundStart = runStart; break; }
        } else {
            run = 0;
        }
    }

    std::vector<TLine*> durationLines;
    if (foundStart >= 0) {
        double xStart = foundStart * nsPerSample;
        double xEnd   = (foundStart + minStableSamples) * nsPerSample;

        TLine* lStart = new TLine(xStart, yMin, xStart, yMax);
        lStart->SetLineColor(kOrange + 7);
        lStart->SetLineStyle(3);
        lStart->SetLineWidth(1);
        lStart->Draw();
        durationLines.push_back(lStart);

        TLine* lEnd = new TLine(xEnd, yMin, xEnd, yMax);
        lEnd->SetLineColor(kOrange + 7);
        lEnd->SetLineStyle(3);
        lEnd->SetLineWidth(1);
        lEnd->Draw();
        durationLines.push_back(lEnd);
    }

    if (!extraLabel.empty()) {
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.15, 0.85, extraLabel.c_str());
    }

    c->Print(pdfName.c_str());

    delete threshLine;
    for (auto* l : durationLines) delete l;
    delete h;
}