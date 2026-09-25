// Waveform.h
#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <vector>
#include <string>

class TCanvas;

class Waveform {
public:
    // Confirmed via direct register readout (0x1028 = 0x0, matching the
    // manual's documented 2 Vpp default) -- see run-plan Sec. 3.1.1/3.2.
    static constexpr double FSR_VOLTS = 2.0;
    static constexpr int    ADC_BITS  = 14;
    static constexpr double ADC_COUNTS = 16384.0; // 2^14
    static constexpr double MV_PER_COUNT = (FSR_VOLTS * 1000.0) / ADC_COUNTS; // ~0.1221 mV/count

    Waveform(const int* wf, int nSamples, bool keepFullWaveform = true);

    double getBaseline() const;

    double getMaxAmpADC()   const; // raw ADC counts
    double getMaxAmpVolts() const; // converted, expressed in mV

    const std::vector<double>& getSubtractedADC()   const; // raw ADC counts
    std::vector<double>        getSubtractedVolts() const; // converted, expressed in mV

    // Finds distinct peaks in the baseline-subtracted waveform, matching
    // the method in arXiv:2511.04724 Sec 3.2: a sample counts as a peak
    // if it's a local max AND its rising slope exceeds derivativeThreshold
    // AND its amplitude exceeds amplitudeThreshold. Peaks closer together
    // than minSeparation samples are merged (only the first counts).
    //
    // amplitudeThreshold = 50 ADC (the real hardware self-trigger
    // threshold -- everything from 50 upward counts).
    //
    // derivativeThreshold: not given a number in the paper (setup-specific).
    // Derived here from the known noise level: per-sample noise sigma ~5 ADC
    // (datasheet/measured), so the difference of two adjacent noise samples
    // has sigma ~5*sqrt(2) ~7 ADC. Default below uses 3x that (~21 ADC/sample)
    // so real rising edges trigger it but noise shouldn't. Adjust if this
    // doesn't match what you see in practice.
    std::vector<double> findPeaks(double amplitudeThreshold = 50.0,
                                   double derivativeThreshold = 21.0,
                                   int minSeparation = 10);

    // Scans the whole record for any stretch of at least minStableSamples
    // consecutive samples above amplitudeThreshold, wherever it occurs
    // (doesn't assume the first rise is the real one -- a small unrelated
    // blip can trigger the record before the real sustained pulse begins).
    // Returns true if such a stretch exists anywhere.
    //
    // amplitudeThreshold = 50 ADC (hardware self-trigger threshold).
    // minStableSamples = 38 (76 ns) -- matches the neutron's dominant
    // 76 ns / 75%-weight decay component from the physics-grounding paper,
    // and ~6x longer than a gamma's single 12 ns decay component, so a
    // gamma pulse alone shouldn't stay elevated long enough to pass this.
    bool isNeutronLike(double amplitudeThreshold = 50.0, int minStableSamples = 38) const;

    // Same as isNeutronLike(), but also reports where the qualifying
    // stretch starts (outStartSample) and its peak amplitude
    // (outPeakAmplitude). If no qualifying stretch is found, both are
    // left at -1 and 0 respectively. Meant for anchoring PSD gates at the
    // actual detected pulse instead of a fixed calculated position.
    bool findNeutronLikeStretch(double amplitudeThreshold, int minStableSamples,
                                 int& outStartSample, double& outPeakAmplitude) const;

    // Checks whether the record, from startSample to the end, shows an
    // overall DECLINING envelope -- distinct from isNeutronLike(), which
    // only checks that the signal stays elevated, not that it trends
    // downward. Real neutron/gamma decays should trend down overall (a
    // sum of decaying exponentials), even with individual bumps riding on
    // top; some false-positive events (e.g. a large gamma interaction
    // with an electronics/PMT recovery artifact, or continuously-busy
    // pile-up) can pass isNeutronLike() while never actually declining --
    // confirmed directly on two real events with the same isNeutronLike()
    // classification but visibly different envelope behavior.
    //
    // Splits [startSample, n) into numChunks equal pieces, takes the
    // MEDIAN of each (not mean, so one large bump can't skew a chunk),
    // then checks how many consecutive chunk-pairs decline. Returns true
    // if at least minDecliningPairs out of (numChunks-1) pairs decline --
    // default requires 3 of 4 (one exception allowed, for a single late
    // bump that doesn't reflect the overall trend).
    bool hasDecliningEnvelope(int startSample, int numChunks = 5, int minDecliningPairs = 3) const;

    // Counts genuinely SEPARATE excursions above threshold -- unlike
    // findPeaks() above (which counts every local max and fails badly on
    // bumpy real events), this only starts a NEW count after the signal
    // has been quiet (at/below threshold) for at least minQuietSamples in
    // a row. A bumpy-but-continuously-elevated event (real neutron or
    // gamma-artifact tail) still counts as ONE excursion; only a genuine
    // sustained return to baseline followed by a new rise counts as a
    // second, separate one.
    //
    // Meant for splitting the not-neutron-like sample (gammas + pile-up)
    // into single-pulse (keep) vs multi-pulse (pile-up, exclude):
    // countSeparatePulses().size() == 1 -> single pulse, keep.
    // size() >= 2 -> pile-up, exclude.
    //
    // amplitudeThreshold = 50 ADC. minQuietSamples = 15 (30 ns) by
    // default -- verified against synthetic single-bumpy-event (1),
    // two-separate-pulses (2), and pure-noise (0) cases.
    std::vector<double> countSeparatePulses(double amplitudeThreshold = 50.0,
                                             int minQuietSamples = 15) const;

    // Same as countSeparatePulses(), but also reports the start sample of
    // each accepted excursion (parallel array to the returned amplitudes).
    // Meant for anchoring PSD gates at the actual detected pulse instead
    // of a fixed calculated position -- most useful when there's exactly
    // one pulse (outStartSamples[0]).
    std::vector<double> countSeparatePulses(double amplitudeThreshold, int minQuietSamples,
                                             std::vector<int>& outStartSamples) const;

    // Baseline window width, as a fraction of the record. Set once before
    // looping over events; applies to every Waveform constructed after.
    //
    // Default is 0.10 (~103 samples for a 1030-sample record). Widened
    // from the earlier 0.05 because a 51-sample window left too much
    // statistical error on the baseline estimate (std err ~ sigma/sqrt(51)
    // ~ 0.7 ADC), and that error accumulates over a long integration gate
    // (~1.3 ADC/sample x 130 samples = ~166 ADC of fake charge, seen
    // directly on Event 47).
    //
    // Upper limit: the real pulse rise sits around sample 128 (~12.4% of
    // the record), so going much past 0.10 risks the window overrunning
    // into the rising edge.
    static void   setPretriggerFraction(double fraction);
    static double getPretriggerFraction();

    // Draws this event's baseline-subtracted waveform with:
    //   - a horizontal dashed line at amplitudeThreshold
    //   - if a sustained stretch of >= minStableSamples above threshold is
    //     found, two vertical dashed lines marking where it starts and
    //     where it reaches the required duration (i.e. the actual window
    //     that satisfies the isNeutronLike() condition)
    //   - if extraLabel is non-empty, draws it as a text label in the
    //     upper-left corner (e.g. "PSD = 0.823")
    // Prints to the currently-open multi-page PDF (c and pdfName are
    // managed by the caller).
    void drawWithThreshold(double amplitudeThreshold, int minStableSamples, long eventIndex,
                            double nsPerSample, TCanvas* c, const std::string& pdfName,
                            const std::string& extraLabel = "") const;

private:
    // Baseline is computed via iterative 3-sigma clipping, but restricted
    // to only the first PRETRIGGER_FRACTION of the record -- not the mean
    // (which had no protection against a spike landing inside the window
    // itself, confirmed directly on an event with a large pulse at t~0-20ns
    // that dragged the whole-window mean off, shifting the entire
    // subtracted trace by ~-50 to -70 ADC), and not the OLD whole-record
    // sigma-cut (which was removed: it couldn't converge for events where
    // a large fraction of the ENTIRE record is disturbed, e.g. Event 82032).
    // Restricting the clip to just the short pretrigger window lets it
    // successfully exclude a small number of spike samples and converge
    // on the true quiet baseline from what's left, without the whole-record
    // failure mode.
    static double pretriggerFraction; // see setPretriggerFraction(), default 0.10
    static constexpr double PRETRIGGER_SIGMA_CUT_A     = 3.0;
    static constexpr int    PRETRIGGER_SIGMA_CUT_ITERS = 5;
    static constexpr double PRETRIGGER_SIGMA_CUT_TOL   = 0.01;

    int n;
    double baseline;
    double maxAmp;
    std::vector<double> subtracted;

    static double computePretriggerBaseline(const int* wf, int n, double fraction,
                                             double A, int maxIter, double tol);
};

#endif // WAVEFORM_H