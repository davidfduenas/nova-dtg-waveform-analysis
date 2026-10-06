// DecayAnalysis.h
// Averages many single-pulse waveforms and fits the average, to extract the
// physical shape of the pulse: its rise time, and the decay time and weight of
// each exponential component.
//
//  - Every waveform is lined up at the point where it first reaches
//    alignFraction of its height, and divided by its area over the first areaNs.
//    (The area is much less noisy than the peak, so scaling by it does not
//    distort the early shape.) The scaled waveforms are averaged sample by
//    sample, and the error of the average at each sample is kept.
//  - The average is fitted with
//        c * step(t) + sum_i A_i * exp(-t / tau_i)
//    c is an optional constant offset that stands for a component too slow to
//    measure in the record.
//  - The fit starts at the peak of the average, not on the rising edge: the
//    decay times then do not depend on how the rising edge is modelled (a fit
//    that includes the edge was found to bias them). If fitRise is true, t0 and
//    sigma are fitted too: the step and the exponentials start at t0 and are
//    smoothed by a Gaussian of width sigma, which describes the rounding of the
//    pulse near its peak. If fitRise is false the exponentials are not smoothed,
//    which is right only if the fit starts well after the peak.
//
// All numbers come in through Settings; none are fixed in the class.
#ifndef DECAYANALYSIS_H
#define DECAYANALYSIS_H

#include <vector>
#include <string>

class DecayAnalysis {
public:
    struct Settings {
        double nsPerSample;
        double alignFraction;    // waveforms are lined up where they first reach this fraction of their height
        int    searchBefore;     // samples before the given onset where that crossing is looked for
        int    searchAfter;      // ... and after it
        int    peakSamples;      // the height is the highest point within this many samples after the onset
        double areaNs;           // each waveform is divided by its area over this long after the alignment point
        double preNs;            // the average starts this long before the alignment point
        double postNs;           // ... and ends this long after it
        double fitEndNs;         // the fits use the points up to this time
        double systematicFloor;  // the error of each point is at least this fraction of the average's maximum
        std::vector<double> sigmaStarts;   // starting values (ns) of the rise width sigma
        std::vector<double> fastStarts;    // starting decay times (ns) of the first, second and third
        std::vector<double> mediumStarts;  // component. Every combination is tried and the best fit kept.
        std::vector<double> slowStarts;    // (slowStarts is used only when three components are fitted)
        int    maxEvaluations;   // limit for each minimisation
        double tolerance;        // a minimisation stops when chi2 changes by less than this fraction
        long   minEvents;        // an average with fewer events than this is not fitted
    };

    struct Average {
        std::vector<double> t;      // ns after the alignment point
        std::vector<double> mean;   // normalised amplitude (1 / ns)
        std::vector<double> error;  // error of the mean
        long events;
    };

    struct Component {
        double tau, tauErr;                 // decay time (ns)
        double amplitude, amplitudeErr;     // A_i (1 / ns)
        double fraction, fractionErr;       // A_i / (c + sum of all A): the weights, as in the literature
        double share, shareErr;             // share of the light in the first areaNs
    };

    struct FitResult {
        bool   ok;                // false if the fit could not be done
        bool   fitRise;           // true if t0 and sigma were fitted (otherwise they are fixed and not reported)
        double fitStartNs, fitEndNs;
        bool   errorsOk;          // false if the errors could not be computed
        bool   withOffset;
        int    nComponents;
        long   events;
        double t0, t0Err;
        double sigma, sigmaErr;
        double rise1090, rise1090Err;
        double offset, offsetErr;
        double offsetFraction, offsetFractionErr;
        double offsetShare, offsetShareErr;
        double shareSum;          // sum of all shares: should be close to 1
        std::vector<Component> components;   // ordered by decay time, shortest first
        double chi2;
        int    ndf;
        double errorScale;        // the errors are multiplied by this (sqrt(chi2/ndf), at least 1)
        bool   allPositive;       // false if a component came out with a negative amplitude
        // The points that were fitted, for plotting.
        std::vector<double> t, y, err, model, pull;
    };

    explicit DecayAnalysis(const Settings& settings);

    // Add one waveform (baseline-subtracted, pulse upward). onsetSample is where
    // its pulse starts, from the pulse finder. Returns false if the waveform
    // could not be lined up or does not fit in the record.
    bool add(const std::vector<double>& x, int onsetSample);

    long numEvents() const { return events; }

    // Time (ns after the alignment point) at which the average is highest.
    double peakTimeNs() const;

    Average average() const;

    // Fit the average of the waveforms added so far, using the points from
    // startNs to the end of the fit range. If fitRise is false the rising edge is
    // not modelled, so startNs should be after it.
    FitResult fit(int nComponents, bool withOffset, double startNs, bool fitRise) const;

    // Fit any average (the same machinery, so it can be tested on known shapes).
    FitResult fitData(const std::vector<double>& t, const std::vector<double>& y,
                      const std::vector<double>& err, int nComponents, bool withOffset,
                      double startNs, bool fitRise, long eventCount) const;

    // Print a fit result on the console.
    static void print(const FitResult& r, const std::string& label);

    // Write one plot to a PDF: the average pulse (normalised to its maximum) in
    // grey, with the fit as a black dashed line, against time.
    void writePdf(const std::string& pdfPath, const FitResult& fit) const;

private:
    Settings cfg;
    long events;
    std::vector<double> sumY;
    std::vector<double> sumYY;

    int preSamples() const;
    int postSamples() const;
};

#endif // DECAYANALYSIS_H
