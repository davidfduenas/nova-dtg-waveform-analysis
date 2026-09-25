// PSD.h
//
// Owns the whole PSD pipeline: pulse-start derivation, gate integration,
// PSD ratio calculation, histogram accumulation, and saving plots. The
// macro that uses this just feeds events in one at a time and calls
// savePlots() once at the end -- nothing plot-related happens outside
// this class.
//
// All configuration values (post-trigger percent, ConstantLatency, gate
// widths, amplitude cutoff, etc.) are passed in through the constructor.
// If a file was acquired with different settings, the CALLER needs to
// pass different values in -- this class does not assume or hardcode
// anything about a specific run.

#ifndef PSD_H
#define PSD_H

#include <vector>
#include <string>

class TH1D;
class TH2D;
class TCanvas;

class PSD {
public:
    struct Result {
        int    pulseStartSample;
        int    shortGateEndSample;
        int    longGateEndSample;
        double qShort;
        double qLong;
        double psdRatio; // (qLong - qShort) / qLong
    };

    // postTriggerPercent    : the --post value actually used for this file (e.g. 80.0)
    // nCoefficient          : 8 for 730 family, 4 for 725 family (UM5118 Sec 1.24)
    // constantLatencySamples: measured value, see measureConstantLatency.cpp
    //                         (depends on trigger source/firmware, NOT on
    //                         record length or post-trigger percent -- but
    //                         must be re-measured if trigger source or
    //                         firmware changes)
    // shortGateNs, longGateNs: gate widths from pulse start
    // nsPerSample           : 2.0 for this 500 MS/s digitizer
    // amplitudeCutoff       : only events with max amplitude above this
    //                         get included in the histograms
    // startShiftNs          : the pre-gate -- shifts BOTH gates earlier by
    //                         this many ns, applied uniformly whether the
    //                         anchor comes from the fixed formula
    //                         (analyze()) or an explicit detected location
    //                         (analyzeAt()). Matches the reference paper's
    //                         explicit pre-gate before charge integration
    //                         starts. Needed because a fast-rising pulse's
    //                         first threshold crossing can land almost on
    //                         the peak itself, leaving no real baseline
    //                         inside the gate otherwise.
    PSD(double postTriggerPercent,
        double nCoefficient,
        double constantLatencySamples,
        double shortGateNs,
        double longGateNs,
        double nsPerSample,
        double amplitudeCutoff,
        double startShiftNs = 0.0);

    ~PSD();

    // Real pulse-start sample for a record of this length, using this
    // instance's configured post-trigger percent/N/ConstantLatency.
    int getPulseStartSample(int recordLength) const;

    // shape = baseline-subtracted samples for one event (Waveform::getSubtractedADC()).
    Result analyze(const std::vector<double>& shape, int recordLength) const;

    // Same as analyze(), but uses explicitStartSample directly instead of
    // computing one from the acquisition-timing formula. Use this when the
    // caller already knows where the real pulse is (e.g. from
    // Waveform::findNeutronLikeStretch() or countSeparatePulses()'s
    // start-sample output) -- keeps amplitude and PSD anchored to the same
    // detected event instead of two independently-derived locations.
    Result analyzeAt(const std::vector<double>& shape, int recordLength,
                      int explicitStartSample) const;

    // Feed one event in. Computes PSD internally and fills the histograms,
    // but only if maxAmp is above the configured amplitude cutoff.
    void addEvent(const std::vector<double>& shape, int recordLength, double maxAmp);

    // Same as addEvent(), but anchored at explicitStartSample (see analyzeAt()).
    void addEventAt(const std::vector<double>& shape, int recordLength,
                     double maxAmp, int explicitStartSample);

    // Draw and save the accumulated histograms as PDFs:
    //   outputPath + "psdRatio_" + label + ".pdf"
    //   outputPath + "psdVsAmplitude_" + label + ".pdf"
    void savePlots(const std::string& outputPath, const std::string& label) const;

    // Draw one event's baseline-subtracted waveform with vertical lines
    // marking pulse start, short gate end, and long gate end, then print
    // it to the currently-open multi-page PDF (c and pdfName are managed
    // by the caller -- this just draws one page and prints it).
    void drawEventWaveform(const std::vector<double>& shape, int recordLength,
                           long eventIndex, TCanvas* c, const std::string& pdfName) const;

    // Same as drawEventWaveform(), but anchored at explicitStartSample (see analyzeAt()).
    void drawEventWaveformAt(const std::vector<double>& shape, int recordLength,
                              int explicitStartSample, long eventIndex,
                              TCanvas* c, const std::string& pdfName) const;

    long getCandidateCount() const;
    long getTotalCount() const;

    // Draws both samples' PSD-vs-amplitude histograms on one canvas:
    // denseSample as a log-z COLZ background, sparseSample overlaid as a
    // distinct-colored scatter on top -- since the two samples can differ
    // by hundreds of times in event count, a single shared color scale
    // would make the smaller one invisible.
    // Saves to outputPath + "psdVsAmplitude_combined_" + label + ".pdf"
    static void saveCombinedPlot(const PSD& denseSample, const PSD& sparseSample,
                                  const std::string& outputPath, const std::string& label);

private:
    double postTriggerPercent;
    double nCoefficient;
    double constantLatencySamples;
    double shortGateNs;
    double longGateNs;
    double nsPerSample;
    double amplitudeCutoff;
    double startShiftNs;

    long totalCount;
    long candidateCount;

    TH1D* hPSD;
    TH2D* hPSDvsAmp;

    // Shared implementation: both analyze() and analyzeAt() funnel through
    // this once a start sample is decided (either computed or explicit).
    Result analyzeFromStart(const std::vector<double>& shape, int recordLength,
                             int startSample) const;

    void drawFromStart(const std::vector<double>& shape, int recordLength, int startSample,
                        long eventIndex, TCanvas* c, const std::string& pdfName) const;
};

#endif // PSD_H