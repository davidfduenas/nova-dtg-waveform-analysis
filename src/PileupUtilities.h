// PileupUtilities.h
// Helpers for the pile-up study:
//  - collects events and writes them to PDFs with each detected pulse marked
//    (1 pulse, 2 pulses and 3 or more in separate files, split into parts of
//    pagesPerPdf pages each);
//  - records the height of every jump found after the first pulse, to show
//    where the jump limit should be.
// All numbers come in through Settings; none are fixed in the class.
#ifndef PILEUPUTILITIES_H
#define PILEUPUTILITIES_H

#include "PileupFinder.h"

#include <vector>
#include <string>

class TCanvas;
class TH1D;
class TH2D;

class PileupUtilities {
public:
    struct Settings {
        double nsPerSample;
        int    pagesPerPdf;      // pages in one PDF file before a new part is started
        long   maxStored;        // most events kept in memory for plotting
        int    canvasWidth;
        int    canvasHeight;
        double headroom;         // y axis of an event page runs to headroom x its maximum
        int    maxPulsesListed;  // pulses written as text on an event page

        int    heightBins;       // jump height histogram: 0 ... heightMax (ADC)
        double heightMax;
        int    ratioBins;        // jump height / first pulse height: 0 ... ratioMax
        double ratioMax;

        // Statistics plots (pileupStatistics.pdf)
        int    fracBins;         // climb height / event peak: 0 ... fracMax
        double fracMax;
        int    sigmaBins;        // climb height / noise before it: 0 ... sigmaMax
        double sigmaMax;
        double cutFraction;      // the finder's cuts, drawn as lines on the fraction-noise plot
        double cutSigmas;
        int    firstBins;        // height of the first climb: 0 ... firstMax (ADC)
        double firstMax;
        int    laterBins;        // height of the biggest later climb: 0 ... laterMax (ADC)
        double laterMax;
        int    gapBins;          // time from the first climb to a later one: 0 ... gapMax (ns)
        double gapMax;
        double gapFlatEnd;       // latest gap that is possible: record length minus where the
                                 // first pulse sits (ns); ends the flat reference line

        std::vector<double> tableEdges; // height edges (ADC) for the console table
    };

    PileupUtilities(const std::string& outDir, const Settings& settings);
    ~PileupUtilities();
    PileupUtilities(const PileupUtilities&) = delete;
    PileupUtilities& operator=(const PileupUtilities&) = delete;

    // Remember one event and its pulses (waveform kept as float to save memory).
    void add(long eventIndex, const std::vector<double>& shape,
             const std::vector<PileupFinder::Pulse>& pulses, double psd);

    long numStored()  const { return (long)events.size(); }
    long numDropped() const { return nDropped; }

    // Write every stored event: first those with 1 pulse, then 2, then 3 or more.
    void write();

    // Statistics for one selected event (call it for every selected event, also
    // those with one pulse).
    //  wideClimbs: every climb, from a wide-open finder run on the event. Those
    //     found at least guardSamples after the first one go into the height study.
    //  pulses: the pulses the finder counted; they give the gap and fraction plots.
    //  eventPeak: the biggest peak in the event (ADC).
    void addEvent(const std::vector<PileupFinder::Pulse>& wideClimbs,
                  const std::vector<PileupFinder::Pulse>& pulses,
                  double eventPeak, int guardSamples);

    // Five pages: first climb against the biggest later climb; number of events
    // with a later climb of at least a given height; time from the first pulse
    // to later pulses, against a flat (random) reference; fraction of events with
    // a later pulse against the height of the first pulse; and, for every climb
    // after the first, its height as a fraction of the event's peak against its
    // height divided by the noise before it (with the finder's cuts marked), which
    // shows where real pulses and tail bumps separate.
    void writeStatisticsPlots(const std::string& pdfPath);

    // Two pages: jump heights, and heights relative to the first pulse.
    void writeHeightPlots(const std::string& pdfPath);

    // Console table of jump counts by height.
    void printHeightSummary() const;

private:
    struct Stored {
        long   index;
        double psd;
        std::vector<float> shape;
        std::vector<PileupFinder::Pulse> pulses;
    };

    std::string dir;
    Settings cfg;
    long nDropped;
    std::vector<Stored> events;

    TH1D* hHeight;
    TH1D* hRatio;
    TH2D* hFirstLater;
    TH1D* hLaterMax;
    TH1D* hGap;
    TH1D* hFirstAll;
    TH1D* hFirstWithLater;
    TH2D* hFracNoise;
    double guardNs;
    std::vector<long> countByHeight;

    int heightBin(double height) const;
    void writeCategory(TCanvas* c, const std::string& tag, int minP, int maxP);
    void drawEvent(TCanvas* c, const Stored& s, const std::string& pdfName);
};

#endif // PILEUPUTILITIES_H
