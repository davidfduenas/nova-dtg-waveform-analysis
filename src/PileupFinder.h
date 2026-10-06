// PileupFinder.h
// Finds pulses in a baseline-subtracted waveform in two steps.
//  1. A climb is a place where the average of the next windowSamples samples is
//     at least 'jump' ADC above the average of the previous windowSamples
//     samples. Each climb is counted once.
//  2. The height of a climb is its peak (within peakSamples samples) minus the
//     average just before it. The first climb in the record is always a pulse
//     (it is what triggered the digitizer). A later climb is a pulse if its
//     height is at least minFraction of the biggest peak in the event AND at
//     least noiseSigmas times the noise just before it (the scatter of the
//     signal around a straight line over the previous noiseSamples samples).
//     Climbs that start within blankSamples of a pulse are ignored, because the
//     top of a pulse is jagged.
// All numbers come in through Settings; none are fixed in the class.
#ifndef PILEUPFINDER_H
#define PILEUPFINDER_H

#include <vector>

class PileupFinder {
public:
    struct Settings {
        double jump;          // smallest climb to look at (ADC)
        double minFraction;   // a later pulse is at least this fraction of the event's peak
        double noiseSigmas;   // ... and at least this many times the noise before it
        int    noiseSamples;  // the noise is measured over this many samples before the climb
        int    windowSamples; // samples averaged on each side of a point
        int    peakSamples;   // the peak is the highest point within this many samples
                              // after the start of the climb
        int    blankSamples;  // climbs this soon after a pulse are ignored
        int    maxPulses;     // stop after this many pulses in one record
    };

    struct Pulse {
        int    onsetSample; // first sample of the climb
        double height;      // peak above the level before it (ADC)
        double noise;       // scatter of the signal just before the climb (ADC)
    };

    explicit PileupFinder(const Settings& settings) : cfg(settings) {}

    std::vector<Pulse> find(const std::vector<double>& x) const;

private:
    Settings cfg;

    // Scatter of x around a straight line over the noiseSamples samples before i.
    double noiseBefore(const std::vector<double>& x, int i) const;
};

#endif // PILEUPFINDER_H
