// Waveform.h
#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <vector>

class Waveform {
public:
    static constexpr double SIGMA_CUT_A        = 3.0;
    static constexpr int    SIGMA_CUT_MAX_ITER = 5;
    static constexpr double SIGMA_CUT_TOL      = 0.01;

    // Confirmed via direct register readout (0x1028 = 0x0, matching the
    // manual's documented 2 Vpp default) -- see run-plan Sec. 3.1.1/3.2.
    static constexpr double FSR_VOLTS = 2.0;
    static constexpr int    ADC_BITS  = 14;
    static constexpr double ADC_COUNTS = 16384.0; // 2^14
    static constexpr double MV_PER_COUNT = (FSR_VOLTS * 1000.0) / ADC_COUNTS; // ~0.1221 mV/count

    Waveform(const int* wf, int nSamples, bool keepFullWaveform = true);

    double getBaseline() const;
    double getSigma()    const;

    double getMaxAmpADC()   const; // raw ADC counts
    double getMaxAmpVolts() const; // converted, expressed in mV

    const std::vector<double>& getSubtractedADC()   const; // raw ADC counts
    std::vector<double>        getSubtractedVolts() const; // converted, expressed in mV

    std::vector<double> findPeaks(double threshold, int minSeparation);

private:
    int n;
    double baseline;
    double sigma;
    double maxAmp;
    std::vector<double> subtracted;

    static double computeBaselineSigmaCut(const int* wf, int n,
                                           double A, int maxIter,
                                           double tol, double* sigmaOut = nullptr);
};

#endif // WAVEFORM_H