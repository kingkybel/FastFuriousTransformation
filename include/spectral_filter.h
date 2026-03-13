/*
 * File Name:   spectral_filter.h
 * Description: Frequency-domain filters that reuse the FFT helper to mask
 *              bands before reconstructing the time-domain signal.
 *
 * Copyright (C) 2026 Dieter J Kybelksties
 */

#ifndef NS_UTIL_SPECTRAL_FILTER_H_INCLUDED
#define NS_UTIL_SPECTRAL_FILTER_H_INCLUDED

#include "fft.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace util
{
class SpectralFilter
{
  public:
    using FFTType       = FFT;
    using FloatType     = FFTType::FloatType;
    using FloatVector   = FFTType::FloatVector;
    using ComplexVector = FFTType::ComplexVector;
    using IntType       = FFTType::IntType;

    enum class Mode
    {
        LowPass,
        HighPass,
        BandPass
    };

    /**
     * @param logPoints Log2(number of points) used by the internal FFT plan.
     * @param sampleRate Sample rate used for the frequency axis.
     */
    SpectralFilter(IntType logPoints = 10, IntType sampleRate = 1'024)
        : fft_(logPoints, sampleRate)
        , mode_(Mode::LowPass)
        , sampleRate_(static_cast<FloatType>(sampleRate))
    {
        assert(sampleRate > 0);
        setNyquist();
        configureLowPass(sampleRate_ / 4.0L);
    }

    void configureLowPass(FloatType cutoffHz)
    {
        mode_   = Mode::LowPass;
        lowHz_  = 0.0L;
        highHz_ = std::clamp(cutoffHz, FloatType(0.0L), nyquist_);
    }

    void configureHighPass(FloatType cutoffHz)
    {
        mode_   = Mode::HighPass;
        lowHz_  = std::clamp(cutoffHz, FloatType(0.0L), nyquist_);
        highHz_ = nyquist_;
    }

    void configureBandPass(FloatType lowHz, FloatType highHz)
    {
        mode_   = Mode::BandPass;
        lowHz_  = std::clamp(lowHz, FloatType(0.0L), nyquist_);
        highHz_ = std::clamp(highHz, lowHz_, nyquist_);
    }

    /**
     * @brief Apply the configured mask in the frequency domain.
     */
    FloatVector apply(FloatVector const &samples)
    {
        assert(static_cast<IntType>(samples.size()) == fft_.numberOfPoints());
        fft_.loadFloatVector(samples);
        ComplexVector spectrum = fft_.transform();

        for (size_t bin = 0; bin < spectrum.size(); ++bin)
        {
            FloatType magnitudeFreq = binFrequencyFromMagnitude(bin);

            if (!shouldKeepFrequency(magnitudeFreq))
            {
                spectrum[bin] = ComplexVector::value_type(0.0L, 0.0L);
            }
        }

        ComplexVector filtered = fft_.inverseTransformFromSpectrum(spectrum);
        FloatVector   result(filtered.size());

        for (size_t i = 0; i < filtered.size(); ++i)
        {
            result[i] = filtered[i].real();
        }

        return result;
    }

    Mode mode() const
    {
        return mode_;
    }

    FloatType lowCutoff() const
    {
        return lowHz_;
    }

    FloatType highCutoff() const
    {
        return highHz_;
    }

  private:
    void setNyquist()
    {
        nyquist_ = sampleRate_ / 2.0L;
        binStep_ = sampleRate_ / static_cast<FloatType>(fft_.numberOfPoints());
    }

    FloatType binFrequencyFromMagnitude(size_t bin) const
    {
        size_t half     = fft_.numberOfPoints() / 2;
        size_t mirrored = (bin <= half) ? bin : fft_.numberOfPoints() - bin;
        return static_cast<FloatType>(mirrored) * binStep_;
    }

    bool shouldKeepFrequency(FloatType frequency) const
    {
        switch (mode_)
        {
            case Mode::LowPass:
                return frequency <= highHz_;
            case Mode::HighPass:
                return frequency >= lowHz_;
            case Mode::BandPass:
                return frequency >= lowHz_ && frequency <= highHz_;
        }
        return false;
    }

    FFTType   fft_;
    Mode      mode_;
    FloatType sampleRate_;
    FloatType nyquist_;
    FloatType binStep_;
    FloatType lowHz_;
    FloatType highHz_;
};
} // namespace util

#endif // NS_UTIL_SPECTRAL_FILTER_H_INCLUDED
