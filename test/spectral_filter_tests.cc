/*
 * File:		SpectralFilter_tests.cc
 * Description:	Unit tests for the FFT-backed spectral filter helpers.
 */

#include "spectral_filter.h"

#include <cassert>
#include <cmath>
#include <gtest/gtest.h>

using namespace util;

namespace
{
FFT::FloatVector makeTone(size_t length, double sampleRate, double freq)
{
    FFT::FloatVector result(length);
    for (size_t i = 0; i < length; ++i)
    {
        result[i] = std::cos(2.0L * M_PI * freq * static_cast<double>(i) / sampleRate);
    }
    return result;
}

FFT::FloatVector addVectors(FFT::FloatVector const &a, FFT::FloatVector const &b)
{
    assert(a.size() == b.size());
    FFT::FloatVector out(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        out[i] = a[i] + b[i];
    }
    return out;
}
} // namespace

TEST(SpectralFilterTest, low_pass_removes_high_frequency)
{
    SpectralFilter filter(10, 1'024);
    filter.configureLowPass(30.0L);

    constexpr size_t size     = 1'024;
    auto             lowTone  = makeTone(size, 1'024.0, 5.0);
    auto             highTone = makeTone(size, 1'024.0, 200.0);
    auto             mixture  = addVectors(lowTone, highTone);

    auto   filtered = filter.apply(mixture);
    double maxDiff  = 0.0;
    for (size_t i = 0; i < size; ++i)
    {
        double diff = static_cast<double>(std::abs(filtered[i] - lowTone[i]));
        maxDiff     = std::max(maxDiff, diff);
    }

    EXPECT_LT(maxDiff, 2e-6);
}

TEST(SpectralFilterTest, high_pass_removes_low_frequency)
{
    SpectralFilter filter(10, 1'024);
    filter.configureHighPass(60.0L);

    constexpr size_t size     = 1'024;
    auto             lowTone  = makeTone(size, 1'024.0, 5.0);
    auto             highTone = makeTone(size, 1'024.0, 200.0);
    auto             mixture  = addVectors(lowTone, highTone);

    auto   filtered = filter.apply(mixture);
    double maxDiff  = 0.0;
    for (size_t i = 0; i < size; ++i)
    {
        double diff = static_cast<double>(std::abs(filtered[i] - highTone[i]));
        maxDiff     = std::max(maxDiff, diff);
    }

    EXPECT_LT(maxDiff, 2e-6);
}

TEST(SpectralFilterTest, band_pass_isolation)
{
    SpectralFilter filter(10, 1'024);
    filter.configureBandPass(35.0L, 80.0L);

    constexpr size_t size     = 1'024;
    auto             lowTone  = makeTone(size, 1'024.0, 5.0);
    auto             bandTone = makeTone(size, 1'024.0, 60.0);
    auto             highTone = makeTone(size, 1'024.0, 200.0);
    auto             mixture  = addVectors(addVectors(lowTone, bandTone), highTone);

    auto   filtered = filter.apply(mixture);
    double maxDiff  = 0.0;
    for (size_t i = 0; i < size; ++i)
    {
        double diff = static_cast<double>(std::abs(filtered[i] - bandTone[i]));
        maxDiff     = std::max(maxDiff, diff);
    }

    EXPECT_LT(maxDiff, 3e-6);
}
