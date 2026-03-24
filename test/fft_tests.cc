/*
 * File:		FFTTest.cc
 * Description:         Unit tests for Fast Fourier Transform.
 *
 * Copyright (C) 2023 Dieter J Kybelksties <github@kybelksties.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * @date: 2023-08-28
 * @author: Dieter J Kybelksties
 */

#include "fft.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include <ranges>
#include <string>

using namespace std;
using namespace util;

class FFTTest : public ::testing::Test
{
  protected:
    using FloatMatrix = std::vector<std::vector<FFT::FloatType>>;

    void SetUp() override
    {
        // just in case
    }

    void TearDown() override
    {
        // just in case
    }

    vector<FFT::FloatType> randomSample()
    {
        std::default_random_engine                     generator;
        std::uniform_real_distribution<FFT::FloatType> distribution(0.0L, 255.0);
        sampleVec.clear();
        for (int i = 0; i < 1'024; i++)
        {
            FFT::FloatType val = distribution(generator);
            sampleVec.emplace_back(val);
        }

        return sampleVec;
    }

    vector<FFT::FloatType> initConstant(FFT::FloatType value)
    {
        sampleVec.clear();
        for (int i = 0; i < 1'024; i++)
        {
            sampleVec.emplace_back(value);
        }
        return sampleVec;
    }

    vector<FFT::FloatType> initAlternate(FFT::FloatType value1, FFT::FloatType value2)
    {
        sampleVec.clear();
        for (int i = 0; i < 512; i++)
        {
            sampleVec.emplace_back(value1);
            sampleVec.emplace_back(value2);
        }
        return sampleVec;
    }

    vector<FFT::FloatType> circularShift(vector<FFT::FloatType> const& vec, int shift) const
    {
        vector<FFT::FloatType> shifted_vec = vec;
        std::ranges::rotate(shifted_vec, shifted_vec.begin() + shift);
        return shifted_vec;
    }

    vector<FFT::FloatType> initImpulse(int position)
    {
        sampleVec.clear();
        sampleVec.resize(1'024, 0.0);
        sampleVec[position] = 1.0;
        return sampleVec;
    }

    vector<FFT::FloatType> initCosine(FFT::FloatType k)
    {
        sampleVec.clear();
        int const N = 1'024;
        for (int i = 0; i < N; i++)
        {
            sampleVec.emplace_back(cos(k * 2 * M_PI * i / N));
        }
        return sampleVec;
    }

    vector<FFT::FloatType> initSine(FFT::FloatType k)
    {
        sampleVec.clear();
        int const N = 1'024;
        for (int i = 0; i < N; i++)
        {
            sampleVec.emplace_back(sin(k * 2 * M_PI * i / N));
        }
        return sampleVec;
    }

    vector<FFT::FloatType> naiveIntensity(vector<FFT::FloatType> const& samples) const
    {
        size_t const              n     = samples.size();
        const FFT::FloatType      twoPi = 2.0L * M_PI;
        vector<FFT::ComplexValue> spectrum(n, FFT::ComplexValue(0.0, 0.0));

        for (size_t k = 0; k < n; ++k)
        {
            auto sum = FFT::ComplexValue(0.0, 0.0);
            for (size_t m = 0; m < n; ++m)
            {
                const FFT::FloatType angle = -twoPi * FFT::FloatType(k) * FFT::FloatType(m) / FFT::FloatType(n);
                sum += samples[m] * std::polar(FFT::FloatType(1.0L), angle);
            }
            spectrum[k] = sum;
        }

        vector<FFT::FloatType> intensities(n);
        const FFT::FloatType   sqrtN = std::sqrt(static_cast<FFT::FloatType>(n));
        for (size_t k = 0; k < n; ++k)
        {
            intensities[k] = std::abs(spectrum[k]) / sqrtN;
        }

        return intensities;
    }

    FloatMatrix buildImpulseGrid(FFT::IntType height, FFT::IntType width, FFT::IntType row, FFT::IntType col) const
    {
        FloatMatrix grid(height, std::vector<FFT::FloatType>(width, FFT::FloatType(0)));
        assert(row >= 0 && row < height);
        assert(col >= 0 && col < width);
        grid[row][col] = FFT::FloatType(1);
        return grid;
    }

    FloatMatrix buildPatternGrid(FFT::IntType height, FFT::IntType width) const
    {
        FloatMatrix grid(height, std::vector<FFT::FloatType>(width, FFT::FloatType(0)));
        for (FFT::IntType row = 0; row < height; ++row)
        {
            for (FFT::IntType col = 0; col < width; ++col)
            {
                FFT::FloatType rowAngle = (row + 1) * 2.0L * M_PI / static_cast<FFT::FloatType>(height);
                FFT::FloatType colAngle = (col + 1) * 2.0L * M_PI / static_cast<FFT::FloatType>(width);
                grid[row][col]          = std::cos(rowAngle) + 0.5L * std::sin(colAngle);
            }
        }
        return grid;
    }

  protected:
    vector<FFT::FloatType> sampleVec; // NOSONAR S3656
    FFT                    fft;       // NOSONAR S3656
};

// 1.Input random data
TEST_F(FFTTest, random_data_test)
{
    auto samples = randomSample();
    fft.loadFloatVector(samples);
    auto transformed = fft.transform();
    auto intensity   = fft.intensityVector();

    EXPECT_EQ(transformed.size(), static_cast<size_t>(fft.numberOfPoints()));
    EXPECT_EQ(intensity.size(), transformed.size());
    for (auto const& value: intensity)
    {
        EXPECT_GE(value, 0.0L);
        EXPECT_FALSE(std::isnan(value));
        EXPECT_FALSE(std::isinf(value));
    }
}

// 2.Inputs are all zeros
TEST_F(FFTTest, all_zeros_test)
{
    fft.loadFloatVector(initConstant(0.0));
    auto transformed = fft.transform();
    auto intensity   = fft.intensityVector();
    for (auto const& value: intensity)
    {
        ASSERT_FLOAT_EQ(value, 0.0);
    }
}

// 3.Inputs are all ones (or some other nonzero value)
TEST_F(FFTTest, all_ones_test)
{
    fft.loadFloatVector(initConstant(1.0));
    auto transformed = fft.transform();
    auto intensity   = fft.intensityVector();
    ASSERT_GE(intensity[0], 1.0);
    for (size_t i = 1UL; i < intensity.size(); i++)
    {
        ASSERT_FLOAT_EQ(intensity[i], 0.0);
    }
}

// 4.Inputs alternate between +1 and -1.
TEST_F(FFTTest, alternate_1_and_minus_1_test)
{
    fft.loadFloatVector(initAlternate(1.0, -1.0));
    auto transformed = fft.transform();
    auto intensity   = fft.intensityVector();
    ASSERT_GT(intensity[512], 0.0);
    for (size_t i = 0; i < intensity.size(); i++)
    {
        if (i == 512)
        {
            continue;
        }
        ASSERT_FLOAT_EQ(intensity[i], 0.0);
    }
}

// 5.Input is e^(8*j*2*pi*i/N) for i = 0,1,2, ...,N-1. (j = sqrt(-1))
TEST_F(FFTTest, e_8_test)
{
    auto samples = initSine(8.0);
    fft.loadFloatVector(samples);
    auto transformed = fft.transform();

    auto const           numPoints         = static_cast<int>(fft.numberOfPoints());
    auto const           positiveIndex     = 8;
    auto const           negativeIndex     = numPoints - positiveIndex;
    const FFT::FloatType halfPoints        = static_cast<FFT::FloatType>(numPoints) / 2.0L;
    const FFT::FloatType expectedIntensity = std::sqrt(static_cast<FFT::FloatType>(numPoints)) / 2.0L;

    EXPECT_NEAR(fft.realAt(positiveIndex), 0.0L, 1e-9L);
    EXPECT_NEAR(fft.realAt(negativeIndex), 0.0L, 1e-9L);

    EXPECT_NEAR(fft.imagAt(positiveIndex), -halfPoints, 1e-9L);
    EXPECT_NEAR(fft.imagAt(negativeIndex), halfPoints, 1e-9L);

    EXPECT_NEAR(fft.getIntensityAt(positiveIndex), expectedIntensity, 1e-9L);
    EXPECT_NEAR(fft.getIntensityAt(negativeIndex), expectedIntensity, 1e-9L);

    ASSERT_FLOAT_EQ(std::abs(transformed[positiveIndex]), halfPoints);
    ASSERT_FLOAT_EQ(std::abs(transformed[negativeIndex]), halfPoints);
}

// 6.Input is cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1.
TEST_F(FFTTest, cos_8_test)
{
    fft.loadFloatVector(initCosine(8.0));
    auto transformed = fft.transform();
    auto intensity   = fft.intensityVector();

    ASSERT_GT(intensity[8], 0.0);
    auto const peak_intensity = intensity[8];

    for (size_t i = 0; i < intensity.size(); i++)
    {
        if (i == 8 || i == 1'024 - 8)
        {
            continue;
        }
        ASSERT_LT(intensity[i], peak_intensity);
    }
}

// 7.Input is e^((43/7)*j*2*pi*i/N) for i = 0,1,2, ...,N-1. (j = sqrt(-1))
TEST_F(FFTTest, e_43_7th_test)
{
    auto samples = initSine(43.0 / 7.0);
    fft.loadFloatVector(samples);
    auto transformed = fft.transform();
    auto intensity   = fft.intensityVector();

    auto const   numPoints   = static_cast<size_t>(fft.numberOfPoints());
    size_t const positiveBin = 6;
    size_t const mirrorBin   = numPoints - positiveBin;
    auto const   max_it      = std::max_element(intensity.begin(), intensity.begin() + numPoints / 2);
    size_t const max_index   = std::distance(intensity.begin(), max_it);

    EXPECT_EQ(max_index, positiveBin);
    EXPECT_GT(intensity[positiveBin], intensity[positiveBin - 1]);
    EXPECT_GT(intensity[positiveBin], intensity[positiveBin + 1]);
    EXPECT_NEAR(intensity[positiveBin], intensity[mirrorBin], 1e-9L);

    ASSERT_LT(transformed[positiveBin].imag(), 0.0L);
    ASSERT_GT(transformed[mirrorBin].imag(), 0.0L);

    auto conjValue = std::conj(transformed[positiveBin]);
    EXPECT_NEAR(transformed[mirrorBin].real(), conjValue.real(), 1e-9L);
    EXPECT_NEAR(transformed[mirrorBin].imag(), conjValue.imag(), 1e-9L);

    auto const neighbor_peak = std::max(intensity[positiveBin - 1], intensity[positiveBin + 1]);
    EXPECT_GT(intensity[positiveBin], 5.0L * neighbor_peak);
}

// 8.Input is cos((43/7)*2*pi*i/N) for i = 0,1,2, ...,N-1.
TEST_F(FFTTest, cos_43_7th_test)
{
    fft.loadFloatVector(initCosine(43.0 / 7.0));
    auto transformed = fft.transform();
    auto intensity   = fft.intensityVector();

    auto   max_it    = std::max_element(intensity.begin(), intensity.begin() + 512);
    size_t max_index = std::distance(intensity.begin(), max_it);

    EXPECT_EQ(max_index, 6);
}

// inverse FFT: If the FFT implementation provides an inverse transform (IFFT), it can be tested by applying the IFFT to
// the output of the FFT and verifying that the original input data is recovered (within numerical precision limits).
TEST_F(FFTTest, inverse_fft_recovers_original)
{
    auto original = randomSample();
    fft.loadFloatVector(original);
    auto transformed = fft.transform();

    auto const numPoints = static_cast<int>(fft.numberOfPoints());
    ASSERT_EQ(static_cast<size_t>(numPoints), original.size());

    vector<FFT::FloatType> reconstructed(numPoints, FFT::FloatType(0));
    const FFT::FloatType   twoPi = 2.0L * M_PI;

    for (int n = 0; n < numPoints; ++n)
    {
        auto sum = FFT::ComplexValue(0.0, 0.0);

        for (int k = 0; k < numPoints; ++k)
        {
            const FFT::FloatType angle = twoPi * FFT::FloatType(k) * FFT::FloatType(n) / FFT::FloatType(numPoints);
            sum += transformed[k] * std::polar(FFT::FloatType(1.0L), angle);
        }

        sum /= static_cast<FFT::FloatType>(numPoints);
        reconstructed[n] = sum.real();
    }

    for (size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_NEAR(original[i], reconstructed[i], 1e-5L);
    }
}

//
// B.Multi FFT tests - run continuous sets of random data
// 1.Data sets start at times 0, N, 2N, 3N, 4N, ....
// 2.Data sets start at times 0, N+1, 2N+2, 3N+3, 4N+4, ....
TEST_F(FFTTest, multi_test)
{
    constexpr size_t       windowSize = 1'024;
    constexpr size_t       repeats    = 5;
    vector<FFT::FloatType> longSignal(windowSize * repeats + repeats, 0.0);

    for (size_t i = 0; i < longSignal.size(); ++i)
    {
        auto const t = static_cast<FFT::FloatType>(i);
        longSignal[i] =
            std::cos(2.0L * M_PI * 4.0L * t / windowSize) + 0.5L * std::sin(2.0L * M_PI * 9.0L * t / windowSize);
    }

    auto verifyWindow = [&](size_t startIndex) {
        ASSERT_LE(startIndex + windowSize, longSignal.size());
        vector<FFT::FloatType> chunk(longSignal.begin() + startIndex, longSignal.begin() + startIndex + windowSize);
        auto                   expected = naiveIntensity(chunk);

        fft.loadFloatVector(chunk);
        fft.transform();
        auto intensity = fft.intensityVector();

        for (size_t bin = 0; bin < windowSize; ++bin)
        {
            EXPECT_NEAR(intensity[bin], expected[bin], 1e-6L) << "bin=" << bin << " start=" << startIndex;
        }
    };

    for (size_t start = 0; start + windowSize <= windowSize * repeats; start += windowSize)
    {
        verifyWindow(start);
    }

    for (size_t period = 0; period < repeats - 1; ++period)
    {
        verifyWindow(period * windowSize + period);
    }
}

//
// Linearity: The DFT (along with its other cousin transforms in the Fourier
// family) is a linear operator, so for all values of a1,a2,x1[n],x2[n]
//
//    , the following equation must hold:
//
// FFT(a1x1[n]+a2x2[n])=a1FFT(x1[n])+a2FFT(x2[n])
//
TEST_F(FFTTest, linearity_test)
{
    FFT            fft1;
    FFT            fft2;
    FFT            fft_combined;
    auto           x1 = initCosine(3.0);
    auto           x2 = initCosine(7.0);
    FFT::FloatType a1 = 0.5;
    FFT::FloatType a2 = 2.0;

    vector<FFT::FloatType> x_combined_vec;
    for (size_t i = 0; i < x1.size(); ++i)
    {
        x_combined_vec.push_back(a1 * x1[i] + a2 * x2[i]);
    }

    fft_combined.loadFloatVector(x_combined_vec);
    auto result_combined = fft_combined.transform();

    fft1.loadFloatVector(x1);
    auto result1 = fft1.transform();

    fft2.loadFloatVector(x2);
    auto result2 = fft2.transform();

    for (size_t i = 0; i < result_combined.size(); ++i)
    {
        auto linear_comb = a1 * result1[i] + a2 * result2[i];
        EXPECT_NEAR(result_combined[i].real(), linear_comb.real(), 1e-9);
        EXPECT_NEAR(result_combined[i].imag(), linear_comb.imag(), 1e-9);
    }
}

// DFT of the unit impulse: A time-domain signal equal to the Kronecker delta
// function is applied to the input of the FFT algorithm and the output is checked
// against the known DFT of the unit impulse function (it transforms to a constant
// value in all output bins). If the FFT algorithm provides an IFFT, it can be
// tested in reverse to show that it yields the unit impulse function again.
TEST_F(FFTTest, impulse_test)
{
    fft.loadFloatVector(initImpulse(0));
    auto transformed = fft.transform();
    auto intensity   = fft.intensityVector();
    ASSERT_GT(intensity[0], FFT::FloatType{});
    auto first_intensity = intensity[0];
    for (size_t i = 1; i < intensity.size(); i++)
    {
        ASSERT_FLOAT_EQ(intensity[i], first_intensity);
    }
}

//    Time shift: Two sets of data are applied to the input of the FFT algorithm;
// the only difference between the two in the time domain is a constant time shift.
// Based on the known properties of the DFT, this should effect a known linear
// phase shift between the two signals' frequency domain representations, where
// the slope of the phase shift is proportional to the time shift.
// A.Single FFT tests - N inputs and N outputs
TEST_F(FFTTest, time_shift_test)
{
    FFT  fft1;
    FFT  fft2;
    auto x         = initCosine(5.0);
    auto x_shifted = circularShift(x, 10);

    fft1.loadFloatVector(x);
    fft1.transform();
    auto intensity1 = fft1.intensityVector();

    fft2.loadFloatVector(x_shifted);
    fft2.transform();
    auto intensity2 = fft2.intensityVector();

    for (size_t i = 0; i < intensity1.size(); ++i)
    {
        EXPECT_NEAR(intensity1[i], intensity2[i], 1e-9);
    }
}

TEST_F(FFTTest, fft2d_inverse_roundtrip)
{
    util::FFT2D fft2d(3, 3);
    auto        grid          = buildPatternGrid(fft2d.height(), fft2d.width());
    auto        spectrum      = fft2d.transform(grid);
    auto        reconstructed = fft2d.transform(spectrum, true);

    for (size_t row = 0; row < spectrum.size(); ++row)
    {
        for (size_t col = 0; col < spectrum[row].size(); ++col)
        {
            EXPECT_NEAR(grid[row][col], reconstructed[row][col].real(), 1e-9L) << "row=" << row << " col=" << col;
            EXPECT_NEAR(0.0L, reconstructed[row][col].imag(), 1e-9L) << "row=" << row << " col=" << col;
        }
    }
}

TEST_F(FFTTest, fft2d_impulse_constant_spectrum)
{
    util::FFT2D fft2d(2, 2);
    auto        impulse  = buildImpulseGrid(fft2d.height(), fft2d.width(), 0, 0);
    auto        spectrum = fft2d.transform(impulse);

    auto const rows = std::views::iota(size_t{0}, spectrum.size());
    auto const cols = std::views::iota(size_t{0}, spectrum[0].size());
    for (auto const& [row, col]: std::views::cartesian_product(rows, cols))
    {
        EXPECT_NEAR(1.0L, spectrum[row][col].real(), 1e-12L);
        EXPECT_NEAR(0.0L, spectrum[row][col].imag(), 1e-12L);
    }
}
