/*
 * File Name:   fft.h
 * Description: Fast Fourier implementation
 * Copyright (C) 2026 Dieter J Kybelksties <github@kybelksties.com>
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
 * @date: 2026-03-23
 * @author: Dieter J Kybelksties
 */

#ifndef NS_UTIL_FFT_H_INCLUDED
#define NS_UTIL_FFT_H_INCLUDED

#include "constants.h"

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

namespace util
{
/**
 * @brief Fast Fourier transform implementation.
 *
 * Performs a radix-2 Cooley–Tukey transform over a fixed-size buffer,
 * stores the complex spectrum, and exposes magnitude/phase helpers.
 */
class FFT
{
  public:
    using IntType   = int64_t;
    using FloatType = long double;

    using IntVector     = std::vector<IntType>;
    using FloatVector   = std::vector<FloatType>;
    using ComplexValue  = std::complex<FloatType>;
    using ComplexVector = std::vector<ComplexValue>;
    using ComplexMatrix = std::vector<ComplexVector>;

  private:
    constexpr static FloatType PI = std::numbers::pi_v<FloatType>;

  public:
    /**
     * @brief Construct an FFT instance with precomputed twiddle factors.
     * @param logOfNumOfPoints log2 of the FFT size (must reference a power of two entry in powTwo).
     * @param sampleRate Sample rate in Hz used for frequency helper conversions.
     * @param calibrate If true, initialize the recording tape with a 1 kHz calibration tone.
     */
    explicit FFT(IntType logOfNumOfPoints = 10, IntType sampleRate = 1'024, bool calibrate = false)
        : logOfPoints_(logOfNumOfPoints)
        , numOfPoints_(powTwo[logOfNumOfPoints])
        , sampleRate_(sampleRate)
    {
        tapeOfDoubles_.resize(numOfPoints_);

        if (calibrate)
        {
            // 1 kHz calibration wave
            FloatType TWO_KILO_PI        = 2.0 * FloatType(PI) * 1000.0L;
            FloatType CALIBRATION_FACTOR = 1600.0L;

            for (IntType i = 0; i < numOfPoints_; i++)
            {
                tapeOfDoubles_[i] =
                    CALIBRATION_FACTOR * std::sin((TWO_KILO_PI * i) / static_cast<FloatType>(sampleRate_));
            }
        }
        else
        {
            for (IntType i = 0; i < numOfPoints_; i++)
            {
                tapeOfDoubles_[i] = FloatType(0);
            }
        }

        sqrtOfPoints_ = std::sqrt(static_cast<FloatType>(numOfPoints_));

        bitReverseVector_.resize(numOfPoints_);
        transformedComplexVector_.resize(numOfPoints_);
        complexExpontials_.resize(logOfPoints_ + 1);

        // Pre-compute complex exponentials
        for (IntType l = 1; l <= logOfPoints_; l++)
        {
            (complexExpontials_[l]).resize(numOfPoints_);

            for (IntType i = 0; i < numOfPoints_; i++)
            {
                static FloatType const TWO_PI = 2. * PI;
                FloatType              re     = std::cos(TWO_PI * i / FloatType(powTwo[l]));
                FloatType              im     = -std::sin(TWO_PI * i / FloatType(powTwo[l]));

                complexExpontials_[l][i] = ComplexValue(re, im);
            }
        }

        // set up bit reverse mapping
        IntType rev        = 0;
        IntType halfPoints = numOfPoints_ / 2;

        for (IntType i = 0; i < numOfPoints_ - 1; i++)
        {
            bitReverseVector_[i] = rev;
            IntType mask         = halfPoints;

            // add 1 backwards
            while (rev >= mask)
            {
                rev -= mask; // turn off this bit
                mask >>= 1;
            }

            rev += mask;
        }

        bitReverseVector_[numOfPoints_ - 1] = numOfPoints_ - 1;
    }

    FFT(const FFT &lhs)                = default;
    FFT &operator=(const FFT &lhs)     = default;
    FFT(FFT &&rhs) noexcept            = default;
    FFT &operator=(FFT &&rhs) noexcept = default;

    IntType numberOfPoints() const
    {
        return numOfPoints_;
    }

    /**
     * @brief Load new samples into the recording tape and reinitialize the FFT buffer.
     * @param sampleVector Samples to append at the end of the current tape.
     */
    void loadFloatVector(FloatVector sampleVector)
    {
        IntType cSample = sampleVector.size();
        if (cSample > numOfPoints_)
        {
            return;
        }

        auto start  = tapeOfDoubles_.begin();
        auto middle = start + cSample;
        auto end    = tapeOfDoubles_.end();

        // make space for samples at the end of the tape
        // shifting previous samples towards the beginning
        std::rotate(start, middle, end);

        // copy samples from iterator to tail end of tape
        IntType iTail = numOfPoints_ - cSample;
        for (IntType i = 0; i < cSample; i++)
        {
            tapeOfDoubles_[i + iTail] = sampleVector[i];
        }
        // Initialize the FFT buffer
        for (IntType i = 0; i < numOfPoints_; i++)
        {
            set(i, tapeOfDoubles_[i]);
        }
    }

    /**
     * @brief Execute the FFT in either forward or inverse direction.
     * @param inverse when true, run the inverse FFT; false runs the forward FFT.
     * @return ComplexVector Resulting spectrum (forward) or reconstructed samples (inverse).
     */
    ComplexVector transform(bool inverse = false)
    {
        if (inverse)
        {
            return inverseTransform();
        }

        IntType step = 1;
        for (IntType level = 1; level <= logOfPoints_; level++)
        {
            IntType increm = step * 2;

            for (IntType j = 0; j < step; j++)
            {
                // U = exp ( - 2 PI j / 2 ^ level )
                ComplexValue U = complexExpontials_[level][j];

                for (IntType i = j; i < numOfPoints_; i += increm)
                {
                    // butterfly
                    ComplexValue T = U;

                    T *= transformedComplexVector_[i + step];

                    transformedComplexVector_[i + step] = transformedComplexVector_[i];
                    transformedComplexVector_[i + step] -= T;
                    transformedComplexVector_[i] += T;
                }
            }

            step *= 2;
        }

        return transformedComplexVector_;
    }

    /**
     * @brief Execute the inverse FFT to recover time-domain samples.
     *
     * The implementation performs bit-reversal on the stored spectrum, runs the
     * radix-2 algorithm with conjugated twiddle factors, and divides by N to
     * normalize the result.
     */
    ComplexVector inverseTransform() const
    {
        return inverseTransformFromSpectrum(transformedComplexVector_);
    }

    /**
     * @brief Run the inverse FFT on an arbitrary spectrum vector.
     */
    ComplexVector inverseTransformFromSpectrum(ComplexVector const &spectrum) const
    {
        return inverseTransformImpl(spectrum);
    }

    /**
     * @brief Return the normalized magnitude of each FFT bin.
     */
    FloatVector intensityVector() const
    {
        FloatVector reval;

        for (auto const &v: transformedComplexVector_)
        {
            reval.push_back(std::abs(v / sqrtOfPoints_));
        }

        return reval;
    }

    /**
     * @brief Get the magnitude of a single bin normalized by √N.
     * @param index Frequency bin index.
     */
    FloatType getIntensityAt(IntType index) const
    {
        assert(index < numOfPoints_);

        return (std::abs(transformedComplexVector_[index]) / sqrtOfPoints_);
    }

    /**
     * @brief Read the real component of a bin.
     */
    FloatType realAt(IntType index) const
    {
        assert(index < numOfPoints_);

        return (transformedComplexVector_[index].real());
    }

    /**
     * @brief Read the imaginary component of a bin.
     */
    FloatType imagAt(IntType index) const
    {
        assert(index < numOfPoints_);

        return (transformedComplexVector_[index].imag());
    }

    /**
     * @brief Convert a bin index to the corresponding frequency in Hz.
     */
    IntType getFrequencyOfSampleAt(IntType point) const
    {
        assert(point < numOfPoints_);

        // return frequency in Hz of a given point
        IntType x = sampleRate_ * point;

        return (x / numOfPoints_);
    }

    /**
     * @brief Convert a frequency in Hz to the closest FFT bin index.
     */
    IntType HzToPoint(IntType freq) const
    {
        return ((numOfPoints_ * freq) / sampleRate_);
    }

    /**
     * @brief Return the highest representable frequency (equal to the sample rate).
     */
    IntType maxFrequency() const
    {
        return sampleRate_;
    }

    /**
     * @brief Access the raw recording tape for diagnostics.
     */
    IntType tapeOfDoublesAt(IntType i) const
    {
        assert(i < numOfPoints_);

        return (static_cast<IntType>(tapeOfDoubles_[i]));
    }

  private:
    ComplexVector inverseTransformImpl(ComplexVector const &spectrum) const
    {
        ComplexVector buffer(numOfPoints_);
        assert(static_cast<IntType>(spectrum.size()) == numOfPoints_);
        for (IntType i = 0; i < numOfPoints_; i++)
        {
            buffer[bitReverseVector_[i]] = spectrum[i];
        }

        IntType step = 1;
        for (IntType level = 1; level <= logOfPoints_; level++)
        {
            IntType increm = step * 2;

            for (IntType j = 0; j < step; j++)
            {
                ComplexValue U = std::conj(complexExpontials_[level][j]);

                for (IntType i = j; i < numOfPoints_; i += increm)
                {
                    ComplexValue T    = U * buffer[i + step];
                    ComplexValue temp = buffer[i];

                    buffer[i + step] = temp - T;
                    buffer[i]        = temp + T;
                }
            }

            step *= 2;
        }

        for (auto scale = static_cast<FloatType>(numOfPoints_); auto &value: buffer)
        {
            value /= scale;
        }

        return buffer;
    }

    void set(IntType i, FloatType val)
    {
        transformedComplexVector_[bitReverseVector_[i]] = ComplexValue(val);
    }

    IntType       logOfPoints_;
    IntType       numOfPoints_;
    IntType       sampleRate_;
    FloatType     sqrtOfPoints_;
    IntVector     bitReverseVector_;         // bit reverse vector
    ComplexVector transformedComplexVector_; // in-place FFT array
    ComplexMatrix complexExpontials_;        // exponentials
    FloatVector   tapeOfDoubles_;            // recording tape
};

/*!
 * @brief 2-dimensional FFT helper that reuses the radix-2 logic from `FFT`.
 *
 * Breaks a 2-D transform into two sequences of 1-D transforms (rows first, then
 * columns) and exposes convenience helpers for real and complex grids.
 */
class FFT2D
{
  public:
    using FFT1D         = FFT;
    using IntType       = FFT1D::IntType;
    using FloatType     = FFT1D::FloatType;
    using FloatVector   = FFT1D::FloatVector;
    using FloatMatrix   = std::vector<FloatVector>;
    using ComplexValue  = FFT1D::ComplexValue;
    using ComplexVector = FFT1D::ComplexVector;
    using ComplexMatrix = FFT1D::ComplexMatrix;

  private:
    constexpr static FloatType PI = std::numbers::pi_v<FloatType>;

    struct Plan
    {
        IntType          logPoints;
        IntType          numPoints;
        ComplexMatrix    complexExponentials;
        FFT1D::IntVector bitReverse;
    };

  public:
    /**
     * @brief Prepare a 2-D transform by planning two independent radix-2 FFTs.
     * @param logWidth log2 of the number of columns.
     * @param logHeight log2 of the number of rows.
     */
    FFT2D(IntType logWidth, IntType logHeight)
        : rowPlan_(buildPlan(logWidth))
        , columnPlan_(buildPlan(logHeight))
        , width_(rowPlan_.numPoints)
        , height_(columnPlan_.numPoints)
    {
    }

    IntType width() const
    {
        return width_;
    }

    IntType height() const
    {
        return height_;
    }

    /**
     * @brief Run a separable transform on a complex matrix.
     * @param input Matrix stored in row-major order.
     * @param inverse Set to true to compute the inverse transform.
     */
    ComplexMatrix transform(ComplexMatrix const &input, bool inverse = false) const
    {
        assert(static_cast<IntType>(input.size()) == height());

        ComplexMatrix rowTransformed;
        rowTransformed.reserve(height());

        for (auto const &row: input)
        {
            assert(static_cast<IntType>(row.size()) == width());
            rowTransformed.push_back(run1DPlan(row, rowPlan_, inverse));
        }

        ComplexMatrix result(height_, ComplexVector(width_));
        ComplexVector column(height_);

        for (IntType col = 0; col < width(); col++)
        {
            for (IntType row = 0; row < height(); row++)
            {
                column[row] = rowTransformed[row][col];
            }

            ComplexVector columnTransformed = run1DPlan(column, columnPlan_, inverse);

            for (IntType row = 0; row < height(); row++)
            {
                result[row][col] = columnTransformed[row];
            }
        }

        return result;
    }

    /**
     * @brief Run a separable transform on a real-valued grid.
     */
    ComplexMatrix transform(FloatMatrix const &input, bool inverse = false) const
    {
        assert(static_cast<IntType>(input.size()) == height());

        ComplexMatrix complexInput;
        complexInput.reserve(height());

        for (auto const &row: input)
        {
            assert(static_cast<IntType>(row.size()) == width());
            ComplexVector converted(width_);

            for (IntType col = 0; col < width(); col++)
            {
                converted[col] = ComplexValue(row[col], FloatType(0));
            }

            complexInput.push_back(std::move(converted));
        }

        return transform(complexInput, inverse);
    }

  private:
    ComplexVector run1DPlan(ComplexVector const &input, Plan const &plan, bool inverse) const
    {
        ComplexVector buffer(plan.numPoints);
        for (IntType i = 0; i < plan.numPoints; i++)
        {
            buffer[plan.bitReverse[i]] = input[i];
        }

        IntType step = 1;
        for (IntType level = 1; level <= plan.logPoints; level++)
        {
            IntType increm = step * 2;

            for (IntType j = 0; j < step; j++)
            {
                ComplexValue U =
                    inverse ? std::conj(plan.complexExponentials[level][j]) : plan.complexExponentials[level][j];

                for (IntType i = j; i < plan.numPoints; i += increm)
                {
                    ComplexValue T    = U * buffer[i + step];
                    ComplexValue temp = buffer[i];

                    buffer[i + step] = temp - T;
                    buffer[i]        = temp + T;
                }
            }

            step *= 2;
        }

        if (inverse)
        {
            for (auto scale = static_cast<FloatType>(plan.numPoints); auto &value: buffer)
            {
                value /= scale;
            }
        }

        return buffer;
    }

    Plan buildPlan(IntType logPoints) const
    {
        assert(logPoints >= 0);
        constexpr IntType POW_TWO_COUNT = std::size(powTwo) / sizeof(powTwo[0]);
        assert(logPoints < POW_TWO_COUNT);

        Plan plan;
        plan.logPoints = logPoints;
        plan.numPoints = powTwo[logPoints];
        plan.bitReverse.resize(plan.numPoints);
        plan.complexExponentials.resize(logPoints + 1);

        FloatType const TWO_PI = FloatType(2.0) * PI;

        for (IntType level = 1; level <= logPoints; level++)
        {
            plan.complexExponentials[level].resize(plan.numPoints);

            for (IntType i = 0; i < plan.numPoints; i++)
            {
                FloatType re = std::cos(TWO_PI * i / static_cast<FloatType>(powTwo[level]));
                FloatType im = -std::sin(TWO_PI * i / static_cast<FloatType>(powTwo[level]));

                plan.complexExponentials[level][i] = ComplexValue(re, im);
            }
        }

        IntType rev        = 0;
        IntType halfPoints = plan.numPoints / 2;

        for (IntType i = 0; i < plan.numPoints - 1; i++)
        {
            plan.bitReverse[i] = rev;
            IntType mask       = halfPoints;

            while (rev >= mask)
            {
                rev -= mask;
                mask >>= 1;
            }

            rev += mask;
        }

        plan.bitReverse[plan.numPoints - 1] = plan.numPoints - 1;

        return plan;
    }

    Plan    rowPlan_;
    Plan    columnPlan_;
    IntType width_;
    IntType height_;
};
}; // namespace util

#endif // !defined(NS_UTIL_FFT_H_INCLUDED)
