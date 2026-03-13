/*
 * File Name:   hough_transform.h
 * Description: Compact line-detection helper using a discretized Hough accumulator.
 *
 * Copyright (C) 2026 Dieter J Kybelksties
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
 * @author: Codex
 */

#ifndef NS_UTIL_HOUGH_TRANSFORM_H_INCLUDED
#define NS_UTIL_HOUGH_TRANSFORM_H_INCLUDED

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>
#include <vector>

namespace util
{
class HoughTransform
{
  public:
    using IntType     = int;
    using FloatType   = double;
    using Accumulator = std::vector<std::vector<IntType>>;

    struct Line
    {
        FloatType rho;
        FloatType theta; // radians
        IntType   votes;
    };

    /**
     * @brief Build an accumulator for lines on a width×height grid.
     * @param width Number of columns in the source image.
     * @param height Number of rows in the source image.
     * @param rhoStep Distance quantization for the rho axis.
     * @param thetaBins Number of discrete theta values (0..π).
     */
    HoughTransform(IntType width, IntType height, FloatType rhoStep = 1.0L, IntType thetaBins = 180)
        : width_(width)
        , height_(height)
        , rhoStep_(rhoStep > 0 ? rhoStep : 1.0L)
        , thetaBins_(thetaBins > 0 ? thetaBins : 180)
        , diag_(std::hypot(static_cast<FloatType>(width_ - 1), static_cast<FloatType>(height_ - 1)))
    {
        assert(width_ > 0 && height_ > 0);
        assert(thetaBins_ > 0);

        FloatType maxRho   = diag_ + rhoStep_ * 0.5L;
        FloatType rhoRange = maxRho * 2.0L;
        rhoBins_           = static_cast<IntType>(std::ceil(rhoRange / rhoStep_));
        if (rhoBins_ <= 0)
        {
            rhoBins_ = 1;
        }
        halfRhoRange_ = (static_cast<FloatType>(rhoBins_) * rhoStep_) / 2.0L;
        thetaStep_    = std::numbers::pi_v<FloatType> / static_cast<FloatType>(thetaBins_);

        accumulator_.assign(rhoBins_, std::vector<IntType>(thetaBins_, 0));
        cosTheta_.resize(thetaBins_);
        sinTheta_.resize(thetaBins_);

        for (IntType t = 0; t < thetaBins_; ++t)
        {
            FloatType theta = static_cast<FloatType>(t) * thetaStep_;
            cosTheta_[t]    = std::cos(theta);
            sinTheta_[t]    = std::sin(theta);
        }
    }

    IntType width() const
    {
        return width_;
    }

    IntType height() const
    {
        return height_;
    }

    IntType rhoBins() const
    {
        return rhoBins_;
    }

    IntType thetaBins() const
    {
        return thetaBins_;
    }

    FloatType thetaStep() const
    {
        return thetaStep_;
    }

    Accumulator const &accumulator() const
    {
        return accumulator_;
    }

    void reset()
    {
        for (auto &row: accumulator_)
        {
            std::fill(row.begin(), row.end(), IntType(0));
        }
    }

    /**
     * @brief Accumulate votes from a dense boolean/image grid.
     * @param input Rows are expected in [0,height) and columns in [0,width).
     * @param threshold Pixel values greater than or equal to `threshold` are interpreted as edges.
     */
    template <typename PixelType>
    void accumulateEdges(std::vector<std::vector<PixelType>> const &input, PixelType threshold = PixelType(1))
    {
        assert(static_cast<IntType>(input.size()) == height_);

        for (IntType row = 0; row < height_; ++row)
        {
            auto const &scan = input[row];
            assert(static_cast<IntType>(scan.size()) == width_);

            for (IntType col = 0; col < width_; ++col)
            {
                if (scan[col] < threshold)
                {
                    continue;
                }

                accumulatePoint(col, row);
            }
        }
    }

    /**
     * @brief Accumulate votes for a sparse list of edge coordinates.
     */
    void accumulateEdgePoints(std::vector<std::pair<IntType, IntType>> const &points)
    {
        for (auto const &item: points)
        {
            assert(item.first >= 0 && item.first < width_);
            assert(item.second >= 0 && item.second < height_);
            accumulatePoint(item.first, item.second);
        }
    }

    /**
     * @brief Read back the most prominent lines from the accumulator.
     * @param maxLines Maximum number of peaks to return.
     * @param minVotes Minimum votes required for a line to appear in the result.
     */
    std::vector<Line> detectLines(std::size_t maxLines = 5, IntType minVotes = 1) const
    {
        std::vector<Line> lines;
        lines.reserve(std::min(static_cast<std::size_t>(rhoBins_ * thetaBins_), maxLines));

        for (IntType rhoIdx = 0; rhoIdx < rhoBins_; ++rhoIdx)
        {
            for (IntType thetaIdx = 0; thetaIdx < thetaBins_; ++thetaIdx)
            {
                IntType votes = accumulator_[rhoIdx][thetaIdx];
                if (votes < minVotes)
                {
                    continue;
                }

                lines.push_back(Line{rhoFromIndex(rhoIdx), static_cast<FloatType>(thetaIdx) * thetaStep_, votes});
            }
        }

        std::sort(lines.begin(), lines.end(), [](auto const &lhs, auto const &rhs) { return lhs.votes > rhs.votes; });

        if (lines.size() > maxLines)
        {
            lines.resize(maxLines);
        }

        return lines;
    }

  private:
    void accumulatePoint(IntType x, IntType y)
    {
        FloatType fx = static_cast<FloatType>(x);
        FloatType fy = static_cast<FloatType>(y);

        for (IntType thetaIdx = 0; thetaIdx < thetaBins_; ++thetaIdx)
        {
            FloatType rho        = fx * cosTheta_[thetaIdx] + fy * sinTheta_[thetaIdx];
            FloatType normalized = (rho + halfRhoRange_) / rhoStep_;
            IntType   rhoIdx     = static_cast<IntType>(std::round(normalized));

            if (rhoIdx < 0 || rhoIdx >= rhoBins_)
            {
                continue;
            }

            ++accumulator_[rhoIdx][thetaIdx];
        }
    }

    FloatType rhoFromIndex(IntType rhoIndex) const
    {
        return (static_cast<FloatType>(rhoIndex) * rhoStep_) - halfRhoRange_;
    }

    IntType                width_;
    IntType                height_;
    FloatType              rhoStep_;
    IntType                thetaBins_;
    IntType                rhoBins_;
    FloatType              halfRhoRange_;
    FloatType              thetaStep_;
    FloatType              diag_;
    Accumulator            accumulator_;
    std::vector<FloatType> cosTheta_;
    std::vector<FloatType> sinTheta_;
};
} // namespace util

#endif // NS_UTIL_HOUGH_TRANSFORM_H_INCLUDED
