/*
 * File:		HoughTransform_tests.cc
 * Description:	Unit tests that exercise the line accumulator.
 *
 * Copyright (C) 2026 Dieter J Kybelksties
 */

#include "hough_transform.h"

#include <cmath>
#include <gtest/gtest.h>
#include <utility>
#include <vector>

using namespace util;

namespace
{
using GridMatrix = std::vector<std::vector<HoughTransform::IntType>>;
} // namespace

TEST(HoughTransformTest, detects_horizontal_line)
{
    constexpr int  width  = 20;
    constexpr int  height = 20;
    HoughTransform hough(width, height, 1.0L, 180);

    GridMatrix    grid(height, std::vector<int>(width, 0));
    constexpr int targetRow = 7;
    for (int col = 0; col < width; ++col)
    {
        grid[targetRow][col] = 255;
    }

    hough.accumulateEdges(grid, 128);
    auto lines = hough.detectLines(3, width);
    ASSERT_FALSE(lines.empty());

    auto const &best = lines.front();
    EXPECT_EQ(best.votes, width);
    EXPECT_NEAR(best.theta, static_cast<double>(M_PI / 2.0), 1.5 * hough.thetaStep());
    EXPECT_NEAR(best.rho, static_cast<double>(targetRow), 1.0);
}

TEST(HoughTransformTest, detects_diagonal_sparse_line)
{
    constexpr int  size = 12;
    HoughTransform hough(size, size, 1.0L, 180);

    std::vector<std::pair<int, int>> points;
    for (int i = 0; i < size; ++i)
    {
        points.emplace_back(i, i);
    }

    hough.accumulateEdgePoints(points);
    auto lines = hough.detectLines(1, size);
    ASSERT_FALSE(lines.empty());

    auto const &best = lines.front();
    EXPECT_EQ(best.votes, size);
    EXPECT_NEAR(best.theta, static_cast<double>(3.0 * M_PI / 4.0), 3.0 * hough.thetaStep());
    EXPECT_NEAR(best.rho, 0.0, 1.0);
}

TEST(HoughTransformTest, reset_clears_accumulator)
{
    constexpr int  width  = 8;
    constexpr int  height = 8;
    HoughTransform hough(width, height, 1.0L, 90);

    GridMatrix grid(height, std::vector<int>(width, 1));
    hough.accumulateEdges(grid, 1);

    auto nonEmpty = hough.detectLines(1, 1);
    ASSERT_FALSE(nonEmpty.empty());

    hough.reset();
    auto empty = hough.detectLines(1, 1);
    EXPECT_TRUE(empty.empty());
}
