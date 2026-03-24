![Fast Furious Transformation banner](assets/banners/fast_furious_transformation_banner.svg)

# Fast Furious Transformation

Lightweight header-only FFT toolkit for embedded-friendly C++ projects. `util::FFT` implements a Cooley–Tukey radix-2 transform, exposes helpers to read intensity, real/imaginary slices, and maps bin indices back to frequencies with optional calibration support.

## Highlights

- `include/fft.h` delivers `util::FFT` plus `util::FFT2D`, both of which reuse the Cooley–Tukey radix-2 plan for single- and two-dimensional transforms while exposing helpers such as `intensityVector()`, `realAt()`, and `imagAt()`.
- `include/spectral_filter.h` provides a frequency-domain filtering helper with low-, high-, and band-pass presets that mask the FFT spectrum before reconstructing the signal.
- `include/hough_transform.h` builds an accumulator for (rho, theta) votes, precomputes the relevant sin/cos tables, and exposes `detectLines()` so the same signal utilities can feed downstream analytics.
- Frequency helpers such as `HzToPoint()` and `getFrequencyOfSampleAt()` translate between sample-rate bins and Hertz.
- Test coverage under `test/run_tests` now encompasses FFT, spectral-filter, and Hough suites; `fft_tests.cc`, `spectral_filter_tests.cc`, and `hough_transform_tests.cc` each target canonical patterns to validate the helpers.
- `cmake-common/` contains shared company-wide build defaults; this project enforces static linking, C++23, and a clang-tidy hook.

## Repository Layout

-- `include/` — public headers (`fft.h`, `spectral_filter.h`, `hough_transform.h`, `constants.h`, etc.).
- `test/` — GoogleTest suite (`FFT_tests.cc`, `SpectralFilter_tests.cc`, `HoughTransform_tests.cc`) and the `run_tests` launcher wired to CMake/CTest.
- `cmake-common/` — shared build settings applied via `include(cmake-common/cmake/DkybBuildSettings.cmake)`.
- `assets/banners/` — project banner SVG that also appears on this README.

## Build & Install

Requires CMake `>= 3.26` and a C++23 compiler (`g++`, `clang++`, or MSVC 2022+).

```bash
git submodule update --init --recursive   # fetch cmake-common and DebugTrace
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --parallel
```

The interface library `FastFuriousTransformation` installs its headers under `include/dkyb`.

## Usage

```c++
#include <dkyb/fft.h>
#include <vector>

using namespace util;

int main()
{
    FFT fft(10, 1024);
    std::vector<FFT::FLOATTYPE> samples(1024, 0.0L);
    // fill samples with your waveform
    fft.loadFloatVector(samples);

    auto spectrum = fft.transform();
    auto intensity = fft.intensityVector();

    std::cout << "bin 5 magnitude: " << intensity[5] << "\n";
    std::cout << "bin 5 real: " << fft.realAt(5) << ", imag: " << fft.imagAt(5) << "\n";
    return 0;
}
```

```c++
#include <dkyb/spectral_filter.h>

util::SpectralFilter lowPass(10, 1024);
lowPass.configureLowPass(40.0L);

auto filtered = lowPass.apply(samples);
// filtered now contains the time-domain signal with high frequencies suppressed
```

Clang-tidy runs automatically if `$CMAKE_CXX_CLANG_TIDY` is set (see the top-level `CMakeLists.txt`).

## Testing

```bash
cmake --build build --target run_tests
ctest --test-dir build
```

`run_tests` links against GoogleTest and marks the binary as a CTest target so you can run the full suite via `ctest`.

## Powered by
Reduce the smells, keep on top of code-quality. Sonar Qube is run on every push to the `main` branch on GitHub.

[![SonarQubeCloud](assets/icons/logo-sonarqube-cloud-small.png)](https://sonarcloud.io/project/overview?id=kingkybel_TypeTraits)
