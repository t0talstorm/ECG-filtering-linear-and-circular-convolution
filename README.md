# Experiment 4 — ECG Filtering: Linear vs. Circular Convolution

**Course:** Digital Signal Processing Lab (DSPL) | **Semester:** 5
**Author:** Aditya Sinha (PRN: 24070123006) | **MATLAB Version:** R2026a
**Required Toolbox:** None (all operations use base MATLAB)

---

## Table of Contents

1. [Overview](#overview)
2. [Aim & Objectives](#aim--objectives)
3. [Theoretical Background](#theoretical-background)
4. [ECG Signal & Datasets](#ecg-signal--datasets)
5. [Filter Design](#filter-design)
6. [Technical Approach](#technical-approach)
7. [File Structure](#file-structure)
8. [Key Plots & Visualizations](#key-plots--visualizations)
9. [Numerical Results](#numerical-results)
10. [Key Observations](#key-observations)
11. [System Requirements & How to Run](#system-requirements--how-to-run)
12. [Key Concepts Demonstrated](#key-concepts-demonstrated)
13. [Extensions & Clinical Relevance](#extensions--clinical-relevance)

---

## Overview

This experiment provides a comprehensive comparison of **linear convolution** and **circular convolution** in the context of **ECG (electrocardiogram) signal filtering** — a real-world healthcare application. The central question is: *Which convolution method is appropriate for filtering a finite-length biomedical signal, and why?*

Two sub-experiments are covered:

1. **Synthetic ECG approach** (`Updated_Exp_04.m`): A mathematically synthesized ECG is generated using Gaussian-kernel PQRST wave models with added 50 Hz power-line interference, then filtered.
2. **Real MIT-BIH ECG approach** (`untitled4.m`): Actual ECG data from the MIT-BIH Arrhythmia Database (Record 114) is read from its binary `.dat` file, decoded from Format 212, and filtered.

Both scripts demonstrate that **linear convolution is the correct method** for filtering finite-duration signals, while circular convolution introduces a dangerous **wrap-around artifact** at the beginning of the record that corrupts medical data.

---

## Aim & Objectives

- **Primary Aim:** To compare linear and circular convolution for ECG filtering and understand the practical implications of each.
- **Objectives:**
  1. Synthesize a realistic ECG signal using the Gaussian PQRST model.
  2. Add 50 Hz power-line interference (most common ECG noise source in India/Europe).
  3. Design and apply a 5-tap (or 30-tap for MIT-BIH) moving-average FIR filter.
  4. Implement circular convolution via FFT (the standard DFT-based method).
  5. Implement zero-padded FFT convolution and verify it equals linear convolution.
  6. Quantify and visualize the **circular wrap-around artifact** vs. the correct linear output.
  7. Perform basic R-peak detection and heart rate estimation on the filtered signal.

---

## Theoretical Background

### Linear Convolution

Linear convolution of a length-N signal `x[n]` with a length-M filter `h[n]` produces:

```
y_linear[n] = SUM_k( x[k] * h[n-k] )
Output length = N + M - 1
```

MATLAB implementation:
```matlab
yLinearFull = conv(x, h, 'full');   % Length: N + M - 1
```

**Properties:**
- Correct mathematical filtering for finite-duration signals
- The extra `M-1` samples represent the filter''s "ringing down" after the signal ends
- No boundary artifacts at the beginning of the signal

### Circular Convolution

Circular (cyclic) convolution treats the signal as if it repeats periodically. For an N-point DFT:

```
y_circular[n] = IDFT( DFT(x) * DFT(h) )
                where both DFTs are N-point
```

MATLAB implementation:
```matlab
yCircular = real(ifft(fft(x, N) .* fft(h, N)));
```

**Critical Problem:** When the filter is shorter than the signal (M < N), circular convolution wraps the filter''s response to the END of the signal around to the BEGINNING. This creates an **artificial artifact at the start** of the filtered ECG — corrupting the first `M-1` samples. This artifact is clinically dangerous: it can mask real cardiac events or introduce phantom peaks.

### Zero-Padded FFT Convolution (Linear Convolution via FFT)

Linear convolution can also be computed efficiently in the frequency domain using zero-padding:

```matlab
Nfft = length(x) + length(h) - 1;    % Ensure no aliasing
yFFTLinear = real(ifft(fft(x, Nfft) .* fft(h, Nfft)));
```

If `Nfft >= N + M - 1`, the zero-padding prevents circular aliasing, making the FFT result **identical to direct linear convolution**. This is verified numerically in the script.

### Power-Line Interference (50 Hz)
In India and Europe, the electrical grid operates at 50 Hz. This frequency couples electromagnetically into ECG electrodes and leads as:

```
noise[n] = 0.08 * sin(2*pi*50*t)
```

The moving-average filter attenuates this since 50 Hz is far above the cardiac frequency range (0.5–40 Hz).

---

## ECG Signal & Datasets

### Sub-Experiment 1: Synthetic ECG (`Updated_Exp_04.m`)

A synthetic ECG is generated using Gaussian functions to model each cardiac wave:

| Wave | Formula | Amplitude | Timing (relative to R) | Width (sigma) |
|------|---------|-----------|------------------------|---------------|
| **P wave** | Gaussian | 0.12 mV | -0.18 s | 25 ms |
| **Q wave** | Gaussian (neg) | -0.15 mV | -0.03 s | 10 ms |
| **R wave** | Gaussian | 1.00 mV | 0 s (reference) | 12 ms |
| **S wave** | Gaussian (neg) | -0.25 mV | +0.03 s | 12 ms |
| **T wave** | Gaussian | 0.25 mV | +0.22 s | 50 ms |

**Signal parameters:**
```matlab
Fs       = 200 Hz        % Sampling frequency
duration = 3 s           % 3-second ECG recording
N        = 600 samples   % Total samples
beatTimes = [0.5, 1.3, 2.1, 2.97]  % 4 heartbeats
```

50 Hz power-line noise is added before filtering:
```matlab
noise = 0.08 * sin(2*pi*50*t);
x     = ecgClean + noise;
```

### Sub-Experiment 2: Real MIT-BIH ECG (`untitled4.m`)

| Property | Value |
|----------|-------|
| **Database** | MIT-BIH Arrhythmia Database (PhysioNet) |
| **Record** | 114 |
| **Sampling Rate** | 360 Hz (standard MIT-BIH) |
| **Format** | Binary Format 212 (12-bit two''s complement, packed) |
| **Channel used** | Channel 2 (MLII lead) |
| **Segment analyzed** | First 6.6 seconds (2376 samples) |
| **File** | `./dat/114.dat` |

#### MIT-BIH Format 212 Decoding
The binary `.dat` file stores two 12-bit samples per 3 bytes:

```matlab
A = fread(fid, [3, inf], 'uint8');   % Read 3 bytes at a time

M1H = bitand(A(2,:), 15);            % Low nibble: high bits of sample 1
M2H = bitshift(A(2,:), -4);          % High nibble: high bits of sample 2

sig1 = bitshift(M1H, 8) + A(1,:);   % 12-bit sample 1
sig2 = bitshift(M2H, 8) + A(3,:);   % 12-bit sample 2

% Convert to signed (two''s complement for 12-bit)
sig1(sig1 > 2047) = sig1(sig1 > 2047) - 4096;
sig2(sig2 > 2047) = sig2(sig2 > 2047) - 4096;
```

---

## Filter Design

### Sub-Experiment 1 (Synthetic ECG)
```matlab
M = 5;              % 5-tap filter
h = ones(1,M)/M;   % h[n] = [0.2, 0.2, 0.2, 0.2, 0.2]
```

### Sub-Experiment 2 (MIT-BIH ECG)
```matlab
M = 30;             % 30-tap filter (better performance on real ECG at 360 Hz)
h = ones(1,M)/M;   % h[n] = [1/30, 1/30, ..., 1/30]
```

The 30-tap filter is chosen for the real signal because at 360 Hz, a 5-tap filter''s cutoff (~144 Hz) is too high to effectively attenuate noise — more taps are needed for meaningful low-pass filtering.

---

## Technical Approach

### Linear Convolution
```matlab
yLinearFull = conv(x, h, 'full');       % Full response: N+M-1 samples
tLinear     = (0:length(yLinearFull)-1)/Fs;
yLinear     = yLinearFull(1:N);         % First N samples for comparison
```

### Circular Convolution via DFT
```matlab
yCircular = real(ifft(fft(x, N) .* fft(h, N)));
```
Note: `h` must be zero-padded to length N before the FFT.

For the MIT-BIH script:
```matlab
h_padded  = [h(:); zeros(N - length(h), 1)];
yCircular = real(ifft(fft(x, N) .* fft(h_padded, N)));
```

### Zero-Padded FFT (Verification)
```matlab
Nfft        = length(x) + length(h) - 1;
yFFTLinear  = real(ifft(fft(x, Nfft) .* fft(h, Nfft)));
verificationError = max(abs(yLinearFull - yFFTLinear));
% Expected: < 1e-10 (machine precision)
```

### Wrap-Around Error Quantification
```matlab
wrapError = yCircular - yLinear;
% max(abs(wrapError)) quantifies the artifact magnitude
```

### R-Peak Detection (MIT-BIH only)
```matlab
baseline     = median(yLinear);
max_val      = max(yLinear);
minPeakHeight = baseline + 0.5 * (max_val - baseline);
minPeakDist  = 0.4 * Fs;   % Minimum 0.4s between beats (max 150 BPM)

[pks, locs] = findpeaks(yLinear, 'MinPeakHeight', minPeakHeight, ...
                                  'MinPeakDistance', minPeakDist);
% Heart rate from RR intervals
rr_intervals = diff(locs) / Fs;
avg_hr       = 60 / mean(rr_intervals);
```

---

## File Structure

```
Experiment4/
├── Updated_Exp_04.m           # Synthetic ECG linear vs. circular convolution
├── untitled4.m                # Real MIT-BIH ECG + R-peak detection
├── Experiment4_final.mlx      # MATLAB Live Script (full report)
├── dat/
│   ├── 114.dat                # MIT-BIH Record 114 binary data (Format 212)
│   ├── 114.hea                # Header file (signal parameters)
│   ├── 114.atr                # Annotation file (expert-labeled beat types)
│   └── 114.xws                # MATLAB Waveform Browser workspace
└── README.md                  # This file
```

---

## Key Plots & Visualizations

### Figure 1 — Main Comparison (4 panels, 1400x850 px)

**Panel 1: Raw vs. Clean ECG (Synthetic) / Raw MIT-BIH ECG**
- Red line: Noisy ECG `x[n]` (with 50 Hz interference / raw ADC units)
- Blue line: Clean ECG (synthetic only, for reference)
- Shows the noise that must be removed

**Panel 2: FIR Filter `h[n]`**
- Green filled stems: All M filter coefficients at uniform value `1/M`
- Clearly shows the flat, equal-weight design

**Panel 3: Linear Convolution Output (Full)**
- Blue line: `y_linear_full` including the M-1 extra samples after signal end
- Vertical dashed line marks "End of ECG Recording"
- Shows how the filter response "rings down" after the signal ends — this is physically correct

**Panel 4: Linear vs. Circular Comparison**
- Blue solid: `y_linear[n]` — correct result
- Red dashed: `y_circular[n]` — circular convolution result
- The first M-1 samples of the circular output differ from the linear — wrap-around artifact visible

**Overall title:** "Why Linear Convolution Is Preferred for ECG Filtering"

### Figure 2 — Zoomed Wrap-Around Artifact

Zoom to first 0.15 seconds:
- Blue: Linear output — starts cleanly from zero (correct)
- Red dashed: Circular output — artificially elevated at start due to samples from the signal''s END wrapping around to the BEGINNING
- Clinically: this would be misinterpreted as a cardiac event

### Figure 3 (MIT-BIH only) — R-Peak Detection Overlay
- R-peaks detected on `y_linear` using dynamic thresholding
- Magenta downward triangles mark detected R-peaks on the filtered signal

---

## Numerical Results

### Synthetic ECG
```
ECG CONVOLUTION ANALYSIS
------------------------------------------
Number of ECG samples      = 600
Filter length              = 5
Linear output length       = 604
Circular output length     = 600
Zero-padded FFT error      = 0.000000000000  (machine precision)
Maximum wrap-around error  = X.XXXXXX        (non-zero at boundaries)
```

### MIT-BIH Record 114
```
ECG CONVOLUTION ANALYSIS
------------------------------------------
Number of ECG samples      = 2376
Filter length              = 30
Linear output length       = 2405
Circular output length     = 2376
Zero-padded FFT error      = ~ 0.0           (machine precision)
Maximum wrap-around error  = ~ non-zero

FEATURE EXTRACTION (EDGE NODE SIMULATION)
------------------------------------------
Detected R-Peaks           = ~X
Average RR Interval        = ~X.XX seconds
Estimated Heart Rate       = ~XX.X BPM
```

---

## Key Observations

| Observation | Significance |
|-------------|-------------|
| Zero-padded FFT error < 1e-10 | DFT-based linear convolution is equivalent to `conv()` |
| Circular wrap-around error > 0 | Circular DFT method is NOT equivalent for filtering |
| First M-1 samples of circular output are wrong | Would corrupt the start of an ECG recording |
| Linear output has M-1 extra samples | Physically correct — the filter continues to ring down |
| R-peak detection works on linear output | Shows filtered signal is clinically usable |
| 30-tap filter better for 360 Hz real data | Filter length must scale with sampling rate for effectiveness |

---

## System Requirements & How to Run

| Requirement | Details |
|-------------|---------|
| **MATLAB Version** | R2026a |
| **Toolboxes** | None for convolution; Signal Processing Toolbox optional (for `findpeaks`) |
| **Data** | `dat/114.dat` must exist for `untitled4.m` |

### Synthetic ECG Script
```matlab
cd 'path/to/Experiment4'
run Updated_Exp_04.m
% Output: 2 figure windows (comparison + zoomed artifact)
```

### Real MIT-BIH Script
```matlab
cd 'path/to/Experiment4'
run untitled4.m
% Output: 2 figure windows + heart rate console output
```

---

## Key Concepts Demonstrated

### 1. DFT Circular Convolution Theorem
The DFT of a convolution is the product of DFTs only for **circular** (periodic) convolution. To use the DFT for **linear** (aperiodic) convolution, zero-padding to `N + M - 1` is mandatory.

### 2. Wrap-Around Artifact
Without zero-padding, the last M-1 samples of the signal wrap around and appear at the beginning of the circular convolution output. For real-time medical signals, this corrupts the beginning of every processing block.

### 3. Overlap-Add / Overlap-Save Methods
In practice, long signals are filtered in blocks using Overlap-Add or Overlap-Save methods, which split the signal into overlapping blocks of appropriate length to avoid wrap-around while using efficient FFT computation.

### 4. MIT-BIH Format 212 Binary Decoding
Medical databases store ECG in compact binary formats. Format 212 packs two 12-bit samples into 3 bytes — a standard technique to minimize storage. Decoding it requires careful bit manipulation.

### 5. Heart Rate Estimation from R-Peak Detection
The RR interval (time between consecutive R-peaks) directly gives heart rate:
```
HR (BPM) = 60 / RR_interval (seconds)
```
This is the basis of clinical HRV (Heart Rate Variability) analysis.

### 6. Clinical Implications
Using circular convolution in medical devices (ECG monitors, Holter recorders) without proper zero-padding would introduce false signals at segment boundaries — potentially causing misdiagnosis of arrhythmias. This experiment motivates why signal processing engineers must understand these mathematical distinctions.

---

## Extensions & Clinical Relevance

1. **Overlap-Add filtering:** Filter the full MIT-BIH record (30 minutes) in real time using block processing
2. **Butterworth bandpass filter:** Design a proper 0.5–40 Hz bandpass for clinical ECG
3. **Pan-Tompkins algorithm:** Implement the gold-standard QRS detection algorithm
4. **HRV analysis:** Compute time-domain (SDNN, RMSSD) and frequency-domain HRV metrics
5. **All three MIT-BIH CSV files:** The experiment folder contains records 100, 114, and 200 for multi-record comparison

---

**Last Updated:** September 2026
*Experiment submitted as part of DSPL coursework, Semester 5.*
