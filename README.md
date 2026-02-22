# FMCW Radar System Simulation

**Institution:** Cairo University – Faculty of Engineering  
**Course:** Signals and Systems  
**Authors:** Salah Waheed & Islam Ahmed  

---

##  Project Overview

This project presents the design and simulation of a complete **Frequency-Modulated Continuous-Wave (FMCW) radar system** implemented in **MATLAB**.

The system models an **automotive-grade mmWave radar (76–77 GHz)** capable of detecting multiple moving targets and estimating both **range** and **velocity** simultaneously using 2D signal processing techniques.

The implementation follows the full radar signal chain, including signal generation, target modeling, range estimation, Doppler processing, and adaptive detection under noisy conditions.

---

## System Specifications

The radar parameters were selected to achieve high-precision detection with a **range resolution of 15 cm**.

| Parameter | Value |
|---|---|
| Carrier Frequency ($f_c$) | 76.5 GHz |
| Bandwidth ($B$) | 1 GHz |
| Chirp Duration ($T_c$) | 3 μs |
| Sampling Frequency ($F_s$) | 2 GHz |
| Range Resolution ($\Delta R$) | 0.15 m |
| Velocity Resolution ($\Delta v$) | 0.45 m/s |
| Max Unambiguous Range ($R_{max}$) | ~450 m |
| Max Unambiguous Velocity ($V_{max}$) | ~116 m/s |

---

##  Implementation Methodology

The simulation is divided into four main stages representing the complete FMCW radar processing pipeline.

### 🔹 Task A — Signal Generation

- Generate Linear Frequency Modulated (LFM) chirp signal.
- Model multiple moving targets with configurable:
  - Initial range
  - Relative velocity
- Add **Additive White Gaussian Noise (AWGN)** to simulate realistic environments.
- Tested SNR levels:
  - 5 dB
  - 10 dB
  - 15 dB

### 🔹 Task B — Range Detection

**Processing Steps:**

1. Mixing (De-chirping)  
   - Multiply transmitted signal by conjugate of received echo.
   - Extract beat frequency ($f_b$).

2. 1D FFT Processing  
   - Apply FFT on beat signal.
   - Obtain range spectrum.

**Result:**  
Targets detected at approximately **40 m** and **120 m** with sub-centimeter error.

### 🔹 Task C — Velocity Estimation

- Analyze phase variation across **512 consecutive chirps**.
- Perform Doppler FFT along slow-time dimension.
- Estimate relative target velocities.

### 🔹 Task D — Range-Doppler Map (RDM)

- Apply full **2D FFT**:
  - Range FFT (fast time)
  - Doppler FFT (slow time)
- Generate Range-Doppler Map for joint detection.
- Normalize spectrum to **0 dB** for consistent visualization across SNR levels.

---

##  Advanced Features (Bonus)

### ✔ Cell-Averaging CFAR (CA-CFAR)

Adaptive thresholding algorithm that:

- Estimates local noise floor using training cells.
- Uses guard cells to prevent target leakage.
- Maintains constant false alarm rate under varying noise.

### ✔ Automated Peak Extraction

- Statistical local maxima detection.
- Automatic velocity estimation.
- Eliminates manual cursor-based measurement.

---

##  Simulation Results

Performance evaluated under **high-noise conditions (5 dB SNR)**.

### Range Accuracy
- Detected ranges matched ground truth with error:
  - **< 1 cm**

### Velocity Accuracy
- Estimated velocities within:
  - **0.5 – 0.7 m/s**
- Consistent with theoretical FFT resolution limits.

---

##  References

1. M. Jankiraman, *FMCW Radar Design*, Artech House, 2018.  
2. S. Rao, *Introduction to mmWave Sensing: FMCW Radars*, Texas Instruments, 2017.  
3. M. A. Richards, *Fundamentals of Radar Signal Processing*, 2nd ed., McGraw-Hill, 2014.

---

## 👨‍💻 Authors

- **Salah Waheed**  
- **Islam Ahmed**
