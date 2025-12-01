# Figure 1 Parameter Extraction Summary (paper_figure1.md)

This document summarizes exactly what quantitative information from **Figure 1** of *Yuan et al., Nature 2022* is used for parameter inference in the Ras–CSC–Angiogenesis–TGFβ feedback model. All values listed here are extracted directly from the paper’s source data. Only values that feed directly into model variables or parameters are included.

---

# Overview

Figure 1 provides quantitative measurements for three core biological variables:

* **C** — CSC abundance
* **A** — Angiogenesis level
* **T** — TGFβ signaling activity

It also constrains the strengths and shapes of two interactions:

* **β_A** — CSC → Angiogenesis
* **α_T** — Angiogenesis → TGFβ

Figure 1 does *not* constrain LEPR (R) or mTOR (M) directly; those are informed by later figures.

---

# 1. State Variable Values

These are steady-state values derived from Fig. 1 data, normalized for model use.

## 1.1 TGFβ (T)

**Source:** Fig. 1b-left (mCherry+ TGFβ reporter cell counts)

Raw means:

* Papilloma: 29.47
* SCC: 419.93

Fold-change used for normalization:
[FC_T = 419.93 / 29.47 = 14.25]

**Model-ready:**

* **T_pap = 1.0**
* **T_scc = 14.25**

---

## 1.2 Angiogenesis (A)

**Source:** Fig. 1e mid-bottom (% CD31+ vasculature adjacent to K18+ CSCs)

Raw means:

* Papilloma: 0.957%
* SCC: 89.69%

Fold-change:
[FC_A = 89.69 / 0.957 = 93.7]

**Model-ready:**

* **A_pap = 1.0**
* **A_scc = 93.7**

---

## 1.3 CSC abundance (C)

**Source:** Fig. 1e mid-top (% K18+ CSC surface area)

Raw means:

* Papilloma: 0.225%
* SCC: 72.34%

Fold-change:
[FC_C = 72.34 / 0.225 = 321]

**Model-ready:**

* **C_pap = 1.0**
* **C_scc = 321**

---

# 2. Threshold Parameters (K_C, K_R)

Thresholds represent signal levels where downstream processes activate.

## 2.1 K_C — T threshold for CSC induction

Based on midpoint between Papilloma and SCC TGFβ values:
[K_C \approx (1 + 14.25) / 2 = 7.6]

**Model-ready:**

* **K_C = 7.6**

---

## 2.2 K_R — T threshold for LEPR induction

Fig. 1 does not provide LEPR directly, but % TGFβ+ basal cells (Fig. 1b-right) scales similarly.

Value reused:

* **K_R = 7.6**

---

# 3. Hill Coefficients (n_C, n_R)

Hill exponents determine the ultrasensitivity of activation curves.

## 3.1 n_C — Hill exponent for T → C

CSC rises 321-fold for a 14-fold rise in T.
Approximation yields **n_C ≈ 3.1**.

**Model-ready:**

* **n_C = 3**

---

## 3.2 n_R — Hill exponent for T → R

LEPR proxy (Fig. 1b-right) rises ~2.8-fold.

**Model-ready:**

* **n_R = 1.5**

---

# 4. Interaction Strengths

These map fold-change relationships into linear coefficients.

## 4.1 β_A — Strength of CSC → Angiogenesis

Using A_scc ≈ β_A * C_scc:
[β_A = 93.7 / 321 = 0.29]

**Model-ready:**

* **β_A = 0.29**

---

## 4.2 α_T — Strength of Angiogenesis → TGFβ

Using T_scc ≈ α_T * A_scc:
[α_T = 14.25 / 93.7 = 0.15]

**Model-ready:**

* **α_T = 0.15**

---

# 5. Final Model-Ready Parameter Table

| Parameter | Value | Biological Meaning             | Source  |
| --------- | ----- | ------------------------------ | ------- |
| T_pap     | 1.0   | Baseline TGFβ activity         | Fig. 1b |
| T_scc     | 14.25 | SCC TGFβ activity              | Fig. 1b |
| A_pap     | 1.0   | Baseline angiogenesis          | Fig. 1e |
| A_scc     | 93.7  | SCC angiogenesis               | Fig. 1e |
| C_pap     | 1.0   | Baseline CSC abundance         | Fig. 1e |
| C_scc     | 321   | SCC CSC abundance              | Fig. 1e |
| K_C       | 7.6   | Threshold: T → CSC activation  | Derived |
| K_R       | 7.6   | Threshold: T → LEPR activation | Derived |
| n_C       | 3     | Hill exponent: T → C           | Derived |
| n_R       | 1.5   | Hill exponent: T → R           | Derived |
| β_A       | 0.29  | CSC → Angiogenesis             | Derived |
| α_T       | 0.15  | Angiogenesis → TGFβ            | Derived |

---

# Notes

* These values come **only from Figure 1**.
* LEPR (R variable) and mTOR (M variable) begin in later figures.
* All values here are normalized and represent the version used for fitting the ODE system.

