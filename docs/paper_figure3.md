# Figure 3 Parameter Extraction Summary (paper_figure3.md)

This document summarizes all quantitative constraints extracted from **Figure 3** of *Yuan et al., Nature 2022*, focused on how **LEPR (R)** enhances **CSC clonogenicity, tumour-initiating capacity, and in-vivo tumour progression**.
These measurements directly constrain the downstream LEPR→mTOR→CSC portions of the Ras–CSC feedback model.

---

# **Overview**

Figure 3 provides explicit numeric values (colony assays, limiting dilutions, tumour volumes).
These produce **hard quantitative model constraints** on:

* CSC self-renewal boost from LEPR
* Tumour-initiating cell fraction
* Multiplicative effect of LEPR signalling on tumour growth
* Necessity of *signalling-competent* LEPR (ΔSig rescue fails)

These map to the model as parameters **γ_colony**, **γ_TIC**, **γ_growth**, and **γ_R_signalling**.

---

# **1. CSC Colony Formation (Fig. 3a)**

## **Raw values (provided)**

* LEPR− colony numbers: **6, 2, 5**
  → Mean = **4.33**

* LEPR+ colony numbers: **13, 15, 11**
  → Mean = **13.0**

* **Fold-change in colony number:**
  [
  \gamma_{\text{colony_number}} = 13.0 / 4.33 \approx 3.0
  ]

## **Colony area (size)**

Means computed across all entries:

* LEPR− mean area = **0.0895 cm²**

* LEPR+ mean area = **0.1357 cm²**

* **Fold-change in colony size:**
  [
  \gamma_{\text{colony_size}} = 0.1357 / 0.0895 \approx 1.52
  ]

## **Model implication**

LEPR increases **effective CSC proliferative potential and self-renewal**.
These map to the CSC growth coefficient (e.g., scaling **k_C,M**).

## **Model-ready parameters**

* **γ_colony_number = 3.0**
* **γ_colony_size = 1.52**

---

# **2. Tumour-Initiating Cell Frequency (TIC) (Fig. 3b)**

Limiting dilution yields:

| Cells injected | Lepr^ctrl PDV | Lepr^null PDV |
| -------------- | ------------- | ------------- |
| 1×10³          | 25%           | 0%            |
| 1×10⁴          | 75%           | 25%           |
| 1×10⁵          | 100%          | 100%          |

Using a single-hit Poisson model:

[
P(\text{tumour}) = 1 - e^{-fN}
]

We obtain:

* **f_TIC_ctrl ≈ 1.6 × 10⁻⁴**
  (≈1 TIC per **6.3×10³** cells)

* **f_TIC_null ≈ 3.2 × 10⁻⁵**
  (≈1 TIC per **3.1×10⁴** cells)

## **Fold-difference in tumour initiation capacity**

[
\gamma_{\text{TIC}} = f_{\text{TIC_ctrl}} / f_{\text{TIC_null}} \approx 5.0
]

## **Model implication**

LEPR increases the probability that a CSC forms a tumour by **≈5×**.
This scales the conversion from CSC number → tumour mass in the model.

## **Model-ready parameters**

* **f_TIC_ctrl = 1.6e−4**
* **f_TIC_null = 3.2e−5**
* **γ_TIC = 5.0**

---

# **3. In-Vivo Growth Advantage (Lepr^ctrl vs Lepr^null) (Fig. 3c)**

Mean tumour volumes (averaged across the 8 curves given):

| Week | Lepr^null PDV | Lepr^ctrl PDV | Fold |
| ---- | ------------- | ------------- | ---- |
| 3    | 18.9          | 41.6          | 2.20 |
| 4    | 50.2          | 194.0         | 3.87 |
| 5    | 141.7         | 329.9         | 2.33 |

Average fold across all timepoints: **≈ 2.8**
Terminal (week 5): **2.33**

## **Model implication**

LEPR enhances CSC-driven tumour growth by **~2–3×** even after accounting for TIC differences.
This is an **amplitude multiplier**, not just a shift in initial conditions.

## **Model-ready parameter**

* **γ_growth_ctrl_null = 2.3**
  (use terminal fold as the constraint)

---

# **4. Requirement for LEPR Signalling (FL vs ΔSig) (Fig. 3d)**

Mean tumour volumes:

| Week | FL    | ΔSig  | Fold |
| ---- | ----- | ----- | ---- |
| 3    | 183.4 | 41.6  | 4.41 |
| 4    | 316.4 | 77.1  | 4.10 |
| 5    | 397.9 | 117.4 | 3.39 |

Average fold = **~4×**

## **Interpretation**

* LEPR protein alone is **insufficient**.
* **Signalling-competent LEPR** enhances tumour growth by **~4×** relative to signalling-dead LEPR.

## **Model implication**

This constrains the **R → M (LEPR → mTOR)** interaction strength.

## **Model-ready value**

* **γ_R_signalling = 4.0**

---

# **Summary Table of All Extracted Parameters (Figure 3)**

| Parameter          | Value  | Origin  | Biological meaning                                  |
| ------------------ | ------ | ------- | --------------------------------------------------- |
| γ_colony_number    | 3.0    | Fig. 3a | LEPR triples CSC colony-forming events              |
| γ_colony_size      | 1.52   | Fig. 3a | LEPR increases colony size by ~50%                  |
| f_TIC_ctrl         | 1.6e−4 | Fig. 3b | TIC frequency in LEPR-intact cells                  |
| f_TIC_null         | 3.2e−5 | Fig. 3b | TIC frequency in LEPR-KO cells                      |
| γ_TIC              | 5.0    | Fig. 3b | LEPR produces 5× greater tumour-initiating capacity |
| γ_growth_ctrl_null | 2.3    | Fig. 3c | LEPR increases macroscopic tumour growth by ~2–3×   |
| γ_R_signalling     | 4.0    | Fig. 3d | Active LEPR signalling increases tumour volume 4×   |

---

# **Model Integration Notes**

Figure 3 constrains the **downstream arm** of the Ras–CSC network:

[
T \rightarrow R \rightarrow M \rightarrow C
]

Specifically:

* **R increases both CSC intrinsic proliferation (γ_colony)**
* **R increases the probability a CSC forms a tumour (γ_TIC)**
* **R increases tumour growth amplitude (γ_growth)**
* **Only signalling-competent R activates M (γ_R_signalling)**

Together, these define the **LEPR–mTOR module**.


