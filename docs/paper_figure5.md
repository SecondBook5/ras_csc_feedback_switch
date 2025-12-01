# Figure 5 Parameter Extraction Summary (paper_figure5.md)

This document summarizes all model-relevant information from **Figure 5** of *Yuan et al., Nature 2022*, focusing on how **LEPR activates PI3K–AKT and mTORC1** to promote malignant SCC expansion. Figure 5 provides **pathway-level enrichment**, **biochemical activation dependencies**, and **functional inhibitor experiments**, all of which tightly constrain the structure and parameterization of the LEPR→PI3K→mTOR→growth module of the model.

---

# 1. Overview: What Figure 5 Contributes to the Model

Figure 5 establishes the following fundamental principles:

* LEPR signalling in SCC is routed through **both PI3K–AKT and mTORC1**.
* The **majority** of PI3K- and mTOR-dependent tumour growth is **LEPR-dependent**, not background noise.
* PI3K inhibition leaves ~30% of tumour growth intact; mTOR inhibition leaves ~30% intact.
* LEPR-null SCCs do not activate pAKT or pS6 in response to leptin.
* LEPR and pS6 co-localize in invasive CSC fronts.

Thus, Figure 5 defines **how LEPR funnels into intracellular effectors** and **how much of the growth phenotype relies on those effectors**.

These effect sizes are **directly translatable into ODE parameters** for the PI3K and mTOR branches.

---

# 2. KEGG Enrichment (Fig. 5b)

Both comparisons—LEPR+ CSC signatures and Lepr^ctrl vs Lepr^null tumour comparisons—show:

* **PI3K–Akt signalling** as the top KEGG pathway.
* **mTORC1 signalling** (via pS6/pS6K) as the major downstream node.

**Model implication:**
PI3K and mTOR must be included as **separate but coupled arms** downstream of R(t) (LEPR).

We define branch weights:

* `omega_PI3K = 0.60`
* `omega_mTOR = 0.40`

These reflect the relative prominence of each pathway in KEGG and in biochemical readouts.

---

# 3. Biochemical Activation (Fig. 5c,e,f,g)

The immunoblots and immunofluorescence data show:

* **pAKT activation requires LEPR** (only Lepr^ctrl cells respond to leptin).
* **pS6/pS6K (mTORC1) activation also requires LEPR**.
* In vivo, pS6 co-localizes with LEPR-reporter+ invasive CSCs.

**Model implications:**

* The activation terms for PI3K and mTOR should each be Hill functions of R(t):

  * `PI3K = eta_PI3K * R^n_PI3K / (K_PI3K^n + R^n)`
  * `mTOR = eta_M * R^n_M / (K_M^n + R^n)`
* Basal activation in R-null background must be **non-zero**, matching inhibitor experiments.

These do not give absolute values but strongly constrain *structure* and *relative dependency*.

---

# 4. PI3K Inhibition Experiment (BKM120) — Fig. 5d

Using the provided raw tumour volumes, the mean values at Day 28 are:

| Condition          | Mean Volume |
| ------------------ | ----------- |
| Leprnull + vehicle | ~76.1       |
| Leprnull + BKM120  | ~25.1       |
| Leprctrl + vehicle | ~490.9      |
| Leprctrl + BKM120  | ~137.2      |

### Inference 1: Residual growth when PI3K is blocked

Residual = inhibitor / vehicle

* LEPR-null: 25.1 / 76.1 ≈ **0.33**
* LEPR-ctrl: 137.2 / 490.9 ≈ **0.28**

**Model parameter:**

* `pi3k_residual_fraction ≈ 0.30`

Meaning: **~30% of growth is PI3K-independent.**

### Inference 2: How much PI3K-dependent growth requires LEPR

Compute the decrease:

* Ctrl drop: 490.9 − 137.2 = **353.7**
* Null drop: 76.1 − 25.1 = **51.0**

Fraction attributable to LEPR:

* 353.7 / (353.7 + 51.0) ≈ **0.87**

**Model parameter:**

* `pi3k_lepr_fraction ≈ 0.87`

Meaning: **~87% of PI3K-dependent growth requires LEPR.**

---

# 5. mTOR Inhibition Experiment (Rapamycin) — Fig. 5h

Using week-6 tumour volumes:

| Condition            | Mean Volume |
| -------------------- | ----------- |
| Leprnull + vehicle   | ~372.5      |
| Leprnull + rapamycin | ~137.4      |
| Leprctrl + vehicle   | ~802.0      |
| Leprctrl + rapamycin | ~225.1      |

### Inference 1: Residual growth when mTORC1 is blocked

Residual = inhibitor / vehicle

* Leprnull: 137.4 / 372.5 ≈ **0.37**
* Leprctrl: 225.1 / 802.0 ≈ **0.28**

**Model parameter:**

* `mtor_residual_fraction ≈ 0.32`

### Inference 2: Fraction of mTOR-sensitive growth that requires LEPR

Drops:

* Ctrl: 802.0 − 225.1 = **576.9**
* Null: 372.5 − 137.4 = **235.0**

Fraction attributable to LEPR:

* 576.9 / (576.9 + 235.0) ≈ **0.71**

**Model parameter:**

* `mtor_lepr_fraction ≈ 0.71`

Meaning: **~71% of mTOR-dependent growth requires LEPR.**

---

# 6. Final Model-Ready Parameters from Figure 5

| Parameter                | Value | Meaning                                             |
| ------------------------ | ----- | --------------------------------------------------- |
| `omega_PI3K`             | 0.60  | Fraction of LEPR effects routed through PI3K–AKT    |
| `omega_mTOR`             | 0.40  | Fraction routed through mTORC1                      |
| `pi3k_residual_fraction` | 0.30  | Growth fraction that persists under PI3K inhibition |
| `pi3k_lepr_fraction`     | 0.87  | PI3K-dependent growth requiring LEPR                |
| `mtor_residual_fraction` | 0.32  | Growth fraction that persists under rapamycin       |
| `mtor_lepr_fraction`     | 0.71  | mTOR-dependent growth requiring LEPR                |

---

# 7. Implications for ODE Architecture

Figure 5 forces the following structural decisions:

1. **Two explicit downstream nodes:** PI3K(t), M(t)
2. **Both activated by LEPR(t)** via Hill kinetics
3. **Growth term** must include:

   * LEPR-driven PI3K and mTOR contributions
   * LEPR-independent residual growth (~30%)
4. **Inhibitor simulations** (BKM120, rapamycin) must reproduce the residual/gain ratios above.



