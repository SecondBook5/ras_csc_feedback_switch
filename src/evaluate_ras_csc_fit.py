#!/usr/bin/env python3
"""
evaluate_ras_csc_fit.py

This script connects the Ras-CSC ODE model to the RNA-seq-derived
calibration targets in:

    data/processed/rnaseq/ras_csc_calibration_targets.csv

The goals are:

1. Load the per-condition module targets:
       dataset, condition,
       C_target, A_target, T_target, M_target
   which were built from the RNA-seq module scores.

2. Define a Ras-CSC feedback model with four state variables:

       C(t): CSC-like state
       A(t): Angiogenesis / vessel supply
       T(t): TGFb / Ras signaling activity
       M(t): mTOR activity

   The equations are intentionally simple Hill-type and linear terms.
   They are meant to capture qualitative structure, not detailed kinetics.

3. For each biological condition in the CSV, assign a Ras/TGFb driving
   level u_TGFb (a scalar input) and simulate the ODE to steady state.

4. Compare the simulated steady states [C, A, T, M] with the targets
   [C_target, A_target, T_target, M_target] and compute residuals.

5. Optionally, use SciPy least_squares to adjust the model parameters
   to reduce the sum of squared residuals across all conditions.

6. Save the fitted parameters and a summary of predictions and residuals
   into a JSON file, along with information criteria (AIC, BIC),
   approximate parameter uncertainties, and local sensitivities.

Usage:

    python src/evaluate_ras_csc_fit.py

to just simulate with default parameters and print residuals, or:

    python src/evaluate_ras_csc_fit.py --fit

to attempt a least-squares parameter fit, and:

    python src/evaluate_ras_csc_fit.py --fit \
        --targets data/processed/rnaseq/ras_csc_calibration_targets.csv \
        --output-json data/processed/model_fits/model_calibration_results.json
"""

# Enable postponed evaluation of type hints so annotations are stored as strings
from __future__ import annotations

# Import dataclasses for structured parameter and target containers
from dataclasses import dataclass, asdict

# Import Path for filesystem path handling
from pathlib import Path

# Import argparse for command-line argument parsing
import argparse

# Import csv and json for reading and writing text based data files
import csv
import json

# Import typing tools for clear type hints
from typing import Dict, List, Tuple, Optional

# Import numpy for numerical work
import numpy as np

# Try to import SciPy, handling absence cleanly
try:
    # Import solve_ivp for ODE integration
    from scipy.integrate import solve_ivp
    # Import least_squares for nonlinear least-squares fitting
    from scipy.optimize import least_squares

    # Record that SciPy is available
    SCIPY_AVAILABLE: bool = True
except ImportError:
    # Record that SciPy is not available on this environment
    SCIPY_AVAILABLE = False


# Define the canonical parameter names in the order they are packed into theta
PARAM_NAMES: List[str] = [
    "k_C_base",
    "k_C_M",
    "K_C",
    "n_C",
    "k_C_decay",
    "k_A_C",
    "k_A_decay",
    "k_T_A",
    "k_T_input",
    "k_T_decay",
    "k_M_T",
    "k_M_decay",
]


@dataclass
class ConditionTarget:
    """
    Container for one row of ras_csc_calibration_targets.csv.

    Each instance represents one biological condition (for example Normal,
    Papilloma, SCC, PDV_WT, PDV_LeprKO) within a given dataset, together
    with the four module targets to be matched by the ODE model at steady
    state. This provides the direct link between omics-derived scores and
    the dynamical model output.

    The model fitting procedure operates over a collection of these rows
    and attempts to find a single parameter set that brings the steady
    states close to all of the targets in a least-squares sense.

    Attributes:
        dataset:
            Name of the dataset (for example 'Bl6', 'PAP_SCC', 'PDV').
        condition:
            Biological condition name (for example 'Normal', 'Papilloma', 'SCC').
        C_target:
            Target value for the CSC-like module (CSC_bulk).
        A_target:
            Target value for the angiogenesis module (Angio_bulk).
        T_target:
            Target value for the TGFb / Ras signaling module (TGFb_bulk).
        M_target:
            Target value for the mTOR module (mTOR_bulk).
    """

    # Store the dataset identifier for this condition
    dataset: str

    # Store the condition label for this row
    condition: str

    # Store the target value for CSC-like activity
    C_target: float

    # Store the target value for angiogenesis
    A_target: float

    # Store the target value for TGFb / Ras signaling
    T_target: float

    # Store the target value for mTOR activity
    M_target: float


@dataclass
class RasCSCParameters:
    """
    Parameter set for the Ras-CSC feedback ODE model.

    The model uses four state variables:

        C: CSC-like state
        A: Angiogenesis / vessel density
        T: TGFb / Ras signaling activity
        M: mTOR activity

    The equations are:

        dC/dt = k_C_base + k_C_M * H(M; K_C, n_C) - k_C_decay * C
        dA/dt = k_A_C * C - k_A_decay * A
        dT/dt = k_T_A * A + u_TGFb * k_T_input - k_T_decay * T
        dM/dt = k_M_T * T - k_M_decay * M

    where H(M; K_C, n_C) is a Hill function:

        H(M; K_C, n_C) = M^n_C / (K_C^n_C + M^n_C)

    and u_TGFb is an external input specifying the Ras / TGFb drive
    for each biological condition.

    The default parameter values are initial guesses chosen only to
    produce qualitatively reasonable behavior. They are refined by
    least-squares fitting against the RNA-seq-derived targets.

    Attributes:
        k_C_base:
            Baseline CSC production rate independent of mTOR.
        k_C_M:
            Maximum CSC production rate driven by mTOR via the Hill term.
        K_C:
            Half-activation constant of mTOR for CSC production.
        n_C:
            Hill coefficient for mTOR-driven CSC production.
        k_C_decay:
            Linear decay (loss) rate of the CSC-like state.
        k_A_C:
            Rate at which CSCs drive angiogenesis or vessel growth.
        k_A_decay:
            Linear decay rate of angiogenesis.
        k_T_A:
            Rate at which angiogenesis increases TGFb / Ras signaling.
        k_T_input:
            Gain from the external Ras / TGFb drive u_TGFb into T.
        k_T_decay:
            Linear decay rate of TGFb / Ras signaling activity.
        k_M_T:
            Rate at which TGFb / Ras signaling increases mTOR activity.
        k_M_decay:
            Linear decay rate of mTOR activity.
    """

    # Baseline CSC production independent of mTOR
    k_C_base: float = 0.05

    # Maximum CSC production driven by mTOR through the Hill term
    k_C_M: float = 0.9

    # Half-activation constant for mTOR control of CSC production
    K_C: float = 0.4

    # Hill coefficient for mTOR control of CSC production
    n_C: float = 3.0

    # Linear decay rate for CSC-like cells
    k_C_decay: float = 0.5

    # Rate at which CSC-like cells drive angiogenesis
    k_A_C: float = 0.7

    # Linear decay rate for angiogenesis
    k_A_decay: float = 0.5

    # Rate at which angiogenesis increases Ras / TGFb signaling
    k_T_A: float = 0.5

    # Gain from external Ras / TGFb drive into T
    k_T_input: float = 1.0

    # Linear decay rate for Ras / TGFb signaling
    k_T_decay: float = 0.5

    # Rate at which Ras / TGFb signaling increases mTOR activity
    k_M_T: float = 0.7

    # Linear decay rate for mTOR activity
    k_M_decay: float = 0.5


def load_calibration_targets(csv_path: Path) -> List[ConditionTarget]:
    """
    Load ras_csc_calibration_targets.csv into a list of ConditionTarget objects.

    This function performs strict checks on the presence and types of
    required columns so that any mismatch or malformed input is caught
    immediately. This protects the calibration step from silently
    consuming incorrect or incomplete data.

    Args:
        csv_path:
            Path to the calibration targets CSV produced by
            the RNA-seq pipeline.

    Returns:
        A list of ConditionTarget instances, one per row in the CSV.

    Raises:
        ValueError:
            If the file is missing, if required columns are absent, or
            if any row cannot be parsed correctly.
    """
    # Check that the file exists before trying to read it
    if not csv_path.is_file():
        # Raise a ValueError with a clear message if the file is missing
        raise ValueError(
            f"Calibration CSV not found at: {csv_path}. "
            f"Run the RNA-seq pipeline up to the export step before using this script."
        )

    # Initialize an empty list to store ConditionTarget objects
    targets: List[ConditionTarget] = []

    # Open the CSV file in text mode with UTF-8 encoding
    with csv_path.open("r", encoding="utf-8") as fh:
        # Create a DictReader to access columns by name
        reader = csv.DictReader(fh)

        # Define the required column names for sanity checking
        required = [
            "dataset",
            "condition",
            "C_target",
            "A_target",
            "T_target",
            "M_target",
        ]

        # Verify that all required columns are present in the CSV header
        for col in required:
            # Check for each column name explicitly
            if col not in reader.fieldnames:
                # If a column is missing, raise an informative error
                raise ValueError(
                    f"Calibration CSV is missing required column '{col}'. "
                    f"Columns found: {reader.fieldnames}"
                )

        # Iterate over each row in the CSV
        for row in reader:
            try:
                # Extract and type-cast values from the current row
                dataset = str(row["dataset"])
                condition = str(row["condition"])
                C_target = float(row["C_target"])
                A_target = float(row["A_target"])
                T_target = float(row["T_target"])
                M_target = float(row["M_target"])
            except (TypeError, ValueError) as exc:
                # Raise an error if any value cannot be parsed correctly
                raise ValueError(
                    f"Failed to parse calibration row: {row}"
                ) from exc

            # Construct a ConditionTarget instance from the parsed values
            target = ConditionTarget(
                dataset=dataset,
                condition=condition,
                C_target=C_target,
                A_target=A_target,
                T_target=T_target,
                M_target=M_target,
            )

            # Append the ConditionTarget to the list
            targets.append(target)

    # If no rows were read, raise an error to avoid silent failures
    if not targets:
        # Raise a clear error if the CSV has only a header
        raise ValueError(
            f"Calibration CSV at {csv_path} contained zero data rows."
        )

    # Return the list of ConditionTarget objects
    return targets


def default_u_tgfb_map(targets: List[ConditionTarget]) -> Dict[Tuple[str, str], float]:
    """
    Build a simple mapping from (dataset, condition) to u_TGFb input.

    This uses a heuristic classification that reflects qualitative
    Ras / TGFb activity:

        - Normal tissue: low Ras / TGFb drive.
        - Papilloma: intermediate drive.
        - SCC: high drive.
        - PDV_WT and PDV_LeprKO: high drive, so any mismatch in
          PDV_LeprKO highlights that Ras alone cannot explain the data.

    This can be refined or replaced later, for example by using actual
    TGFb module scores instead of labels.

    Args:
        targets:
            List of ConditionTarget objects that define the existing
            dataset-condition pairs.

    Returns:
        A dictionary mapping (dataset, condition) -> u_TGFb value.

    Raises:
        ValueError:
            If a condition label cannot be mapped to an input value.
    """
    # Initialize an empty dictionary for the mapping
    mapping: Dict[Tuple[str, str], float] = {}

    # Loop over each ConditionTarget to populate the mapping
    for ct in targets:
        # Convert the condition name to lowercase for classification
        cond_lower = ct.condition.lower()

        # Initialize a default u_TGFb value
        u_val: float

        # Assign u_TGFb based on the condition label
        if "normal" in cond_lower:
            # Use low input for Normal
            u_val = 0.2
        elif "pap" in cond_lower:
            # Use intermediate input for Papilloma
            u_val = 0.6
        elif "scc" in cond_lower:
            # Use high input for SCC
            u_val = 1.0
        elif "pdv_wt" in cond_lower:
            # Use high input for PDV wild-type
            u_val = 1.0
        elif "pdv_leprko" in cond_lower:
            # Use high input for LeprKO to highlight structural misfit
            u_val = 1.0
        else:
            # If condition is unknown, raise an explicit error
            raise ValueError(
                f"Do not know how to assign u_TGFb for condition '{ct.condition}'. "
                f"Update default_u_tgfb_map to handle this case."
            )

        # Store the mapping for the (dataset, condition) pair
        mapping[(ct.dataset, ct.condition)] = u_val

    # Return the fully populated mapping
    return mapping


def hill_function(x: float, K: float, n: float) -> float:
    """
    Compute a standard Hill activation function.

    H(x; K, n) = x^n / (K^n + x^n)

    This function models a saturating activation where small x values
    give negligible output and large x values approach 1. Here it is
    used to represent mTOR control of CSC production.

    Args:
        x:
            Input value, for example M (mTOR activity).
        K:
            Half-activation constant that sets the midpoint of the curve.
        n:
            Hill coefficient that sets the steepness of the transition.

    Returns:
        The Hill activation value between 0 and 1.

    Raises:
        ValueError:
            If K or n are non-positive, which would make the function
            ill defined for this use case.
    """
    # Guard against non-positive K which would break the exponentiation
    if K <= 0.0:
        # Raise a clear error if K is invalid
        raise ValueError(
            f"Hill function received non-positive K={K}. K must be positive."
        )

    # Guard against non-positive n which would make the exponent undefined
    if n <= 0.0:
        # Raise a clear error if n is invalid
        raise ValueError(
            f"Hill function received non-positive n={n}. n must be positive."
        )

    # Clamp x to be at least zero to avoid negative powers
    x_clamped = max(x, 0.0)

    # Compute x^n
    x_pow = x_clamped**n

    # Compute K^n
    K_pow = K**n

    # Compute the denominator K^n + x^n and protect against zero
    denom = K_pow + x_pow if (K_pow + x_pow) > 1e-12 else 1e-12

    # Compute the Hill function value
    value = x_pow / denom

    # Return the Hill activation
    return value


def ras_csc_ode(
    t: float,
    y: np.ndarray,
    params: RasCSCParameters,
    u_tgfb: float,
) -> np.ndarray:
    """
    Compute the time derivatives for the Ras-CSC model at a given time.

    The system is autonomous in this formulation, so time enters only
    through the standard ODE solver interface. The model uses a simple
    feedback structure where CSCs support angiogenesis, angiogenesis
    supports Ras / TGFb signaling, Ras / TGFb supports mTOR, and mTOR
    in turn boosts CSCs through a Hill-type activation.

    The state vector y is:

        y[0] = C (CSC-like state)
        y[1] = A (Angiogenesis / vessels)
        y[2] = T (TGFb / Ras signaling)
        y[3] = M (mTOR activity)

    Args:
        t:
            Time variable, present for compatibility with solve_ivp but
            not used explicitly in the right hand side.
        y:
            Current state vector [C, A, T, M].
        params:
            RasCSCParameters instance with model rate constants.
        u_tgfb:
            External Ras / TGFb drive specific to the current condition.

    Returns:
        A numpy array with the derivatives [dCdt, dAdt, dTdt, dMdt].
    """
    # Unpack state variables from y with explicit names
    C = float(y[0])
    A = float(y[1])
    T = float(y[2])
    M = float(y[3])

    # Compute the Hill activation of mTOR for CSC production
    H_M = hill_function(M, params.K_C, params.n_C)

    # Compute dC/dt using baseline plus Hill-driven term minus decay
    dCdt = params.k_C_base + params.k_C_M * H_M - params.k_C_decay * C

    # Compute dA/dt using CSC-driven angiogenesis minus decay
    dAdt = params.k_A_C * C - params.k_A_decay * A

    # Compute dT/dt using angiogenesis-driven term, external input, and decay
    dTdt = params.k_T_A * A + u_tgfb * params.k_T_input - params.k_T_decay * T

    # Compute dM/dt using TGFb / Ras-driven activation and decay
    dMdt = params.k_M_T * T - params.k_M_decay * M

    # Return derivatives as a numpy array
    return np.array([dCdt, dAdt, dTdt, dMdt], dtype=float)


def simulate_to_steady_state(
    params: RasCSCParameters,
    u_tgfb: float,
    t_max: float = 200.0,
    y0: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Simulate the Ras-CSC model for a single condition and return the final state.

    This function integrates the ODE from t = 0 to t = t_max and uses the
    final state as an approximation of the steady state. For this size of
    system and the coarse aims of the course project this is acceptable.
    If needed, more formal steady-state detection could be added later.

    Args:
        params:
            RasCSCParameters with the current model parameters.
        u_tgfb:
            External Ras / TGFb drive for this condition.
        t_max:
            Final time for integration. Must be positive. Larger values allow
            more time for the system to settle but increase compute cost.
        y0:
            Optional initial state vector. If None, uses a small positive
            default [C, A, T, M] = [0.01, 0.01, 0.01, 0.01].

    Returns:
        A numpy array with the final state [C, A, T, M] at t = t_max.

    Raises:
        RuntimeError:
            If SciPy is not available or if the ODE solver fails.
        ValueError:
            If the initial state has an incorrect shape.
    """
    # Ensure SciPy is available before attempting integration
    if not SCIPY_AVAILABLE:
        # Raise a runtime error if SciPy is missing
        raise RuntimeError(
            "SciPy is required for ODE integration but is not installed. "
            "Install it with 'pip install scipy' and re-run this script."
        )

    # If no initial state is provided, define a small positive default
    if y0 is None:
        # Define a small positive starting state for all variables
        y0 = np.array([0.01, 0.01, 0.01, 0.01], dtype=float)

    # Check that the initial state has the correct length
    if y0.shape != (4,):
        # Raise an error if the shape is inconsistent
        raise ValueError(
            f"Initial state y0 must be length 4 [C, A, T, M], got shape {y0.shape}."
        )

    # Define the time span for integration as (0, t_max)
    t_span = (0.0, float(t_max))

    # Construct a time evaluation grid for inspecting the trajectory
    t_eval = np.linspace(t_span[0], t_span[1], num=200)

    # Call SciPy's solve_ivp to integrate the ODE system
    sol = solve_ivp(
        fun=lambda t, y: ras_csc_ode(t, y, params=params, u_tgfb=u_tgfb),
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        vectorized=False,
        rtol=1e-6,
        atol=1e-9,
    )

    # Check if the solver reported success
    if not sol.success:
        # Raise a runtime error with the solver's message
        raise RuntimeError(
            f"ODE solver failed with message: {sol.message}"
        )

    # Extract the final state from the last time point
    final_state = sol.y[:, -1]

    # Return the final state as a numpy array
    return final_state


def state_to_predicted_targets(state: np.ndarray) -> Tuple[float, float, float, float]:
    """
    Map the 4-dimensional model state [C, A, T, M] to predicted targets.

    The current implementation uses an identity mapping:

        C_pred = C
        A_pred = A
        T_pred = T
        M_pred = M

    If needed, this can later be extended with affine transformations
    to better align model output with the z-score scale of the RNA-seq
    modules.

    Args:
        state:
            Numpy array of length 4 with [C, A, T, M].

    Returns:
        Tuple (C_pred, A_pred, T_pred, M_pred).

    Raises:
        ValueError:
            If the state vector does not have length 4.
    """
    # Validate the length of the state vector
    if state.shape != (4,):
        # Raise an error if the length is incorrect
        raise ValueError(
            f"State vector must be length 4 [C, A, T, M], got shape {state.shape}."
        )

    # Unpack the state into named variables
    C, A, T, M = map(float, state)

    # Return the predicted targets as a tuple
    return C, A, T, M


def compute_residuals(
    theta: np.ndarray,
    targets: List[ConditionTarget],
    u_map: Dict[Tuple[str, str], float],
) -> np.ndarray:
    """
    Compute residuals between model predictions and calibration targets.

    The parameter vector theta is mapped into a RasCSCParameters instance.
    For each condition, the model is simulated to steady state and the
    predicted [C, A, T, M] values are compared with the corresponding
    [C_target, A_target, T_target, M_target]. Residuals across all
    conditions and all four modules are concatenated into a single array.

    This function is used both by least_squares during fitting and by
    post-fit diagnostics.

    Args:
        theta:
            Numpy array of parameter values. The order is:

                0: k_C_base
                1: k_C_M
                2: K_C
                3: n_C
                4: k_C_decay
                5: k_A_C
                6: k_A_decay
                7: k_T_A
                8: k_T_input
                9: k_T_decay
               10: k_M_T
               11: k_M_decay

        targets:
            List of ConditionTarget objects giving desired module values.
        u_map:
            Mapping from (dataset, condition) to u_TGFb input.

    Returns:
        A numpy array of residuals stacked across conditions and modules.

    Raises:
        KeyError:
            If any dataset-condition pair is missing a u_TGFb entry.
        RuntimeError:
            If simulation fails for any condition.
    """
    # Create a RasCSCParameters object from the theta vector
    params = RasCSCParameters(
        k_C_base=float(theta[0]),
        k_C_M=float(theta[1]),
        K_C=float(theta[2]),
        n_C=float(theta[3]),
        k_C_decay=float(theta[4]),
        k_A_C=float(theta[5]),
        k_A_decay=float(theta[6]),
        k_T_A=float(theta[7]),
        k_T_input=float(theta[8]),
        k_T_decay=float(theta[9]),
        k_M_T=float(theta[10]),
        k_M_decay=float(theta[11]),
    )

    # Initialize a list to collect residuals
    residuals: List[float] = []

    # Loop over each condition target
    for ct in targets:
        # Build the key for this dataset-condition pair
        key = (ct.dataset, ct.condition)

        # Check that a u_TGFb value exists for this key
        if key not in u_map:
            # Raise a KeyError if the mapping is incomplete
            raise KeyError(
                f"No u_TGFb value found for (dataset='{ct.dataset}', "
                f"condition='{ct.condition}')."
            )

        # Retrieve the Ras / TGFb drive for this condition
        u_tgfb = u_map[key]

        # Simulate the system to steady state for this condition
        try:
            final_state = simulate_to_steady_state(
                params=params,
                u_tgfb=u_tgfb,
            )
        except Exception as exc:
            # If simulation fails, raise an error that includes context
            raise RuntimeError(
                f"Simulation failed for dataset='{ct.dataset}', "
                f"condition='{ct.condition}'."
            ) from exc

        # Map the final state to predicted module values
        C_pred, A_pred, T_pred, M_pred = state_to_predicted_targets(
            final_state
        )

        # Append residuals (predicted - target) for each module
        residuals.append(C_pred - ct.C_target)
        residuals.append(A_pred - ct.A_target)
        residuals.append(T_pred - ct.T_target)
        residuals.append(M_pred - ct.M_target)

    # Convert the list of residuals to a numpy array
    residual_array = np.array(residuals, dtype=float)

    # Return the residual array
    return residual_array


def compute_aic_bic_from_result(result: "OptimizeResult") -> Tuple[float, float]:
    """
    Compute AIC and BIC for a nonlinear least-squares fit.

    This function assumes that the model was fitted by minimizing the sum
    of squared residuals under a Gaussian error model. It uses the
    residual vector and the number of free parameters to derive an
    approximate log-likelihood, and from that computes the Akaike
    Information Criterion (AIC) and Bayesian Information Criterion (BIC).

    These criteria are intended for model comparison on the same dataset:
    a lower AIC or BIC indicates a better trade off between fit quality
    and model complexity.

    Args:
        result:
            SciPy OptimizeResult returned by `scipy.optimize.least_squares`.
            It must contain fields:
              - `fun` : the residual vector at the solution.
              - `x`   : the parameter vector at the solution.

    Returns:
        A tuple `(aic, bic)`:
          - `aic` : Akaike Information Criterion.
          - `bic` : Bayesian Information Criterion.

    Raises:
        ValueError:
            If the residual vector is empty, if there are no parameters,
            or if the residual sum of squares is non-positive.
    """
    # Convert the residual vector to a numpy array
    residuals = np.asarray(result.fun, dtype=float)

    # Convert the parameter vector to a numpy array
    params = np.asarray(result.x, dtype=float)

    # Determine the number of data points from the residual vector length
    n_data = residuals.size

    # Determine the number of fitted parameters from the parameter length
    n_params = params.size

    # Guard against missing or degenerate data
    if n_data <= 0:
        # Raise an error if there are no residuals
        raise ValueError(
            "compute_aic_bic_from_result: residual vector is empty; cannot compute AIC/BIC.")

    # Guard against missing or degenerate parameters
    if n_params <= 0:
        # Raise an error if there are no parameters
        raise ValueError(
            "compute_aic_bic_from_result: parameter vector is empty; cannot compute AIC/BIC.")

    # Compute the residual sum of squares RSS = sum_i r_i^2
    rss = float(np.sum(residuals ** 2))

    # Guard against non-positive RSS which would break the log step
    if rss <= 0.0:
        # Raise an error if RSS is not positive
        raise ValueError(
            f"compute_aic_bic_from_result: residual sum of squares is non-positive (RSS={rss}); "
            "check residuals and upstream code."
        )

    # Compute the unbiased estimate of the error variance sigma^2 = RSS / n_data
    sigma2_hat = rss / float(n_data)

    # Compute the log-likelihood up to an additive constant for Gaussian errors
    log_likelihood = -0.5 * float(n_data) * (
        np.log(2.0 * np.pi * sigma2_hat) + 1.0
    )

    # Compute Akaike Information Criterion: AIC = 2k - 2 logL
    aic = 2.0 * float(n_params) - 2.0 * log_likelihood

    # Compute Bayesian Information Criterion: BIC = k log(n) - 2 logL
    bic = float(n_params) * np.log(float(n_data)) - 2.0 * log_likelihood

    # Return the pair of information criteria
    return aic, bic


def compute_param_uncertainty_from_result(
    result: "OptimizeResult",
    param_names: List[str],
    z_value: float = 1.96,
) -> Dict[str, Dict[str, float]]:
    """
    Compute approximate standard errors and confidence intervals for
    fitted parameters using the Jacobian from a least-squares fit.

    This function uses the standard Gauss Newton approximation for
    nonlinear least squares: the parameter covariance matrix is estimated
    as sigma^2 * (J^T J)^(-1), where J is the residual Jacobian at the
    optimum and sigma^2 is the estimated residual variance. The diagonal
    of this covariance matrix provides the squared standard errors for
    each parameter, which can then be used to derive normal-based
    confidence intervals.

    These intervals are local approximations around the fitted point;
    they are not guaranteed to capture global nonlinearities but they
    are sufficient to show which parameters are tightly constrained and
    which are essentially free given the data.

    Args:
        result:
            SciPy OptimizeResult returned by `scipy.optimize.least_squares`,
            containing the residuals, parameters, and Jacobian.
        param_names:
            A list of parameter names in the same order as `result.x`.
        z_value:
            The z-score for the desired confidence interval. For a
            95 percent CI under a normal approximation, use z = 1.96.

    Returns:
        A dictionary mapping each parameter name to a dictionary with keys:
          - "estimate": point estimate (float)
          - "se":      standard error (float)
          - "ci_lower": lower bound of CI (float)
          - "ci_upper": upper bound of CI (float)

    Raises:
        ValueError:
            If the dimensions of the Jacobian are inconsistent, if there
            are too few data points relative to parameters (N <= k),
            or if the parameter name list does not match the vector length.
    """
    # Convert residual vector to numpy array
    residuals = np.asarray(result.fun, dtype=float)

    # Convert parameter vector to numpy array
    params = np.asarray(result.x, dtype=float)

    # Convert Jacobian to numpy array
    J = np.asarray(result.jac, dtype=float)

    # Determine number of data points and parameters from Jacobian shape
    n_data, n_params = J.shape

    # Guard against mismatched parameter size
    if params.size != n_params:
        # Raise an error if vector and Jacobian disagree
        raise ValueError(
            f"compute_param_uncertainty_from_result: parameter vector length {params.size} "
            f"does not match Jacobian column count {n_params}."
        )

    # Guard against mismatched parameter names length
    if len(param_names) != n_params:
        # Raise an error if the name list length does not match
        raise ValueError(
            f"compute_param_uncertainty_from_result: param_names length {len(param_names)} "
            f"does not match number of parameters {n_params}."
        )

    # Guard against insufficient data points for covariance estimation
    if n_data <= n_params:
        # Raise an error if there are not enough data points
        raise ValueError(
            f"compute_param_uncertainty_from_result: not enough data points (N={n_data}) "
            f"relative to parameters (k={n_params}); cannot estimate covariance reliably."
        )

    # Compute residual sum of squares
    rss = float(np.sum(residuals ** 2))

    # Estimate residual variance using unbiased estimator sigma^2 = RSS / (N - k)
    sigma2_hat = rss / float(n_data - n_params)

    # Compute J^T J for the Gauss Newton approximation
    JTJ = J.T @ J

    # Attempt to invert J^T J directly
    try:
        # Compute the inverse of J^T J
        JTJ_inv = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        # Fall back to pseudo inverse if the matrix is singular or ill conditioned
        JTJ_inv = np.linalg.pinv(JTJ)

    # Estimate covariance matrix as sigma^2 * (J^T J)^(-1)
    cov_matrix = sigma2_hat * JTJ_inv

    # Extract standard errors as square roots of diagonal entries
    se_params = np.sqrt(np.clip(np.diag(cov_matrix), a_min=0.0, a_max=None))

    # Prepare output dictionary
    summary: Dict[str, Dict[str, float]] = {}

    # Iterate over each parameter index and name
    for idx, name in enumerate(param_names):
        # Retrieve the estimate for this parameter
        theta_hat = float(params[idx])

        # Retrieve the standard error for this parameter
        se = float(se_params[idx])

        # Compute the lower confidence bound
        ci_lower = theta_hat - float(z_value) * se

        # Compute the upper confidence bound
        ci_upper = theta_hat + float(z_value) * se

        # Store results in the nested dictionary
        summary[name] = {
            "estimate": theta_hat,
            "se": se,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
        }

    # Return the full parameter uncertainty summary
    return summary


def summarize_local_sensitivity(
    result: "OptimizeResult",
    param_names: List[str],
) -> Dict[str, float]:
    """
    Summarize local parameter sensitivity using the residual Jacobian.

    This function computes a simple root mean square measure of
    how strongly each parameter influences the residuals at the best fit
    point. It uses the Jacobian matrix returned by SciPy's least_squares
    function, where each entry J_ij approximates the derivative of the
    i-th residual with respect to the j-th parameter.

    For each parameter j, the sensitivity score is:

        S_j = sqrt( mean_i( J_ij^2 ) )

    Larger S_j indicates that small changes in parameter j have a larger
    effect on the residuals, suggesting that the parameter is more
    influential and better identified by the data. Very small S_j
    suggests a parameter that is weakly identified.

    Args:
        result:
            SciPy OptimizeResult returned by `scipy.optimize.least_squares`,
            containing the Jacobian at the optimum.
        param_names:
            List of parameter names corresponding to the columns of the
            Jacobian and entries of the parameter vector.

    Returns:
        A dictionary mapping each parameter name to its RMS sensitivity
        score S_j.

    Raises:
        ValueError:
            If the Jacobian dimensions are inconsistent with the number
            of parameters or the provided parameter names.
    """
    # Convert Jacobian to numpy array
    J = np.asarray(result.jac, dtype=float)

    # Determine number of data points and parameters from Jacobian shape
    n_data, n_params = J.shape

    # Guard against mismatch between parameter names and Jacobian columns
    if len(param_names) != n_params:
        # Raise an error if the lengths do not match
        raise ValueError(
            f"summarize_local_sensitivity: param_names length {len(param_names)} "
            f"does not match number of parameters {n_params}."
        )

    # Compute mean of squared Jacobian entries along the data axis
    mean_sq = np.mean(J ** 2, axis=0)

    # Take square root to obtain RMS sensitivity scores
    rms_sensitivity = np.sqrt(mean_sq)

    # Map parameter names to sensitivity scores
    sensitivity_dict: Dict[str, float] = {
        name: float(rms_sensitivity[idx])
        for idx, name in enumerate(param_names)
    }

    # Return the dictionary of sensitivities
    return sensitivity_dict


def fit_parameters(
    initial_params: RasCSCParameters,
    targets: List[ConditionTarget],
    u_map: Dict[Tuple[str, str], float],
) -> Tuple[RasCSCParameters, "OptimizeResult"]:
    """
    Run a least-squares fit to adjust parameters against calibration targets.

    This function uses SciPy's least_squares to modify the RasCSCParameters
    within simple bounds. It is intended as a coarse calibration step that
    finds a parameter set with low residual sum of squares. AIC, BIC,
    approximate confidence intervals, and local sensitivities are computed
    separately based on the returned OptimizeResult.

    The function is careful not to crash the entire pipeline if SciPy
    reports that the maximum number of function evaluations was exceeded:
    in that case it prints a warning and still returns the best parameters
    found so far, along with the full OptimizeResult for inspection.

    Args:
        initial_params:
            Starting RasCSCParameters instance for optimization.
        targets:
            List of ConditionTarget objects.
        u_map:
            Mapping from (dataset, condition) to u_TGFb.

    Returns:
        A tuple (fitted_params, result) where:
          - fitted_params is a new RasCSCParameters instance with fitted values.
          - result is the SciPy OptimizeResult from least_squares.

    Raises:
        RuntimeError:
            If SciPy is not available or if optimization fails in an
            unexpected way that does not yield a result.
    """
    # Ensure SciPy is available for optimization
    if not SCIPY_AVAILABLE:
        # Raise a runtime error if SciPy is missing
        raise RuntimeError(
            "SciPy is required for parameter fitting but is not installed. "
            "Install it with 'pip install scipy' or run without --fit."
        )

    # Pack the initial parameters into a theta vector in the canonical order
    theta0 = np.array(
        [
            initial_params.k_C_base,
            initial_params.k_C_M,
            initial_params.K_C,
            initial_params.n_C,
            initial_params.k_C_decay,
            initial_params.k_A_C,
            initial_params.k_A_decay,
            initial_params.k_T_A,
            initial_params.k_T_input,
            initial_params.k_T_decay,
            initial_params.k_M_T,
            initial_params.k_M_decay,
        ],
        dtype=float,
    )

    # Define simple lower and upper bounds to keep parameters positive
    lower_bounds = np.full_like(theta0, 1e-3)
    upper_bounds = np.full_like(theta0, 10.0)

    # Try to run least_squares with a reasonably generous evaluation limit
    try:
        # Call least_squares with the residual function and bounds
        result = least_squares(
            fun=lambda th: compute_residuals(th, targets=targets, u_map=u_map),
            x0=theta0,
            bounds=(lower_bounds, upper_bounds),
            method="trf",
            max_nfev=2000,
            verbose=1,
        )
    except Exception as exc:
        # Raise a runtime error if optimization fails unexpectedly
        raise RuntimeError(
            f"Least-squares optimization raised an exception: {exc}"
        ) from exc

    # Check whether optimization reported success
    if not result.success:
        # Print a warning but still proceed with the best parameters found
        print(
            f"[WARN] least_squares did not reach SciPy 'success' state: {result.message}"
        )
    else:
        # Print informative message when convergence was achieved
        print(f"[INFO] least_squares converged: {result.message}")

    # Extract the optimized theta vector
    theta_opt = np.asarray(result.x, dtype=float)

    # Construct a RasCSCParameters with fitted values
    fitted_params = RasCSCParameters(
        k_C_base=float(theta_opt[0]),
        k_C_M=float(theta_opt[1]),
        K_C=float(theta_opt[2]),
        n_C=float(theta_opt[3]),
        k_C_decay=float(theta_opt[4]),
        k_A_C=float(theta_opt[5]),
        k_A_decay=float(theta_opt[6]),
        k_T_A=float(theta_opt[7]),
        k_T_input=float(theta_opt[8]),
        k_T_decay=float(theta_opt[9]),
        k_M_T=float(theta_opt[10]),
        k_M_decay=float(theta_opt[11]),
    )

    # Return both the fitted parameters and the OptimizeResult
    return fitted_params, result


def main() -> None:
    """
    Entry point for the evaluate_ras_csc_fit script.

    This function:

    1. Parses command-line arguments.
    2. Loads the calibration targets CSV.
    3. Constructs a default u_TGFb mapping.
    4. Either:
       - Simulates with the default parameter set and prints residuals, or
       - Fits parameters (if --fit is requested) and prints residuals.
    5. For fit mode, computes AIC, BIC, approximate parameter confidence
       intervals, and local parameter sensitivities from the OptimizeResult.
    6. Writes everything to a JSON file with keys:

        - 'initial_params': starting parameters
        - 'fitted_params': optimized parameters (or defaults in no-fit mode)
        - 'conditions': list of dicts with dataset, condition, targets,
                        predictions, and residuals
        - 'fit_diagnostics': optional diagnostics (fit mode only) including
                             RSS, AIC, BIC, parameter uncertainty, and
                             local sensitivities
    """
    # Create an argument parser for command-line interaction
    parser = argparse.ArgumentParser(
        description="Evaluate and optionally fit the Ras-CSC model to RNA-seq calibration targets."
    )

    # Add argument for calibration CSV path with a sensible default
    parser.add_argument(
        "--targets",
        type=str,
        default="data/processed/rnaseq/ras_csc_calibration_targets.csv",
        help="Path to ras_csc_calibration_targets.csv produced by the RNA-seq pipeline.",
    )

    # Add a flag that toggles parameter fitting via least-squares
    parser.add_argument(
        "--fit",
        action="store_true",
        help="If provided, perform a least-squares parameter fit using SciPy.",
    )

    # Add argument for output JSON path, with default under data/processed/model_fits
    parser.add_argument(
        "--output-json",
        type=str,
        default="data/processed/model_fits/model_calibration_results.json",
        help="Path to JSON file where calibration results will be written.",
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Convert targets path to a Path object
    targets_path = Path(args.targets)

    # Convert output JSON path to a Path object
    output_json_path = Path(args.output_json)

    # Load the calibration targets from CSV
    targets = load_calibration_targets(targets_path)

    # Build the default u_TGFb mapping for all dataset-condition pairs
    u_map = default_u_tgfb_map(targets)

    # Initialize the starting parameter set
    initial_params = RasCSCParameters()

    # Set a variable to hold the OptimizeResult if fitting is performed
    opt_result: Optional["OptimizeResult"] = None

    # Decide whether to fit parameters or just evaluate the initial set
    if args.fit:
        # Inform the user that fitting has been requested
        print("[INFO] Performing least-squares parameter fit against RNA-seq targets...")

        # Perform parameter fitting and obtain a new RasCSCParameters instance
        fitted_params, opt_result = fit_parameters(
            initial_params,
            targets=targets,
            u_map=u_map,
        )
    else:
        # If not fitting, just reuse the initial parameters as "fitted"
        print("[INFO] Running model with default parameter values (no fitting).")
        fitted_params = initial_params

    # Prepare a list to hold condition-by-condition summaries
    condition_summaries: List[Dict[str, object]] = []

    # Loop over each condition to compute predictions and residuals
    for ct in targets:
        # Look up the corresponding u_TGFb input
        key = (ct.dataset, ct.condition)
        u_tgfb = u_map[key]

        # Simulate to steady state with the chosen parameters
        final_state = simulate_to_steady_state(
            params=fitted_params,
            u_tgfb=u_tgfb,
        )

        # Map the final state to predicted module values
        C_pred, A_pred, T_pred, M_pred = state_to_predicted_targets(
            final_state
        )

        # Compute residuals for each module
        res_C = C_pred - ct.C_target
        res_A = A_pred - ct.A_target
        res_T = T_pred - ct.T_target
        res_M = M_pred - ct.M_target

        # Build a dictionary summary for this condition
        summary = {
            "dataset": ct.dataset,
            "condition": ct.condition,
            "u_TGFb": u_tgfb,
            "targets": {
                "C_target": ct.C_target,
                "A_target": ct.A_target,
                "T_target": ct.T_target,
                "M_target": ct.M_target,
            },
            "predicted": {
                "C_pred": C_pred,
                "A_pred": A_pred,
                "T_pred": T_pred,
                "M_pred": M_pred,
            },
            "residuals": {
                "C_residual": res_C,
                "A_residual": res_A,
                "T_residual": res_T,
                "M_residual": res_M,
            },
        }

        # Append the summary to the list
        condition_summaries.append(summary)

    # Print a brief summary of condition wise residuals to the console
    print("\n[SUMMARY] Condition-wise residuals (chosen parameters):")
    for s in condition_summaries:
        # Extract a compact string for residuals
        res_str = (
            f"C={s['residuals']['C_residual']:.3f}, "
            f"A={s['residuals']['A_residual']:.3f}, "
            f"T={s['residuals']['T_residual']:.3f}, "
            f"M={s['residuals']['M_residual']:.3f}"
        )
        # Print dataset, condition, and residuals
        print(f"  {s['dataset']:>5s} / {s['condition']:<10s}: {res_str}")

    # Prepare fit diagnostics dictionary for JSON output
    fit_diagnostics: Dict[str, object] = {}

    # Compute information criteria and uncertainty only if we have an OptimizeResult
    if opt_result is not None:
        try:
            # Compute AIC and BIC from the optimization result
            aic, bic = compute_aic_bic_from_result(opt_result)

            # Compute residual sum of squares
            residual_vec = np.asarray(opt_result.fun, dtype=float)
            rss = float(np.sum(residual_vec ** 2))

            # Determine counts of data points and parameters
            n_data = residual_vec.size
            n_params = opt_result.x.size

            # Compute parameter uncertainty from the Jacobian
            param_uncertainty = compute_param_uncertainty_from_result(
                opt_result,
                param_names=PARAM_NAMES,
            )

            # Compute local parameter sensitivities
            local_sensitivity = summarize_local_sensitivity(
                opt_result,
                param_names=PARAM_NAMES,
            )

            # Populate the diagnostics dictionary
            fit_diagnostics = {
                "n_data": int(n_data),
                "n_params": int(n_params),
                "rss": rss,
                "aic": float(aic),
                "bic": float(bic),
                "fit_success": bool(opt_result.success),
                "fit_message": str(opt_result.message),
                "n_function_evals": int(opt_result.nfev),
                "param_uncertainty": param_uncertainty,
                "local_sensitivity": local_sensitivity,
            }

            # Print a compact summary to the console
            print("\n[STATS] Fit diagnostics:")
            print(f"  N_data = {n_data}, N_params = {n_params}")
            print(f"  RSS = {rss:.3f}")
            print(f"  AIC = {aic:.3f}, BIC = {bic:.3f}")

            # Print local sensitivities sorted from largest to smallest
            print("\n[STATS] Local parameter sensitivities (RMS on residuals):")
            for name, value in sorted(
                local_sensitivity.items(),
                key=lambda kv: kv[1],
                reverse=True,
            ):
                print(f"  {name}: {value:.4e}")

            # Print parameter estimates with approximate 95 percent CIs
            print("\n[STATS] Parameter estimates with approximate 95 percent CIs:")
            for name in PARAM_NAMES:
                if name not in param_uncertainty:
                    continue
                info = param_uncertainty[name]
                est = info["estimate"]
                se = info["se"]
                lo = info["ci_lower"]
                hi = info["ci_upper"]
                print(
                    f"  {name}: est={est:.4f}, SE={se:.4f}, "
                    f"CI95=({lo:.4f}, {hi:.4f})"
                )

        except Exception as exc:
            # Print a warning but do not prevent JSON writing
            print(f"[WARN] Could not compute full fit diagnostics: {exc}")

    # Build the directory for the JSON output if it does not exist
    if not output_json_path.parent.exists():
        # Create parent directories as needed
        output_json_path.parent.mkdir(parents=True, exist_ok=True)

    # Prepare the JSON payload with initial and fitted parameters and summaries
    json_payload: Dict[str, object] = {
        "initial_params": asdict(initial_params),
        "fitted_params": asdict(fitted_params),
        "conditions": condition_summaries,
    }

    # Include fit diagnostics only if available
    if fit_diagnostics:
        json_payload["fit_diagnostics"] = fit_diagnostics

    # Write the JSON payload to the specified file
    with output_json_path.open("w", encoding="utf-8") as fh:
        json.dump(json_payload, fh, indent=2)

    # Inform the user where the results were written
    print(f"\n[INFO] Calibration results written to: {output_json_path}")


# Standard guard for script execution
if __name__ == "__main__":
    # Call the main function when the script is executed directly
    main()
