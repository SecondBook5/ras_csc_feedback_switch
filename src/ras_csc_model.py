#!/usr/bin/env python3
"""
Ras–CSC–microenvironment feedback ODE model (calibrated version).

This module implements a small ordinary differential equation (ODE) model
for Ras-driven cancer stem cell (CSC) crosstalk with the microenvironment,
inspired by Yuan et al. (Nature 2022, "Ras drives malignancy through
stem cell crosstalk with the microenvironment").

The goal is not to reproduce every molecular detail, but to capture a minimal
positive feedback loop:

    C (CSC program) -> V (angiogenesis / vasculature)
        -> L (local leptin) -> (L, R) -> M (PI3K–AKT–mTOR)
        -> C (CSC program)

Two model topologies are compared:

- Model A (feedback model):
    mTOR activity M feeds back into CSC expansion through a Hill term.
    This can create a malignant "high C, high V, high L, high M" state
    in addition to a benign or low state.

- Model B (no-feedback model):
    The mTOR-to-CSC feedback is set to zero (k_C_M = 0).
    The same microenvironmental structure is present, but no strong
    positive feedback closes the loop on C.

Both models are encoded through the same ODE right-hand side; they differ
only in parameter values (k_C_M > 0 vs k_C_M = 0).

All variables here are *dimensionless* and interpreted as **normalised
module activities**:

    - C : normalised CSC program (0 ~ Pap, 1 ~ SCC)
    - V : normalised angiogenesis/vascular support
    - L : normalised local leptin
    - R : normalised LEPR competence
    - M : normalised mTOR activity

The C-equation is deliberately aligned with the project plan:

    dC/dt = k_C_base + k_C_M * H_M(M) - k_C_decay * C

where H_M(M) is a Hill function of M with half-saturation K_C and Hill
coefficient n_C. In the original plan, Ras * k_C_baseline appears;
here, that product is absorbed into k_C_base after calibration.

This module provides:

    - RasCSCParams dataclass: container for kinetic parameters.
    - ras_csc_rhs(): ODE right-hand side using the calibrated parameters.
    - simulate_trajectory(): simple fixed-step RK4 integrator.
    - load_calibrated_params(): load Model A/B parameters from JSON.
    - A demo (__main__) that shows:
          * Model A: Pap IC -> low state, SCC IC -> high state
          * Model B: Pap/SCC ICs both collapse to the same low state

The JSON file with calibrated parameters must be created by running:

    python src/calibrate_params.py

prior to using this module.
"""

from __future__ import annotations

# Import standard libraries
import json
import os
import sys

# Import dataclass for structured parameter storage
from dataclasses import dataclass

# Import typing utilities for type hints
from typing import Callable, Tuple

# Import numpy for numerical computations
import numpy as np


@dataclass
class RasCSCParams:
    """
    Container for Ras–CSC feedback model parameters.

    This dataclass groups together all kinetic parameters for the ODE system.
    It is used to cleanly pass parameters into the ODE right-hand side and to
    represent "Model A" (feedback) and "Model B" (no feedback) as simple
    changes in parameter values rather than separate code paths.

    Parameter semantics:

    - CSC dynamics:
        k_C_base : Basal net production term for C. In the original plan,
                   Ras * k_C_baseline appears; here that product is absorbed
                   into this constant.
        k_C_M    : Strength of positive feedback from mTOR M onto CSC
                   production via a Hill function H_M(M).
        K_C      : Half-saturation constant for H_M(M).
        n_C      : Hill coefficient for H_M(M).
        k_C_decay: Linear loss rate for C (turnover/differentiation).

        The C equation is:

            dC/dt = k_C_base + k_C_M * H_M(M) - k_C_decay * C

        with:

            H_M(M) = M^n_C / (K_C^n_C + M^n_C).

        This is the linear form corresponding to the project plan. We do NOT
        use the previous logistic placeholder C * (1 - C); that would conflict
        with the Pap/SCC calibration.

    - Angiogenesis / vasculature:
        k_V_C    : Rate at which CSC activity promotes vascular support V.
        k_V_deg  : Decay rate for V.

        Equation:

            dV/dt = k_V_C * C - k_V_deg * V

    - Leptin:
        k_L0     : Baseline local leptin source term (systemic leak-in).
        k_L_V    : Contribution of angiogenesis V to local leptin.
        k_L_deg  : Decay or clearance rate of leptin.

        Equation:

            dL/dt = k_L0 + k_L_V * V - k_L_deg * L

    - LEPR:
        k_R_TGFb : Rate at which TGFβ (external) induces LEPR competence R.
        u_TGFb   : Dimensionless TGFβ input (external drive).
        k_R_deg  : Decay/turnover rate for LEPR R.

        Equation:

            dR/dt = k_R_TGFb * u_TGFb - k_R_deg * R

    - mTOR:
        k_M_act  : Maximal activation rate of mTOR by leptin–LEPR signaling.
        K_M      : Half-saturation for the Hill dependence on S = L * R.
        m_M      : Hill coefficient in the mTOR equation.
        k_M_deg  : Decay/deactivation rate for mTOR activity.

        Equation:

            S = L * R
            H_LR(S) = S^m_M / (K_M^m_M + S^m_M)
            dM/dt = k_M_act * H_LR(S) - k_M_deg * M

    - Numerical safety:
        clip_state: If True, negative state components are clipped to zero
                    after each integration step.

    In the calibrated workflow:

        - model_calibration_results.json sets these parameters for Model A
          (feedback) and Model B (no feedback).
        - Model B is identical to Model A except that k_C_M is forced to 0
          and k_C_base may be different (as produced by calibration).
    """

    # Basal CSC production term (Ras * k_C_baseline absorbed here)
    k_C_base: float

    # Strength of mTOR-to-CSC positive feedback
    k_C_M: float

    # Half-saturation constant for M in H_M(M)
    K_C: float

    # Hill coefficient for M in H_M(M)
    n_C: float

    # Linear CSC loss rate
    k_C_decay: float

    # CSC-driven angiogenesis gain
    k_V_C: float

    # Angiogenesis decay rate
    k_V_deg: float

    # Baseline leptin influx
    k_L0: float

    # Angiogenesis-driven leptin gain
    k_L_V: float

    # Leptin decay rate
    k_L_deg: float

    # TGFβ-driven LEPR induction rate
    k_R_TGFb: float

    # External TGFβ input level (dimensionless)
    u_TGFb: float

    # LEPR decay rate
    k_R_deg: float

    # Maximal mTOR activation rate by L*R
    k_M_act: float

    # Half-saturation for L*R in the mTOR Hill term
    K_M: float

    # Hill coefficient in the mTOR Hill term
    m_M: float

    # mTOR decay rate
    k_M_deg: float

    # Flag controlling whether state variables are clipped to be non-negative
    clip_state: bool = True


def get_repo_root() -> str:
    """
    Determine the repository root directory based on this script's location.

    The script is expected to live in a 'src/' subdirectory of the project.
    By taking the parent directory of the directory containing this file, we
    obtain the project root in a way that does not depend on the current
    working directory.

    Returns:
        str:
            Absolute path to the repository root directory.

    Raises:
        RuntimeError:
            If the computed repository root does not exist.
    """
    # Compute the absolute path of this file
    script_path: str = os.path.abspath(__file__)
    # Get the directory containing this file (expected to be src/)
    src_dir: str = os.path.dirname(script_path)
    # Move one level up to reach the presumed repo root
    repo_root: str = os.path.dirname(src_dir)

    # Check that the repo root actually exists on disk
    if not os.path.isdir(repo_root):
        # Raise a clear error if the directory is missing
        raise RuntimeError(
            f"[ERROR] Computed repo root does not exist: {repo_root}. "
            f"Check that ras_csc_model.py is located under 'src/'."
        )

    # Return the verified repo root
    return repo_root


def load_calibrated_params() -> Tuple[RasCSCParams, RasCSCParams]:
    """
    Load calibrated parameter sets for Model A and Model B from JSON.

    This function reads the file:

        data/processed/model_calibration_results.json

    relative to the repository root (as determined by get_repo_root()).
    The JSON is expected to contain keys:

        - "model_A_params": dict of parameters for the feedback model
        - "model_B_params": dict of parameters for the no-feedback model

    Each dict must provide numeric values for all fields of RasCSCParams.
    The "clip_state" field may be stored as 0/1 or True/False; it is coerced
    to bool.

    Returns:
        Tuple[RasCSCParams, RasCSCParams]:
            (params_A, params_B) for Model A (feedback) and Model B (no feedback).

    Raises:
        FileNotFoundError:
            If the JSON file does not exist.
        RuntimeError:
            If the JSON is malformed or missing required keys.
    """
    # Determine the repository root directory
    repo_root: str = get_repo_root()

    # Build the full path to the calibration JSON
    json_path: str = os.path.join(
        repo_root,
        "data",
        "processed",
        "model_calibration_results.json",
    )

    # Check that the JSON file exists
    if not os.path.exists(json_path):
        # Raise a clear error instructing the user to run calibrate_params.py
        raise FileNotFoundError(
            f"[ERROR] Calibration file not found at: {json_path}\n"
            f"Run 'python src/calibrate_params.py' first to create it."
        )

    # Attempt to read and parse the JSON file
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            payload = json.load(f)
    except Exception as exc:
        # Wrap and re-raise any parsing error
        raise RuntimeError(
            f"[ERROR] Failed to parse calibration JSON at {json_path}: {exc}"
        ) from exc

    # Extract parameter dictionaries for Model A and Model B
    try:
        raw_A = payload["model_A_params"]
        raw_B = payload["model_B_params"]
    except KeyError as exc:
        # Raise an informative error if required keys are missing
        raise RuntimeError(
            f"[ERROR] Calibration JSON is missing required key {exc!r}. "
            f"Expected keys 'model_A_params' and 'model_B_params'."
        ) from exc

    # Define a helper to coerce numeric dicts into RasCSCParams instances
    def _dict_to_params(raw: dict, label: str) -> RasCSCParams:
        """
        Convert a raw parameter dictionary from JSON into a RasCSCParams object.

        Args:
            raw (dict):
                Dictionary with fields corresponding to RasCSCParams attributes.
            label (str):
                Label used only for error messages (e.g., 'Model A').

        Returns:
            RasCSCParams:
                Dataclass instance constructed from the raw values.

        Raises:
            RuntimeError:
                If required keys are missing.
        """
        # Define the expected keys for RasCSCParams (excluding clip_state)
        required_keys = [
            "k_C_base",
            "k_C_M",
            "K_C",
            "n_C",
            "k_C_decay",
            "k_V_C",
            "k_V_deg",
            "k_L0",
            "k_L_V",
            "k_L_deg",
            "k_R_TGFb",
            "u_TGFb",
            "k_R_deg",
            "k_M_act",
            "K_M",
            "m_M",
            "k_M_deg",
        ]

        # Verify that all required keys are present
        missing = [k for k in required_keys if k not in raw]
        if missing:
            # Raise an error if any keys are missing
            raise RuntimeError(
                f"[ERROR] Calibration parameters for {label} are missing keys: "
                f"{', '.join(missing)}"
            )

        # Extract clip_state, allowing for absence (default to True)
        clip_raw = raw.get("clip_state", True)
        # Coerce clip_state to a proper bool from 0/1 or True/False
        if isinstance(clip_raw, (int, float)):
            clip_state = bool(int(clip_raw))
        else:
            clip_state = bool(clip_raw)

        # Construct the RasCSCParams instance with float casting for safety
        return RasCSCParams(
            k_C_base=float(raw["k_C_base"]),
            k_C_M=float(raw["k_C_M"]),
            K_C=float(raw["K_C"]),
            n_C=float(raw["n_C"]),
            k_C_decay=float(raw["k_C_decay"]),
            k_V_C=float(raw["k_V_C"]),
            k_V_deg=float(raw["k_V_deg"]),
            k_L0=float(raw["k_L0"]),
            k_L_V=float(raw["k_L_V"]),
            k_L_deg=float(raw["k_L_deg"]),
            k_R_TGFb=float(raw["k_R_TGFb"]),
            u_TGFb=float(raw["u_TGFb"]),
            k_R_deg=float(raw["k_R_deg"]),
            k_M_act=float(raw["k_M_act"]),
            K_M=float(raw["K_M"]),
            m_M=float(raw["m_M"]),
            k_M_deg=float(raw["k_M_deg"]),
            clip_state=clip_state,
        )

    # Convert raw dicts into parameter dataclasses
    params_A: RasCSCParams = _dict_to_params(raw_A, "Model A")
    params_B: RasCSCParams = _dict_to_params(raw_B, "Model B")

    # Return both parameter sets
    return params_A, params_B


def ras_csc_rhs(
    t: float,
    y: np.ndarray,
    params: RasCSCParams,
) -> np.ndarray:
    """
    Evaluate the Ras–CSC ODE right-hand side at time t and state y.

    The state vector y has length 5 with entries:

        y[0] = C : CSC program / CSC abundance (normalised)
        y[1] = V : angiogenic / vascular support (normalised)
        y[2] = L : local leptin level (normalised)
        y[3] = R : LEPR competence (normalised)
        y[4] = M : mTOR activity (normalised)

    Dynamics:

        1) CSC program C:

            H_M(M) = M^n_C / (K_C^n_C + M^n_C)

            dC/dt = k_C_base + k_C_M * H_M(M) - k_C_decay * C

           Here k_C_base + k_C_M * H_M(M) is the effective production rate
           (including Ras and feedback), and k_C_decay * C is a linear loss.
           Pap and SCC normalised levels (C ~ 0 and C ~ 1) are achievable
           steady states for suitable parameter values.

        2) Angiogenesis V:

            dV/dt = k_V_C * C - k_V_deg * V

           CSC activity promotes angiogenesis; V decays when not driven.

        3) Leptin L:

            dL/dt = k_L0 + k_L_V * V - k_L_deg * L

           Vascular support increases leptin delivery; there is a baseline
           source k_L0 and clearance k_L_deg.

        4) LEPR competence R:

            dR/dt = k_R_TGFb * u_TGFb - k_R_deg * R

           TGFβ (u_TGFb) induces LEPR; R decays with rate k_R_deg.

        5) mTOR activity M:

            S = L * R
            H_LR(S) = S^m_M / (K_M^m_M + S^m_M)

            dM/dt = k_M_act * H_LR(S) - k_M_deg * M

           Leptin–LEPR signaling activates mTOR, which decays otherwise.

    The distinction between Model A and Model B is entirely through k_C_M:

        - Model A: k_C_M > 0 => positive feedback C→V→L→M→C.
        - Model B: k_C_M = 0 => loop is broken at M→C.

    Args:
        t (float):
            Current time. Not used explicitly (system is autonomous),
            but included for compatibility with standard ODE solver APIs.
        y (np.ndarray):
            Current state vector, shape (5,).
        params (RasCSCParams):
            Parameter set controlling all kinetic terms.

    Returns:
        np.ndarray:
            Time derivative dy/dt, shape (5,).

    Raises:
        ValueError:
            If the input state does not have length 5.
    """
    # Ensure the state vector has the correct dimension
    if y.shape[0] != 5:
        raise ValueError(
            f"State vector y must have length 5 (C, V, L, R, M); "
            f"got length {y.shape[0]} instead."
        )

    # Optionally clip negative values for numerical safety
    if params.clip_state:
        y = np.maximum(y, 0.0)

    # Unpack the state variables
    C: float = float(y[0])
    V: float = float(y[1])
    L: float = float(y[2])
    R: float = float(y[3])
    M: float = float(y[4])

    # Compute Hill function H_M(M) for mTOR feedback on C
    M_power: float = M ** params.n_C
    K_C_power: float = params.K_C ** params.n_C
    H_M: float = 0.0
    if M_power + K_C_power > 0.0:
        H_M = M_power / (M_power + K_C_power)

    # Compute CSC derivative using the linear form (NOT logistic)
    dC_dt: float = params.k_C_base + params.k_C_M * H_M - params.k_C_decay * C

    # Compute angiogenesis derivative: C drives V, which decays
    dV_dt: float = params.k_V_C * C - params.k_V_deg * V

    # Compute leptin derivative: baseline + V-driven - clearance
    dL_dt: float = params.k_L0 + params.k_L_V * V - params.k_L_deg * L

    # Compute LEPR competence derivative: TGFβ-driven induction - decay
    dR_dt: float = params.k_R_TGFb * params.u_TGFb - params.k_R_deg * R

    # Compute leptin–LEPR signal S
    S: float = L * R

    # Compute Hill function H_LR(S) for mTOR activation
    S_power: float = S ** params.m_M
    K_M_power: float = params.K_M ** params.m_M
    H_LR: float = 0.0
    if S_power + K_M_power > 0.0:
        H_LR = S_power / (S_power + K_M_power)

    # Compute mTOR derivative: activation by L*R - decay
    dM_dt: float = params.k_M_act * H_LR - params.k_M_deg * M

    # Pack derivatives into a NumPy array
    dydt: np.ndarray = np.array(
        [dC_dt, dV_dt, dL_dt, dR_dt, dM_dt],
        dtype=float,
    )

    # Return the derivative vector
    return dydt


def simulate_trajectory(
    rhs: Callable[[float, np.ndarray, RasCSCParams], np.ndarray],
    y0: np.ndarray,
    params: RasCSCParams,
    t_start: float,
    t_end: float,
    dt: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate the Ras–CSC ODE system using a fixed-step RK4 integrator.

    This function performs an explicit fourth-order Runge–Kutta integration
    for the autonomous ODE system:

        dy/dt = rhs(t, y, params)

    It is intended for exploratory simulations and parameter sweeps, not for
    stiff systems.

    Args:
        rhs (Callable[[float, np.ndarray, RasCSCParams], np.ndarray]):
            Function that computes dy/dt given (t, y, params).
            For this project, this will usually be ras_csc_rhs.
        y0 (np.ndarray):
            Initial state, shape (5,). Must represent [C, V, L, R, M].
        params (RasCSCParams):
            Parameter set controlling the dynamics.
        t_start (float):
            Starting time for the integration.
        t_end (float):
            Final time for the integration; must satisfy t_end > t_start.
        dt (float):
            Time step for the RK4 scheme; must be positive.

    Returns:
        Tuple[np.ndarray, np.ndarray]:
            (times, traj) where:

                - times is a 1D array of shape (N_steps + 1,)
                - traj  is a 2D array of shape (N_steps + 1, 5)

            containing the time points and corresponding state vectors.

    Raises:
        ValueError:
            If dt is non-positive, t_end <= t_start, or y0 has wrong length.
    """
    # Validate the time interval
    if t_end <= t_start:
        raise ValueError(
            f"t_end must be greater than t_start; "
            f"got t_start={t_start}, t_end={t_end}."
        )

    # Validate the step size
    if dt <= 0.0:
        raise ValueError(f"dt must be positive; got dt={dt}.")

    # Convert initial state to NumPy array
    y0_arr = np.asarray(y0, dtype=float)

    # Check that the initial state has the correct dimension
    if y0_arr.shape[0] != 5:
        raise ValueError(
            f"Initial state y0 must have length 5 (C, V, L, R, M); "
            f"got length {y0_arr.shape[0]} instead."
        )

    # Compute the number of integration steps
    n_steps_float: float = (t_end - t_start) / dt
    n_steps: int = int(np.floor(n_steps_float))

    # Ensure at least one step is taken
    if n_steps < 1:
        raise ValueError(
            f"Number of RK4 steps from t_start={t_start} to t_end={t_end} "
            f"with dt={dt} is less than 1. Adjust t_end or dt."
        )

    # Create time grid
    times: np.ndarray = np.linspace(
        t_start,
        t_start + n_steps * dt,
        n_steps + 1,
    )

    # Allocate trajectory array
    traj: np.ndarray = np.zeros((n_steps + 1, 5), dtype=float)

    # Set initial condition
    traj[0, :] = y0_arr

    # Initialize current state and time
    y_curr: np.ndarray = y0_arr.copy()
    t_curr: float = t_start

    # Perform RK4 integration
    for i in range(n_steps):
        # Compute k1 at current state
        k1: np.ndarray = rhs(t_curr, y_curr, params)

        # Compute k2 at midpoint using k1
        k2: np.ndarray = rhs(
            t_curr + 0.5 * dt,
            y_curr + 0.5 * dt * k1,
            params,
        )

        # Compute k3 at midpoint using k2
        k3: np.ndarray = rhs(
            t_curr + 0.5 * dt,
            y_curr + 0.5 * dt * k2,
            params,
        )

        # Compute k4 at next time using k3
        k4: np.ndarray = rhs(
            t_curr + dt,
            y_curr + dt * k3,
            params,
        )

        # Combine k1–k4 to update the state (RK4)
        delta_y: np.ndarray = (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        # Update the state
        y_next: np.ndarray = y_curr + delta_y

        # Enforce non-negativity if clipping is enabled
        if params.clip_state:
            y_next = np.maximum(y_next, 0.0)

        # Advance time
        t_curr = t_curr + dt

        # Store the new state
        traj[i + 1, :] = y_next

        # Prepare for next step
        y_curr = y_next

    # Return the time grid and trajectory
    return times, traj


if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Demonstration: show how Model A vs Model B behave with calibrated
    # parameters when started from Pap-like and SCC-like initial states.
    # ------------------------------------------------------------------

    try:
        # Load calibrated parameters for Model A (feedback) and Model B (no feedback)
        params_A, params_B = load_calibrated_params()
    except Exception as exc:
        # If anything goes wrong, print the error and exit with non-zero status
        print(str(exc), file=sys.stderr)
        sys.exit(1)

    # Define a Papilloma-like initial condition in normalised space
    # (C, V, L, R, M) all low
    y0_pap: np.ndarray = np.array(
        [0.0, 0.0, 0.0, 0.0, 0.0],
        dtype=float,
    )

    # Define an SCC-like initial condition in normalised space
    # (C, V, L, R, M) all high
    y0_scc: np.ndarray = np.array(
        [1.0, 1.0, 1.0, 1.0, 1.0],
        dtype=float,
    )

    # Set integration parameters for the demo
    t_start: float = 0.0
    t_end: float = 200.0
    dt: float = 0.01

    # Simulate Model A from Pap-like initial condition
    times_A_pap, traj_A_pap = simulate_trajectory(
        rhs=ras_csc_rhs,
        y0=y0_pap,
        params=params_A,
        t_start=t_start,
        t_end=t_end,
        dt=dt,
    )

    # Simulate Model A from SCC-like initial condition
    times_A_scc, traj_A_scc = simulate_trajectory(
        rhs=ras_csc_rhs,
        y0=y0_scc,
        params=params_A,
        t_start=t_start,
        t_end=t_end,
        dt=dt,
    )

    # Simulate Model B from Pap-like initial condition
    times_B_pap, traj_B_pap = simulate_trajectory(
        rhs=ras_csc_rhs,
        y0=y0_pap,
        params=params_B,
        t_start=t_start,
        t_end=t_end,
        dt=dt,
    )

    # Simulate Model B from SCC-like initial condition
    times_B_scc, traj_B_scc = simulate_trajectory(
        rhs=ras_csc_rhs,
        y0=y0_scc,
        params=params_B,
        t_start=t_start,
        t_end=t_end,
        dt=dt,
    )

    # Extract final states
    final_A_pap = traj_A_pap[-1, :]
    final_A_scc = traj_A_scc[-1, :]
    final_B_pap = traj_B_pap[-1, :]
    final_B_scc = traj_B_scc[-1, :]

    # Print a concise summary
    print("=== Ras–CSC feedback model demo (calibrated) ===")
    print("Parameters Model A (feedback):", params_A)
    print("Parameters Model B (no feedback):", params_B)
    print("")
    print("Final state Model A from Pap-like initial condition (C, V, L, R, M):")
    print(final_A_pap)
    print("")
    print("Final state Model A from SCC-like initial condition (C, V, L, R, M):")
    print(final_A_scc)
    print("")
    print("Final state Model B from Pap-like initial condition (C, V, L, R, M):")
    print(final_B_pap)
    print("")
    print("Final state Model B from SCC-like initial condition (C, V, L, R, M):")
    print(final_B_scc)
    print("")
    print("Interpretation:")
    print("  - Model A (feedback) supports a low Pap-like state and a high SCC-like state")
    print("    under the same parameter set (initial conditions select the attractor).")
    print("  - Model B (no feedback) collapses Pap and SCC initial conditions")
    print("    to the same low state, consistent with the structural argument.")
