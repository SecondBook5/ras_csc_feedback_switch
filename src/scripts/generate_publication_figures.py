def create_figure1_calibration() -> None:
    """
    Create calibration / validation plots as SEPARATE panels instead of one bulk grid.

    This function replaces the previous 4-panel composite Figure 1 and instead saves
    four individual figures:
      - Figure1A_Model_vs_Data.{png,pdf}
      - Figure1B_Residuals_by_Condition.{png,pdf}
      - Figure1C_Fit_Quality.{png,pdf}
      - Figure1D_Progression_Trends.{png,pdf}

    The underlying data and visual logic are identical to the original implementation,
    but decoupled so each panel can be placed independently in the manuscript.

    Raises:
        FileNotFoundError: If model_vs_data_CORRECTED.csv cannot be found.
        ValueError: If required columns are missing from the calibration CSV.
    """
    # Print status for tracking in the log
    print("\n[FIGURE 1] Creating SEPARATE calibration panels...")

    # Load the calibration data
    csv_path = RESULTS_DIR / "model_vs_data_CORRECTED.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"[ERROR] Missing calibration file: {csv_path}")

    # Read the CSV into a DataFrame
    df = pd.read_csv(csv_path)

    # Define required columns for safety
    required_cols = {
        "condition",
        "module",
        "data_zscore",
        "model_zscore",
        "model_raw",
        "residual",
    }
    missing_cols = required_cols.difference(df.columns)
    if missing_cols:
        raise ValueError(
            f"[ERROR] model_vs_data_CORRECTED.csv missing columns: {missing_cols}. "
            f"Found: {list(df.columns)}"
        )

    # Define module info
    modules = ["C", "A", "T", "M"]
    module_names = {
        "C": "CSC Fraction",
        "A": "Angiogenesis",
        "T": "TGFβ Activity",
        "M": "mTOR Activity",
    }

    # -------------------------------------------------------------------------
    # Panel A: Model vs Data (scatter)  → Figure1A_Model_vs_Data
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 5))

    # Plot each module as a different marker group
    for mod in modules:
        # Filter rows for this module
        df_mod = df[df["module"] == mod]

        # Scatter with edgecolor
        ax.scatter(
            df_mod["data_zscore"],
            df_mod["model_zscore"],
            s=120,
            alpha=0.8,
            label=module_names[mod],
            edgecolor="black",
            linewidth=1.5,
        )

    # Compute limits for identity line
    lims_min = df[["data_zscore", "model_zscore"]].min().min() - 0.2
    lims_max = df[["data_zscore", "model_zscore"]].max().max() + 0.2
    lims = [lims_min, lims_max]

    # Plot identity line
    ax.plot(lims, lims, "k--", alpha=0.5, linewidth=2, label="Identity")

    # Compute overall Pearson correlation
    from scipy.stats import pearsonr

    r_val, p_val = pearsonr(df["data_zscore"], df["model_zscore"])

    # Annotate correlation in the top-left corner
    ax.text(
        0.05,
        0.95,
        f"Overall: r={r_val:.3f}, p={p_val:.2e}",
        transform=ax.transAxes,
        fontsize=11,
        weight="bold",
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    # Label axes and style
    ax.set_xlabel("Data (z-score)", fontsize=13, weight="bold")
    ax.set_ylabel("Model (z-score)", fontsize=13, weight="bold")
    ax.set_title("Figure 1A: Model vs Data", fontsize=15, weight="bold")
    ax.legend(frameon=True, shadow=True, loc="lower right")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # Save panel A
    for fmt in ["png", "pdf"]:
        out_path = FIG_DIR / f"Figure1A_Model_vs_Data.{fmt}"
        plt.savefig(out_path, dpi=600, bbox_inches="tight")
        print(f"  [SAVED] {out_path}")
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Panel B: Residual heatmap  → Figure1B_Residuals_by_Condition
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 5))

    # Build residual matrix: rows = condition, columns = module
    pivot = df.pivot_table(values="residual", index="condition", columns="module")
    # Ensure column order
    pivot = pivot[modules]

    # Create image plot
    im = ax.imshow(
        pivot.values,
        cmap="RdBu_r",
        aspect="auto",
        vmin=-1.5,
        vmax=1.5,
    )

    # Set tick labels
    ax.set_xticks(range(len(modules)))
    ax.set_xticklabels([module_names[m] for m in modules], rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)

    # Add residual values inside cells
    for i in range(len(pivot.index)):
        for j in range(len(modules)):
            val = pivot.values[i, j]
            text_color = "white" if abs(val) > 0.8 else "black"
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontsize=10,
                weight="bold",
            )

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Residual (z-score)", fontsize=11, weight="bold")

    # Label axes
    ax.set_xlabel("Module", fontsize=13, weight="bold")
    ax.set_ylabel("Condition", fontsize=13, weight="bold")
    ax.set_title("Figure 1B: Residuals by Condition", fontsize=15, weight="bold")

    # Save panel B
    for fmt in ["png", "pdf"]:
        out_path = FIG_DIR / f"Figure1B_Residuals_by_Condition.{fmt}"
        plt.savefig(out_path, dpi=600, bbox_inches="tight")
        print(f"  [SAVED] {out_path}")
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Panel C: Fit quality summary (text box)  → Figure1C_Fit_Quality
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.axis("off")

    # Compute fit metrics
    rss = float((df["residual"] ** 2).sum())
    rmse = float(np.sqrt(rss / len(df)))
    n_data = int(len(df))
    n_params = 11  # your chosen number of free parameters

    metrics_text = f"""
    FIT QUALITY METRICS

    RSS (Residual Sum of Squares):    {rss:.3f}
    RMSE (Root Mean Squared Error):   {rmse:.3f} z-score units

    N data points:                    {n_data}
    N free parameters:                {n_params}
    Degrees of freedom:               {n_data - n_params}

    INTERPRETATION:
    • RMSE < 0.7 indicates good fit for this model size
    • Model captures Normal → Papilloma → SCC progression
    • LeprKO phenotype: reduced CSC and mTOR outputs
    • Threshold-like response supports bistability hypothesis
    """

    ax.text(
        0.05,
        0.95,
        metrics_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        family="monospace",
        bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.8, pad=1),
    )
    ax.set_title("Figure 1C: Fit Quality Summary", fontsize=15, weight="bold")

    # Save panel C
    for fmt in ["png", "pdf"]:
        out_path = FIG_DIR / f"Figure1C_Fit_Quality.{fmt}"
        plt.savefig(out_path, dpi=600, bbox_inches="tight")
        print(f"  [SAVED] {out_path}")
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Panel D: Progression trends (raw model outputs)  → Figure1D_Progression_Trends
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))

    # Fixed condition order to line up with biology
    cond_order = ["Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO"]

    for mod in modules:
        # Filter rows for this module
        df_mod = df[df["module"] == mod]

        # Compute mean raw model output per condition in the desired order
        means = []
        for cond in cond_order:
            subset = df_mod[df_mod["condition"] == cond]
            if len(subset) > 0:
                means.append(float(subset["model_raw"].mean()))
            else:
                means.append(np.nan)

        # Plot the trajectory for this module
        ax.plot(
            cond_order,
            means,
            "o-",
            linewidth=2,
            markersize=8,
            label=module_names[mod],
            alpha=0.8,
        )

    # Horizontal threshold line (e.g., 0.5)
    ax.axhline(
        0.5,
        color="red",
        linestyle="--",
        alpha=0.3,
        linewidth=1.5,
    )

    # Label axes and style
    ax.set_xlabel("Condition", fontsize=13, weight="bold")
    ax.set_ylabel("Model Output (raw)", fontsize=13, weight="bold")
    ax.set_title("Figure 1D: Progression Trends", fontsize=15, weight="bold")
    ax.legend(frameon=True, shadow=True)
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.set_xticklabels(cond_order, rotation=45, ha="right")
    ax.set_ylim([0, 1])

    # Save panel D
    for fmt in ["png", "pdf"]:
        out_path = FIG_DIR / f"Figure1D_Progression_Trends.{fmt}"
        plt.savefig(out_path, dpi=600, bbox_inches="tight")
        print(f"  [SAVED] {out_path}")
    plt.close(fig)
