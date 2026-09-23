"""Render the Q1 high-pass sensitivity figures from frozen analysis tables.

Run after q1_highpass_real and q1_filter_simulation have written their CSVs:
    python q1_highpass_figures.py highpass_analysis
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal


AUDIT_SCRIPTS = os.environ.get("NATURE_FIGURE_SCRIPTS")
if AUDIT_SCRIPTS:
    sys.path.insert(0, AUDIT_SCRIPTS)
    from audit_panel_alignment import require_matplotlib_panel_alignment
else:
    require_matplotlib_panel_alignment = None


HP_VALUES = (0.10, 0.20, 0.25, 0.50)
HP_COLORS = {
    0.10: "#285c9a",
    0.20: "#5386b3",
    0.25: "#b7733b",
    0.50: "#a5393b",
}
SOURCE_UNITS = "source units"
RECORDS = ("A", "B")
TASKS = (1, 2)
CHANNELS = ("F3", "Fz", "F4")

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "savefig.facecolor": "white",
    }
)


def finish(
    fig: plt.Figure,
    directory: Path,
    stem: str,
    qa_directory: Path | None,
    exclude_axes=None,
) -> None:
    fig.canvas.draw()
    if qa_directory is not None:
        if require_matplotlib_panel_alignment is None:
            raise RuntimeError("QA requested but NATURE_FIGURE_SCRIPTS is unavailable")
        require_matplotlib_panel_alignment(
            fig,
            json_out=str(directory / f"{stem}.alignment.json"),
            overlay_svg=str(qa_directory / f"{stem}.alignment.svg"),
            tolerance_pt=1.5,
            gutter_tolerance_pt=1.5,
            exclude_axes=exclude_axes or [],
            strict=True,
        )
    fig.savefig(directory / f"{stem}.png", dpi=600, bbox_inches="tight")
    if qa_directory is not None:
        fig.savefig(qa_directory / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def validate_sources(
    real: pd.DataFrame,
    curves: pd.DataFrame,
    sim: pd.DataFrame,
    example: pd.DataFrame,
) -> None:
    """Reject incomplete or inconsistent tables before any figure is drawn."""
    keys = ["record", "task", "channel", "highpass"]
    expected = {
        (record, task, channel, hp)
        for record in RECORDS
        for task in TASKS
        for channel in CHANNELS
        for hp in HP_VALUES
    }
    for name, frame in (("real", real), ("curves", curves)):
        if not set(keys).issubset(frame.columns):
            raise ValueError(f"{name} table is missing key columns")
        observed = set(map(tuple, frame[keys].drop_duplicates().itertuples(index=False, name=None)))
        if observed != expected:
            raise ValueError(f"{name} table does not cover the complete 12 × 4 matrix")
    if len(real) != 48 or real.duplicated(keys).any():
        raise ValueError("Real-data table must have one row per matrix/cutoff cell")
    if real[["delta_mean_250_500", "split_half_r", "baseline_rms"]].isna().any().any():
        raise ValueError("A plotted real-data measurement is missing")
    if len(curves) != 48 * 411 or curves.groupby(keys).size().ne(411).any():
        raise ValueError("Each real-data waveform must contain exactly 411 samples")
    analysis = curves.loc[curves.timeSec.ge(-0.2)]
    if analysis[["left", "right", "delta"]].isna().any().any():
        raise ValueError("An analysis-interval waveform sample is missing")
    if not np.allclose(analysis.delta, analysis.right - analysis.left, atol=1e-8):
        raise ValueError("Right-minus-left waveform disagrees with condition curves")
    means = (
        curves.loc[curves.timeSec.ge(0.25) & curves.timeSec.lt(0.5)]
        .groupby(keys, as_index=False).delta.mean()
    )
    merged = real.merge(means, on=keys, validate="one_to_one")
    if not np.allclose(merged.delta_mean_250_500, merged.delta, atol=1e-8):
        raise ValueError("Real-data window means disagree with the waveform CSV")

    sim_keys = ["record", "task", "channel", "templateHp", "replicate", "highpass", "condition"]
    if sim.duplicated(sim_keys).any() or len(sim) != 4 * 2 * 12 * 4 * 3 * 3:
        raise ValueError("Simulation recovery table has missing/duplicate combinations")
    if set(sim.condition) != {"Left", "Right", "Delta"}:
        raise ValueError("Simulation condition labels are incomplete")
    if set(sim.replicate) != set(range(1, 13)) or set(sim.templateHp) != {0.1, 0.5}:
        raise ValueError("Simulation repetitions or reference-template families differ")
    sim_metrics = ["waveformRmse", "waveformR", "meanError250_500", "preStimRmse", "postStimRmse"]
    if not np.isfinite(sim[sim_metrics].to_numpy(float)).all():
        raise ValueError("A simulation recovery metric is not finite")
    delta_rows = sim.loc[sim.condition.eq("Delta")]
    if not np.allclose(delta_rows.deltaMeanError250_500, delta_rows.meanError250_500):
        raise ValueError("Simulated Right-minus-left recovery error is inconsistent")
    if len(example) != 2 * 4 * 3 * 257 or example.groupby(
        ["templateHp", "highpass", "condition"]
    ).size().ne(257).any():
        raise ValueError("The recovery example is incomplete")
    if not np.isfinite(example[["reference", "recovered"]].to_numpy(float)).all():
        raise ValueError("The recovery example contains missing waveform values")


def frequency_response(directory: Path, qa_directory: Path | None) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 3.4), layout="constrained")
    frequencies = np.geomspace(0.025, 10, 2500)
    lp = signal.butter(4, 30, btype="lowpass", fs=256, output="sos")
    _, lp_response = signal.sosfreqz(lp, worN=frequencies, fs=256)
    for hp in HP_VALUES:
        hp_sos = signal.butter(2, hp, btype="highpass", fs=256, output="sos")
        _, hp_response = signal.sosfreqz(hp_sos, worN=frequencies, fs=256)
        magnitude = np.abs(hp_response * lp_response) ** 2
        ax.semilogx(
            frequencies,
            20 * np.log10(np.maximum(magnitude, 1e-10)),
            color=HP_COLORS[hp],
            lw=1.8,
            label=f"{hp:g} Hz",
        )
    ax.axhline(-6.0206, lw=0.7, color="#666666", ls=":")
    ax.set(xlim=(0.03, 8), ylim=(-55, 1), xlabel="Frequency (Hz)", ylabel="Two-pass gain (dB)")
    tick_values = (0.05, 0.1, 0.2, 0.5, 1, 2, 5)
    ax.set_xticks(tick_values, [f"{value:g}" for value in tick_values])
    ax.tick_params(axis="x", bottom=False, labelbottom=False, top=True, labeltop=True, pad=5)
    ax.minorticks_off()
    ax.legend(title="High-pass cutoff", ncol=2, loc="lower right")
    finish(fig, directory, "filter_frequency_response", qa_directory)


def impulse_response(directory: Path, qa_directory: Path | None) -> None:
    fs = 256
    time = np.arange(-20 * fs, 20 * fs + 1) / fs
    impulse = np.zeros(time.size)
    impulse[time.size // 2] = 1
    fig, ax = plt.subplots(figsize=(6.8, 3.4), layout="constrained")
    for hp in HP_VALUES:
        bh, ah = signal.butter(2, hp, btype="highpass", fs=fs)
        bl, al = signal.butter(4, 30, btype="lowpass", fs=fs)
        response = signal.filtfilt(bl, al, signal.filtfilt(bh, ah, impulse))
        absolute = np.abs(response)
        visible = np.where(absolute >= 1e-7, absolute, np.nan)
        ax.semilogy(time, visible, lw=1.2,
                    color=HP_COLORS[hp], label=f"{hp:g} Hz")
    ax.set(xlim=(-15, 15), ylim=(1e-7, 1), xlabel="Time relative to impulse (s)",
           ylabel="Absolute impulse response")
    ax.set_yticks([1e-6, 1e-4, 1e-2, 1], ["1e-6", "1e-4", "1e-2", "1"])
    ax.minorticks_off()
    ax.legend(title="High-pass cutoff", ncol=4, loc="lower center",
              bbox_to_anchor=(0.5, 1.02), fontsize=7)
    finish(fig, directory, "filter_impulse_response", qa_directory)


def example_erp(curves: pd.DataFrame, directory: Path, qa_directory: Path | None) -> None:
    frame = curves.loc[
        (curves.record == "A") & (curves.task == 2) & (curves.channel == "F3")
    ].sort_values("timeSec")
    if frame.empty:
        raise ValueError("A Task 2 F3 waveform rows are missing")
    fig, ax = plt.subplots(figsize=(6.8, 3.6), layout="constrained")
    for hp in HP_VALUES:
        part = frame.loc[np.isclose(frame.highpass, hp) & (frame.timeSec >= -0.2)]
        ax.plot(part.timeSec, part.left, color=HP_COLORS[hp], lw=1.3, ls="--")
        ax.plot(part.timeSec, part.right, color=HP_COLORS[hp], lw=1.7, label=f"{hp:g} Hz")
    ax.axvline(0, color="#555555", lw=0.7)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.axvspan(0.25, 0.5, color="#dddddd", alpha=0.3, zorder=0)
    ax.set(xlim=(-0.2, 0.8), xlabel="Time from visual cue (s)", ylabel=f"Response ({SOURCE_UNITS})")
    ax.legend(title="Right: solid; Left: dashed", ncol=4, fontsize=7,
              loc="lower center", bbox_to_anchor=(0.5, 1.01))
    finish(fig, directory, "A_task2_F3_left_right", qa_directory)

    fig, ax = plt.subplots(figsize=(6.8, 3.4), layout="constrained")
    for hp in HP_VALUES:
        part = frame.loc[np.isclose(frame.highpass, hp) & (frame.timeSec >= -0.2)]
        ax.plot(part.timeSec, part.delta, color=HP_COLORS[hp], lw=1.6, label=f"{hp:g} Hz")
    ax.axvline(0, color="#555555", lw=0.7)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.axvspan(0.25, 0.5, color="#dddddd", alpha=0.3, zorder=0)
    ax.set(xlim=(-0.2, 0.8), xlabel="Time from visual cue (s)", ylabel=f"Right − Left ({SOURCE_UNITS})")
    ax.legend(title="High-pass cutoff", ncol=4,
              loc="lower center", bbox_to_anchor=(0.5, 1.01))
    finish(fig, directory, "A_task2_F3_delta", qa_directory)


def full_matrix(real: pd.DataFrame, directory: Path, qa_directory: Path | None) -> None:
    ordered = real.copy()
    ordered["key"] = ordered.record + " · Task " + ordered.task.astype(str) + " · " + ordered.channel
    row_order = [f"{rec} · Task {task} · {channel}" for rec in ("A", "B") for task in (1, 2) for channel in ("F3", "Fz", "F4")]
    matrix = ordered.pivot(index="key", columns="highpass", values="delta_mean_250_500").loc[row_order, list(HP_VALUES)]
    values = matrix.to_numpy(float)
    bound = max(np.nanmax(np.abs(values)), 1)
    fig, ax = plt.subplots(figsize=(6.8, 5.1), layout="constrained")
    im = ax.imshow(values, cmap="RdBu_r", vmin=-bound, vmax=bound, aspect="auto")
    ax.set_xticks(np.arange(4), [f"{hp:g}" for hp in HP_VALUES])
    ax.set_yticks(np.arange(len(row_order)), row_order)
    ax.set_xlabel("High-pass cutoff (Hz)")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            shown = 0.0 if abs(values[i, j]) < 0.05 else values[i, j]
            ax.text(j, i, f"{shown:.1f}", ha="center", va="center", fontsize=6.5,
                    color="white" if abs(values[i, j]) > 0.65 * bound else "#1e293b")
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(f"Right − Left, 250–500 ms ({SOURCE_UNITS})")
    finish(fig, directory, "delta_mean_full_matrix", qa_directory, exclude_axes=[cb.ax])


def reliability_and_baseline(real: pd.DataFrame, directory: Path, qa_directory: Path | None) -> None:
    # Every pale trajectory is a fixed record-task-channel combination; the dark line is their median.
    for column, ylabel, stem in (
        ("split_half_r", "Right − Left split-half correlation", "split_half_by_highpass"),
        ("baseline_rms", f"Residual baseline RMS (−0.2 to 0 s; {SOURCE_UNITS})", "baseline_noise_by_highpass"),
    ):
        fig, ax = plt.subplots(figsize=(6.8, 3.4), layout="constrained")
        for (_, _, _), group in real.groupby(["record", "task", "channel"], sort=True):
            group = group.sort_values("highpass")
            ax.plot(group.highpass, group[column], color="#b4bec7", lw=0.8, alpha=0.8)
        medians = real.groupby("highpass")[column].median().reindex(HP_VALUES)
        ax.plot(HP_VALUES, medians.values, color="#1f4e75", marker="o", lw=2.2, label="Median across 12 combinations")
        ax.set_xticks(HP_VALUES)
        ax.set(xlabel="High-pass cutoff (Hz)", ylabel=ylabel)
        ax.legend(loc="best", fontsize=7)
        finish(fig, directory, stem, qa_directory)


def simulation_plots(
    sim: pd.DataFrame,
    example: pd.DataFrame,
    directory: Path,
    qa_directory: Path | None,
) -> None:
    for metric, ylabel, stem in (
        ("waveformRmse", f"Waveform RMSE ({SOURCE_UNITS})", "simulation_rmse"),
        ("waveformR", "Waveform correlation with reference", "simulation_waveform_correlation"),
        ("deltaMeanError250_500", f"Absolute Right − Left mean error ({SOURCE_UNITS})", "simulation_delta_error"),
    ):
        included_conditions = ["Delta"] if metric.startswith("delta") else ["Left", "Right"]
        relevant = sim.loc[sim.condition.isin(included_conditions)].copy()
        if metric.startswith("delta"):
            relevant[metric] = relevant[metric].abs()
        fig, ax = plt.subplots(figsize=(6.8, 3.4), layout="constrained")
        for template_hp, color in ((0.1, "#285c9a"), (0.5, "#a5393b")):
            rows = relevant.loc[np.isclose(relevant.templateHp, template_hp)]
            if rows.empty:
                continue
            unit_means = rows.groupby(
                ["record", "task", "channel", "condition", "highpass"], as_index=False
            )[metric].median()
            summary = unit_means.groupby("highpass")[metric].agg(
                ["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]
            ).reindex(HP_VALUES)
            center, lo, hi = (summary.iloc[:, k].to_numpy(float) for k in range(3))
            ax.errorbar(HP_VALUES, center, yerr=[center - lo, hi - center], color=color,
                        marker="o", capsize=2.5, lw=1.6, label=f"Reference built at {template_hp:g} Hz")
        ax.set_xticks(HP_VALUES)
        ax.set(xlabel="Applied high-pass cutoff (Hz)", ylabel=ylabel)
        n_units = 12 if metric.startswith("delta") else 24
        ax.legend(title=f"Median across {n_units} combinations; IQR", fontsize=7)
        finish(fig, directory, stem, qa_directory)

    e = example.copy()
    e = e.loc[e.condition.eq("Delta") & np.isclose(e.templateHp, 0.1)]
    if e.empty:
        raise ValueError("The 0.1-Hz template recovery example is missing")
    e = e.loc[(e.record == e.record.iloc[0]) & (e.task == e.task.iloc[0]) & (e.channel == e.channel.iloc[0]) & (e.replicate == e.replicate.iloc[0])]
    fig, ax = plt.subplots(figsize=(6.8, 3.4), layout="constrained")
    reference = e.loc[np.isclose(e.highpass, HP_VALUES[0])].sort_values("timeSec")
    ax.plot(reference.timeSec, reference.reference, color="#242b33", lw=2, label="Semi-synthetic reference")
    ax.plot(reference.timeSec, reference.noisyBeforeHighpass, color="#7c8995", lw=1.2,
            ls="--", label="Noisy before high-pass")
    for hp in HP_VALUES:
        part = e.loc[np.isclose(e.highpass, hp)].sort_values("timeSec")
        ax.plot(part.timeSec, part.recovered, color=HP_COLORS[hp], lw=1.2, label=f"Recovered, {hp:g} Hz")
    ax.axvline(0, color="#666666", lw=0.7)
    ax.set(xlim=(-0.2, 0.8), xlabel="Time from visual cue (s)", ylabel=f"Right − Left ({SOURCE_UNITS})")
    ax.legend(ncol=2, fontsize=7)
    finish(fig, directory, "simulation_recovery_example", qa_directory)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--qa-dir", type=Path, default=None,
                        help="Temporary directory for QA-only PDFs and alignment overlays")
    args = parser.parse_args()
    out = args.output_dir.resolve()
    directory = out / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    qa_directory = args.qa_dir.resolve() if args.qa_dir is not None else None
    if qa_directory is not None:
        qa_directory.mkdir(parents=True, exist_ok=True)
    real = pd.read_csv(out / "q1_highpass_sensitivity_real.csv")
    curves = pd.read_csv(out / "q1_highpass_waveforms.csv")
    sim = pd.read_csv(out / "q1_filter_simulation_recovery.csv")
    example = pd.read_csv(out / "q1_filter_simulation_example_waveform.csv")
    validate_sources(real, curves, sim, example)
    frequency_response(directory, qa_directory)
    impulse_response(directory, qa_directory)
    example_erp(curves, directory, qa_directory)
    full_matrix(real, directory, qa_directory)
    reliability_and_baseline(real, directory, qa_directory)
    simulation_plots(sim, example, directory, qa_directory)
    print(f"Wrote high-pass figures to {directory}")


if __name__ == "__main__":
    main()
