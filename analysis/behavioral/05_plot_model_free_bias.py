"""Create three publication-ready paired plots of model-free temporal bias.

Run from the repository root:
    python analysis/behavioral/06_plot_model_free_bias.py

Reads the participant-level and test outputs of 05_model_free_bias.py.
Writes one PNG and one vector PDF per group to results/behavioral.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "behavioral"
SUBJECT_CSV = OUT / "model_free_bias_subject_level.csv"
TEST_CSV = OUT / "model_free_bias_group_tests.csv"

GROUPS = (
    ("fMRI", "fMRI group"),
    ("replication", "Replication group"),
    ("control", "No-encoding group"),
)
PANELS = (
    ("directional_D", "early_inward", "late_inward", "Directional inward displacement"),
    ("magnitude_D", "early_abs_mean_error", "late_abs_mean_error", "Absolute condition-mean bias"),
)

INK = "#24313D"
GRID = "#E0E5E9"
LINE = "#AAB5BE"
EARLY = "#2E678F"
LATE = "#C67B4A"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def format_p(value: float) -> str:
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def check_test(values: np.ndarray, published: dict[str, str]) -> None:
    result = stats.ttest_1samp(values, popmean=0)
    if not math.isclose(float(result.statistic), float(published["t"]), abs_tol=1e-9):
        raise ValueError("Plotted data do not reproduce the stored t statistic")
    if not math.isclose(float(result.pvalue), float(published["p_two_sided"]), abs_tol=1e-12):
        raise ValueError("Plotted data do not reproduce the stored p value")


def mean_ci(values: np.ndarray) -> tuple[float, float, float]:
    mean = float(np.mean(values))
    sem = float(stats.sem(values))
    bound = float(stats.t.ppf(0.975, len(values) - 1) * sem)
    return mean, mean - bound, mean + bound


def plot_panel(ax, rows: list[dict[str, str]], metric: str, early_col: str,
               late_col: str, title: str, test: dict[str, str]) -> None:
    early = np.array([float(row[early_col]) for row in rows])
    late = np.array([float(row[late_col]) for row in rows])
    differences = early - late
    check_test(differences, test)
    if not np.allclose(differences, [float(row[metric]) for row in rows], atol=1e-10):
        raise ValueError(f"{metric}: participant differences disagree with stored D")

    # One stable horizontal offset per participant prevents overplotting while
    # retaining each person's early-to-late pairing.
    offsets = np.linspace(-0.095, 0.095, len(rows))
    for index in range(len(rows)):
        ax.plot([offsets[index], 1 + offsets[index]],
                [early[index], late[index]], color=LINE, linewidth=0.9,
                alpha=0.58, zorder=1)
    ax.scatter(offsets, early, s=27, facecolor=EARLY, edgecolor="white",
               linewidth=0.35, alpha=0.9, zorder=2)
    ax.scatter(1 + offsets, late, s=27, facecolor=LATE, edgecolor="white",
               linewidth=0.35, alpha=0.9, zorder=2)

    for x, values, color in ((0, early, EARLY), (1, late, LATE)):
        mean, low, high = mean_ci(values)
        ax.errorbar(x, mean, yerr=[[mean - low], [high - mean]], fmt="D",
                    markersize=7.5, color=color, markeredgecolor=INK,
                    markeredgewidth=0.9, ecolor=INK, elinewidth=1.6,
                    capsize=5, capthick=1.6, zorder=5)

    ax.set_title(title, loc="left", fontsize=11.5, color=INK, weight="semibold", pad=12)
    ax.text(0.98, 0.975,
            f"Early − Late = {float(test['mean']):.2f} pp\n"
            f"95% CI [{float(test['ci95_low']):.2f}, {float(test['ci95_high']):.2f}]\n"
            f"t({int(float(test['df']))}) = {float(test['t']):.2f}, p = {format_p(float(test['p_two_sided']))}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9.2,
            color=INK, linespacing=1.4,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 3})
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(0, 32)
    ax.set_xticks([0, 1], ["Early", "Late"])
    ax.set_yticks([0, 5, 10, 15, 20, 25])
    ax.set_ylabel("Percentage points", fontsize=10, color=INK)
    ax.tick_params(axis="both", colors=INK, labelsize=9.5, length=0)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#9AA7B1")


def plot_group(group: str, display_name: str, subjects: list[dict[str, str]],
               tests: dict[tuple[str, str], dict[str, str]]) -> None:
    rows = sorted((row for row in subjects if row["group"] == group),
                  key=lambda row: row["subject"])
    if not rows:
        raise ValueError(f"No participant rows for {group}")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.1), dpi=150)
    fig.patch.set_facecolor("white")
    for ax, (metric, early_col, late_col, title) in zip(axes, PANELS):
        plot_panel(ax, rows, metric, early_col, late_col, title, tests[(group, metric)])

    fig.suptitle(f"{display_name}  ·  n = {len(rows)}", x=0.065, y=0.985,
                 ha="left", fontsize=15, weight="semibold", color=INK)
    fig.text(0.065, 0.905,
             "Raw signed-error analysis · two-sided paired t tests on participant means",
             ha="left", fontsize=9.6, color="#52616D")
    fig.text(0.065, 0.035,
             "Each connected pair is one participant. Diamonds and bars show condition means and 95% CIs; "
             "test CIs refer to paired differences.",
             ha="left", fontsize=8.6, color="#52616D")
    fig.subplots_adjust(left=0.075, right=0.97, top=0.77, bottom=0.15, wspace=0.22)

    stem = OUT / f"model_free_bias_{group}"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.18)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    print(f"Created {stem.with_suffix('.png').name} and {stem.with_suffix('.pdf').name}")


def main() -> None:
    subjects = read_csv(SUBJECT_CSV)
    test_rows = read_csv(TEST_CSV)
    tests = {(row["group"], row["metric"]): row for row in test_rows}
    for group, name in GROUPS:
        plot_group(group, name, subjects, tests)


if __name__ == "__main__":
    main()
