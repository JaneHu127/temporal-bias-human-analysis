"""Model-free early-versus-late temporal bias from raw signed errors.

Run from the repository root with ``python analysis/behavioral/05_model_free_bias.py``.
Requires SciPy. Source data are read-only; outputs go to results/behavioral.

The participant is the unit of inference. Signed error is reported minus
target location, in percentage points. The two prespecified contrasts are:

  directional_D = (mean_e20 + mean_e40 + mean_e60 + mean_e80) / 2
  magnitude_D   = (|mean_e20| + |mean_e40| - |mean_e60| - |mean_e80|) / 2

The first compares inward displacement with its direction retained. The
second parallels the manuscript's early-versus-late absolute-mu contrast,
but uses absolute values of *condition means*, not trial-wise absolute error.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "behavioral"
LOCATIONS = (20, 40, 60, 80)
GROUPS = (
    ("fMRI", "fMRI_group/fmri_group_trials.csv", "temporal_location(%)", "reported_location(%)", 18),
    ("replication", "replication_group/replication_group_trials.csv", "target_location", "reported_location", 30),
    ("control", "control_group/control_group_trials.csv", "temporal_location(%)", "reported_location(%)", 9),
)
METRICS = (
    ("directional_D", "Early minus late inward displacement"),
    ("magnitude_D", "Early minus late absolute condition-mean bias"),
)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_group(group: str, relative: str, location_col: str, response_col: str, expected_n: int):
    path = ROOT / "data" / "behavioral" / relative
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"subject", "trial", "interval_width", location_col, response_col, "signed_error"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing columns {sorted(required - set(reader.fieldnames or []))}")
        rows = list(reader)

    by_cell: dict[tuple[str, int], list[float]] = defaultdict(list)
    seen_trials: set[tuple[str, str]] = set()
    subjects: set[str] = set()
    invalid = 0
    for row in rows:
        subject = row["subject"].strip()
        subjects.add(subject)
        key = (subject, row["trial"].strip())
        if key in seen_trials:
            raise ValueError(f"{path}: duplicate subject/trial {key}")
        seen_trials.add(key)

        valid_flag = row.get("valid_response", "True").strip().lower()
        if valid_flag not in {"true", "false"}:
            raise ValueError(f"{path}: unexpected valid_response {valid_flag!r}")
        if valid_flag == "false":
            if row["signed_error"].strip():
                raise ValueError(f"{path}: invalid response has signed_error at {key}")
            invalid += 1
            continue

        if not row["signed_error"].strip():
            raise ValueError(f"{path}: valid response missing signed_error at {key}")
        location = float(row[location_col])
        if location not in LOCATIONS:
            raise ValueError(f"{path}: unexpected temporal location {location} at {key}")
        width = float(row["interval_width"])
        if width not in (12, 60, 300):
            raise ValueError(f"{path}: unexpected interval width {width} at {key}")
        error = float(row["signed_error"])
        response = float(row[response_col])
        if not all(map(math.isfinite, (location, width, error, response))):
            raise ValueError(f"{path}: non-finite numeric value at {key}")
        if not 0 <= response <= 100:
            raise ValueError(f"{path}: out-of-range response at {key}")
        if not math.isclose(error, response - location, abs_tol=1e-5):
            raise ValueError(f"{path}: signed_error != response - target at {key}")
        by_cell[(subject, int(location))].append(error)

    if len(subjects) != expected_n:
        raise ValueError(f"{path}: expected {expected_n} subjects; found {len(subjects)}")
    if sum(map(len, by_cell.values())) + invalid != len(rows):
        raise ValueError(f"{path}: valid/invalid row accounting failed")

    output = []
    for subject in sorted(subjects):
        cells = {loc: by_cell[(subject, loc)] for loc in LOCATIONS}
        if any(not values for values in cells.values()):
            raise ValueError(f"{path}: {subject} is missing a temporal location")
        means = {loc: sum(values) / len(values) for loc, values in cells.items()}
        early_inward = (means[20] + means[40]) / 2
        late_inward = -(means[60] + means[80]) / 2
        early_magnitude = (abs(means[20]) + abs(means[40])) / 2
        late_magnitude = (abs(means[60]) + abs(means[80])) / 2
        record = {"group": group, "subject": subject, "n_valid_trials": sum(map(len, cells.values()))}
        record.update({f"n_{loc}": len(cells[loc]) for loc in LOCATIONS})
        record.update({f"mean_error_{loc}": means[loc] for loc in LOCATIONS})
        record.update({
            "early_inward": early_inward,
            "late_inward": late_inward,
            "directional_D": early_inward - late_inward,
            "early_abs_mean_error": early_magnitude,
            "late_abs_mean_error": late_magnitude,
            "magnitude_D": early_magnitude - late_magnitude,
            "mirror_20_80": means[20] + means[80],
            "mirror_40_60": means[40] + means[60],
        })
        output.append(record)

    audit = {
        "group": group,
        "source_file": str(path.relative_to(ROOT)).replace("\\", "/"),
        "n_source_rows": len(rows),
        "n_excluded_invalid": invalid,
        "n_valid_rows": len(rows) - invalid,
        "n_subjects": len(subjects),
    }
    return output, audit


def one_sample(group: str, metric: str, values: list[float]) -> dict:
    n = len(values)
    mean = float(stats.tmean(values))
    sd = float(stats.tstd(values))
    sem = sd / math.sqrt(n)
    t_value = mean / sem
    df = n - 1
    p = float(2 * stats.t.sf(abs(t_value), df))
    critical = float(stats.t.ppf(0.975, df))
    return {
        "group": group, "metric": metric, "n": n, "mean": mean, "sd": sd, "sem": sem,
        "ci95_low": mean - critical * sem, "ci95_high": mean + critical * sem,
        "t": t_value, "df": df, "p_two_sided": p, "cohens_dz": mean / sd,
    }


def welch(metric: str, group_a: str, group_b: str, a: list[float], b: list[float]) -> dict:
    mean_a, mean_b = float(stats.tmean(a)), float(stats.tmean(b))
    var_a, var_b = float(stats.tvar(a)), float(stats.tvar(b))
    term_a, term_b = var_a / len(a), var_b / len(b)
    se = math.sqrt(term_a + term_b)
    difference = mean_a - mean_b
    df = (term_a + term_b) ** 2 / (term_a**2 / (len(a) - 1) + term_b**2 / (len(b) - 1))
    t_value = difference / se
    p = float(2 * stats.t.sf(abs(t_value), df))
    critical = float(stats.t.ppf(0.975, df))
    return {
        "metric": metric, "group_a": group_a, "group_b": group_b,
        "n_a": len(a), "n_b": len(b), "mean_difference": difference,
        "ci95_low": difference - critical * se, "ci95_high": difference + critical * se,
        "t": t_value, "df": df, "p_two_sided": p,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    subject_rows: list[dict] = []
    audits: list[dict] = []
    for spec in GROUPS:
        participants, audit = load_group(*spec)
        subject_rows.extend(participants)
        audits.append(audit)

    subject_path = OUT / "model_free_bias_subject_level.csv"
    write_csv(subject_path, subject_rows, list(subject_rows[0]))
    write_csv(OUT / "model_free_bias_data_audit.csv", audits, list(audits[0]))

    tests = []
    comparisons = []
    for metric, _ in METRICS:
        by_group = {
            group: [float(row[metric]) for row in subject_rows if row["group"] == group]
            for group, *_ in GROUPS
        }
        for group in by_group:
            tests.append(one_sample(group, metric, by_group[group]))
        for group in ("fMRI", "replication"):
            comparisons.append(welch(metric, group, "control", by_group[group], by_group["control"]))

    # Holm adjustment for the two encoding-versus-control comparisons within
    # each metric. Both raw and adjusted p values are retained.
    for metric, _ in METRICS:
        items = sorted((r for r in comparisons if r["metric"] == metric), key=lambda r: r["p_two_sided"])
        adjusted = 0.0
        for index, row in enumerate(items):
            adjusted = max(adjusted, min(1.0, (len(items) - index) * row["p_two_sided"]))
            row["p_holm_within_metric"] = adjusted

    write_csv(OUT / "model_free_bias_group_tests.csv", tests, list(tests[0]))
    write_csv(OUT / "model_free_bias_group_comparisons.csv", comparisons, list(comparisons[0]))

    lines = [
        "Model-free asymmetry of raw signed temporal error",
        "================================================",
        "Signed error = reported location - target location (percentage points).",
        "Analysis unit: participant. Each condition mean is calculated from valid trials.",
        "All tests are two-sided. Group tests are one-sample t tests of D = 0.",
        "Group comparisons are Welch t tests; Holm correction covers the two",
        "encoding-versus-control comparisons separately for each metric.",
        "",
        "directional_D = early inward displacement - late inward displacement",
        "              = (mean_e20 + mean_e40 + mean_e60 + mean_e80) / 2",
        "magnitude_D = (|mean_e20| + |mean_e40| - |mean_e60| - |mean_e80|) / 2",
        "The absolute value is applied AFTER averaging signed errors within",
        "participant and condition. It is not mean trial-wise absolute error.",
        "Magnitude_D parallels |mu| but is in raw percentage points, whereas",
        "MemToolbox mu is fitted to transformed angular errors.",
        "",
        "Data audit:",
    ]
    for row in audits:
        lines.append(
            f"  {row['group']}: {row['n_subjects']} participants, {row['n_valid_rows']} valid trials, "
            f"{row['n_excluded_invalid']} excluded invalid trials"
        )
    lines += ["", "Condition-mean sign checks (opposite to predicted inward direction):"]
    for group, *_ in GROUPS:
        members = [r for r in subject_rows if r["group"] == group]
        counts = {
            loc: sum(
                r[f"mean_error_{loc}"] < 0 if loc in (20, 40)
                else r[f"mean_error_{loc}"] > 0
                for r in members
            )
            for loc in LOCATIONS
        }
        lines.append(
            f"  {group}: 20%={counts[20]}, 40%={counts[40]}, "
            f"60%={counts[60]}, 80%={counts[80]} participants"
        )
    for metric, description in METRICS:
        lines += ["", f"{metric}: {description}"]
        for row in (r for r in tests if r["metric"] == metric):
            lines.append(
                f"  {row['group']}: n={row['n']}, mean={row['mean']:.3f} pp, "
                f"SD={row['sd']:.3f}, 95% CI [{row['ci95_low']:.3f}, {row['ci95_high']:.3f}], "
                f"t({row['df']})={row['t']:.3f}, p={row['p_two_sided']:.6g}, "
                f"Cohen's dz={row['cohens_dz']:.3f}"
            )
        lines.append("  Encoding versus no-encoding control:")
        for row in (r for r in comparisons if r["metric"] == metric):
            lines.append(
                f"    {row['group_a']} - control: difference={row['mean_difference']:.3f} pp, "
                f"95% CI [{row['ci95_low']:.3f}, {row['ci95_high']:.3f}], "
                f"Welch t({row['df']:.2f})={row['t']:.3f}, "
                f"p={row['p_two_sided']:.6g}, Holm p={row['p_holm_within_metric']:.6g}"
            )
    lines += [
        "",
        "Interpretation notes:",
        "  Positive directional_D supports stronger early inward displacement.",
        "  Positive magnitude_D supports stronger early absolute mean bias.",
        "  Magnitude_D discards direction; inspect condition-mean signs above.",
        "  A nonsignificant result in the n=9 control group does not establish absence.",
        "  Encoding-versus-control comparisons are exploratory, and the small",
        "  control sample makes their Welch normality approximation uncertain.",
        "  Neither encoding-versus-control comparison reaches p < .05 after",
        "  Holm adjustment within either metric.",
        "  These aggregate tests do not adjust for interval width; the design",
        "  counterbalanced location and width, and within-width checks can supplement them.",
        "",
    ]
    (OUT / "model_free_bias_report.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
