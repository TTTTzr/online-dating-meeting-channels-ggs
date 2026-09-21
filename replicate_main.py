#!/usr/bin/env python3
"""Reproduce Figure 1 from the harmonized GGS-II Wave 1 Stata releases."""

import argparse
import hashlib
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd


RELEASES = [
    ("AR", "1_2", "Argentina (Buenos Aires)", 935),
    ("AT", "1_1", "Austria", 4646),
    ("CZ", "2_1", "Czechia", 3431),
    ("EE", "2_2", "Estonia", 5874),
    ("FR", "1_0", "France", 6546),
    ("HK", "1_0", "Hong Kong", 1238),
    ("ISL", "1_0", "Iceland", 1693),
    ("NL", "1_0", "Netherlands", 3940),
    ("TW", "1_0", "Taiwan", 2418),
    ("UK", "1_2", "United Kingdom", 4632),
    ("UY", "1_4", "Uruguay", 3419),
]
PERIODS = [2002, 2005, 2008, 2011, 2014, 2017, 2020]
LABELS = ["2002–04", "2005–07", "2008–10", "2011–13", "2014–16", "2017–19", "2020–21"]
CHANNELS = ["Online dating", "Other online", "Friends/family", "School/work", "Public/leisure", "Other"]
RECODE = {1: 3, 2: 3, 3: 4, 4: 0, 5: 1, 6: 4, 7: 4, 8: 4, 9: 4, 10: 2, 11: 2, 12: 5}
COLORS = ["#C84C3C", "#7565A7", "#5F7F72", "#4F718B", "#9A7442", "#7A7A7A"]
BOOTSTRAPS = 500
SEED = 8262026


def find_files(data_dir):
    files = []
    for code, version, country, expected_n in RELEASES:
        filename = f"GGSII_Wave1_{code}_V_{version}.dta"
        matches = sorted(data_dir.rglob(filename))
        if len(matches) != 1:
            raise FileNotFoundError(f"Expected exactly one {filename} under {data_dir}; found {len(matches)}.")
        files.append((matches[0], code, country, expected_n))
    return files


def build_context(path, code, country):
    required = ["respid", "weight", "intdatey", "dem02y", "dem30a", "dem30by", "dem22a"]
    with pd.read_stata(path, iterator=True, convert_categoricals=False) as reader:
        available = set(reader.variable_labels())
        missing = set(required) - available
        if missing:
            raise ValueError(f"{path.name}: missing variables {sorted(missing)}")
        histories = [(f"lhi04_y{i}", f"lhi04a_{i}") for i in range(1, 21)]
        histories = [(y, c) for y, c in histories if y in available and c in available]
        if not histories:
            raise ValueError(f"{path.name}: no prior-union meeting histories")
        columns = required + [v for pair in histories for v in pair]
        raw = reader.read(columns=columns, convert_categoricals=False, preserve_dtypes=False)
    if raw.respid.isna().any():
        raise ValueError(f"{path.name}: respondent IDs are missing")
    for col in columns:
        if col != "respid":
            raw[col] = pd.to_numeric(raw[col], errors="coerce").mask(raw[col] > 1_000_000)

    base = pd.DataFrame({
        "country": country,
        "respondent_id": code + "_" + raw.respid.astype(str),
        "weight": raw.weight,
        "interview_year": raw.intdatey,
        "birth_year": raw.dem02y,
    })
    current = base.assign(union_year=raw.dem30by, channel_index=raw.dem22a.map(RECODE))
    rows = [current.loc[raw.dem30a.eq(1)]]
    for year, channel in histories:
        rows.append(base.assign(union_year=raw[year], channel_index=raw[channel].map(RECODE)))
    data = pd.concat(rows, ignore_index=True)
    age = data.union_year - data.birth_year
    keep = (data.union_year.between(2002, 2021) & age.between(15, 69)
            & data.channel_index.notna() & data.weight.gt(0))
    data = data.loc[keep].copy()
    data["period_index"] = ((data.union_year - 2002) // 3).astype(int)
    # Only the 2020–21 period requires interviews in 2022 or later.
    data = data.loc[data.period_index.lt(6) | data.interview_year.ge(2022)].copy()
    data["channel_index"] = data.channel_index.astype(int)
    return data.reset_index(drop=True)


def weighted_shares(cells, weights):
    totals = np.bincount(cells, weights=weights, minlength=42).reshape(7, 6)
    denominators = totals.sum(axis=1, keepdims=True)
    if (denominators == 0).any():
        raise ValueError("A period has zero total weight.")
    return totals / denominators


def estimate_context(data, seed):
    cells = (6 * data.period_index + data.channel_index).to_numpy()
    weights = data.weight.to_numpy(dtype=float)
    point = weighted_shares(cells, weights)
    ids, cluster = np.unique(data.respondent_id.to_numpy(), return_inverse=True)
    rng = np.random.default_rng(seed)
    draws = np.empty((BOOTSTRAPS, 7, 6))
    # One multiplier per respondent, shared by all of their unions and periods.
    for b in range(BOOTSTRAPS):
        multipliers = rng.poisson(1.0, len(ids))
        draws[b] = weighted_shares(cells, weights * multipliers[cluster])
    low, high = np.quantile(draws, [0.025, 0.975], axis=0)
    change = point[5] - point[0]
    change_low, change_high = np.quantile(draws[:, 5] - draws[:, 0], [0.025, 0.975], axis=0)
    country = data.country.iloc[0]
    levels = [{"country": country, "cohort_start": year, "cohort": LABELS[p],
               "channel": channel, "estimate": point[p, c], "lower": low[p, c], "upper": high[p, c]}
              for p, year in enumerate(PERIODS) for c, channel in enumerate(CHANNELS)]
    changes = [{"country": country, "channel": channel, "estimate": change[c],
                "lower": change_low[c], "upper": change_high[c]}
               for c, channel in enumerate(CHANNELS)]
    return levels, changes


def draw_figure(levels, changes, output, sample_n, covid_n):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8.8, "axes.titlesize": 9.2, "axes.labelsize": 9.2,
        "xtick.labelsize": 7.8, "ytick.labelsize": 7.8, "axes.linewidth": 0.6,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    contexts = sorted(levels.country.unique())
    online = levels.loc[levels.channel.eq("Online dating")]
    fig = plt.figure(figsize=(14, 8.4), facecolor="white")
    outer = GridSpec(2, 1, figure=fig, height_ratios=[0.45, 0.55], top=0.845,
                     bottom=0.115, left=0.105, right=0.985, hspace=0.46)
    grid_a = outer[0].subgridspec(2, 6, wspace=0.30, hspace=0.64)
    y_max = np.ceil(max(0.35, online.upper.max() * 1.08) / 0.05) * 0.05
    for i, country in enumerate(contexts):
        ax = fig.add_subplot(grid_a[i // 6, i % 6])
        sub = online.loc[online.country.eq(country)].sort_values("cohort_start")
        x = np.arange(7)
        est, low, high = (sub[col].to_numpy() for col in ["estimate", "lower", "upper"])
        ax.fill_between(x[:6], low[:6], high[:6], color=COLORS[0], alpha=0.13, linewidth=0)
        ax.plot(x[:6], est[:6], color=COLORS[0], linewidth=1.7, zorder=2)
        ax.scatter(x[:6], est[:6], color=COLORS[0], s=15, edgecolor="white", linewidth=0.45, zorder=3)
        ax.plot(x[5:], est[5:], color=COLORS[0], linewidth=1.25, linestyle=(0, (2.2, 1.8)), zorder=3)
        ax.errorbar([x[6]], [est[6]], yerr=[[est[6] - low[6]], [high[6] - est[6]]],
                    fmt="none", ecolor=COLORS[0], elinewidth=0.8, capsize=1.8, zorder=3)
        ax.scatter(x[6], est[6], marker="D", s=28, facecolor="white", edgecolor=COLORS[0], linewidth=1.1, zorder=4)
        ax.set_title(country, loc="left", pad=2, fontweight="semibold")
        ax.set_xlim(-0.18, 6.18)
        ax.set_ylim(0, y_max)
        ax.set_xticks(x, ["02", "05", "08", "11", "14", "17", "20"])
        ticks = np.arange(0, y_max + 0.001, 0.10)
        ax.set_yticks(ticks, [f"{v:.0%}" for v in ticks] if i % 6 == 0 else [""] * len(ticks))
        ax.grid(axis="y", color="#D8D8D8", linewidth=0.45)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#8A8A8A")
        ax.tick_params(length=2.2, color="#8A8A8A", pad=1.5)

    # Other stays in the denominator, but is not displayed in Panel B.
    grid_b = outer[1].subgridspec(1, 5, wspace=0.14)
    visible = changes.loc[changes.channel.ne("Other")]
    bound = max(10, int(np.ceil(np.abs(visible[["lower", "upper"]].to_numpy() * 100).max() / 5) * 5))
    y = np.arange(len(contexts))[::-1]
    titles = ["Online\ndating", "Other\nonline", "Friends /\nfamily", "School /\nwork", "Public /\nleisure"]
    for i, channel in enumerate(CHANNELS[:5]):
        ax = fig.add_subplot(grid_b[0, i])
        sub = visible.loc[visible.channel.eq(channel)].set_index("country").loc[contexts]
        est, low, high = (sub[col].to_numpy() * 100 for col in ["estimate", "lower", "upper"])
        ax.axvline(0, color="#A8A8A8", linewidth=0.8, zorder=0)
        for pos in y:
            ax.axhline(pos, color="#ECECEC", linewidth=0.4, zorder=0)
        ax.hlines(y, low, high, color=COLORS[i], linewidth=1.2, alpha=0.88, zorder=1)
        ax.scatter(est, y, s=25 if i == 0 else 19, color=COLORS[i], edgecolor="white", linewidth=0.45, zorder=2)
        ax.set_xlim(-bound, bound)
        ax.set_ylim(-0.65, len(contexts) - 0.35)
        ax.set_xticks([-bound, 0, bound], [str(-bound), "0", f"+{bound}"])
        ax.set_yticks(y, contexts if i == 0 else [""] * len(contexts))
        ax.set_title(titles[i], color=COLORS[i], pad=5, fontweight="semibold")
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set_color("#8A8A8A")
        ax.tick_params(axis="y", length=0, pad=4)
        ax.tick_params(axis="x", length=2.2, color="#8A8A8A", pad=2)

    fig.text(0.02, 0.965, "As Online Dating Grew, Which Meeting Channels Lost Share among New Unions?",
             fontsize=19.5, fontweight="semibold", ha="left")
    fig.text(0.02, 0.928, "Survey-weighted shares of co-residential unions beginning in each period, by how partners first met",
             fontsize=11.2, ha="left", color="#4D4D4D")
    fig.text(0.02, 0.873, "A", fontsize=13, fontweight="bold", ha="left")
    fig.text(0.047, 0.873, "Share of unions whose partners first met through online dating", fontsize=11, fontweight="semibold", ha="left")
    fig.text(0.985, 0.873, "Period start year (02 = 2002–04)  •  Weighted share", fontsize=8.4, ha="right", color="#595959")
    fig.text(0.02, 0.485, "B", fontsize=13, fontweight="bold", ha="left")
    fig.text(0.047, 0.485, "Weighted percentage-point change in meeting channels, 2017–19 minus 2002–04",
             fontsize=10.8, fontweight="semibold", ha="left")
    fig.text(0.56, 0.079, "Percentage-point change", fontsize=8.7, ha="center", color="#4F4F4F")
    note = (f"Note: N = {sample_n:,} unions in 11 countries and regions; the 2020–21 period contains N = {covid_n:,}. "
            "Period is defined by the start of co-residence. Points are GGS-weighted estimates; bands/bars are "
            "respondent-clustered 95% Poisson-bootstrap intervals. The 2020–21 period is shown only in Panel A; "
            "Panel B uses pre-pandemic periods. Residual ‘Other’ is included in all denominators but omitted from Panel B; "
            "full decomposition appears in the supplement. Online dating is always separate from other online meeting.")
    fig.text(0.02, 0.039, note, fontsize=7.35, ha="left", va="bottom", color="#4F4F4F", wrap=True)
    fig.text(0.02, 0.014, "Source: Generations and Gender Survey, Round 2, Wave 1; author’s calculations.",
             fontsize=7.6, color="#4F4F4F")
    fig.savefig(output / "figure_1.png", dpi=300, facecolor="white")
    fig.savefig(output / "figure_1.pdf", facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True, help="Directory containing the 11 GGS-II .dta releases (subfolders allowed)")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "output")
    args = parser.parse_args()
    files = find_files(args.data_dir)
    samples, levels, changes, checksums = [], [], [], []
    for i, (path, code, country, expected_n) in enumerate(files):
        data = build_context(path, code, country)
        if len(data) != expected_n:
            raise ValueError(f"{country}: expected {expected_n} unions, found {len(data)}; check the data release.")
        if data.loc[data.period_index.eq(6)].shape[0] < 100:
            raise ValueError(f"{country}: fewer than 100 endpoint unions")
        a, b = estimate_context(data, SEED + i * 10007)
        samples.append(data)
        levels.extend(a)
        changes.extend(b)
        checksums.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
        print(f"{country}: {len(data):,} unions")
    data = pd.concat(samples, ignore_index=True)
    levels = pd.DataFrame(levels).sort_values(["country", "cohort_start", "channel"])
    changes = pd.DataFrame(changes).sort_values(["country", "channel"])
    multiplicity = data.groupby("respondent_id").size()
    covid_n = int(data.period_index.eq(6).sum())
    observed = (len(data), len(multiplicity), int((multiplicity > 1).sum()), int(multiplicity[multiplicity > 1].sum()), covid_n)
    if observed != (38772, 31775, 5614, 12611, 3374):
        raise ValueError(f"Sample counts differ from R2: {observed}")
    np.testing.assert_allclose(levels.groupby(["country", "cohort_start"]).estimate.sum(), 1, atol=1e-12)
    np.testing.assert_allclose(changes.groupby("country").estimate.sum(), 0, atol=1e-12)

    data["cohort"] = data.period_index.map(dict(enumerate(LABELS)))
    counts = data.groupby(["country", "cohort"]).agg(unions=("respondent_id", "size"), respondents=("respondent_id", "nunique")).reset_index()
    findings = []
    for country, sub in changes.groupby("country", sort=True):
        sub = sub.set_index("channel")
        offline = sub.loc[["Friends/family", "School/work", "Public/leisure"], "estimate"]
        findings.append({"country": country, "online_dating_change_pp": 100 * sub.loc["Online dating", "estimate"],
                         "largest_offline_decline": offline.idxmin(), "largest_offline_decline_pp": 100 * offline.min()})
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    levels.to_csv(out / "main_figure_shares.csv", index=False)
    changes.to_csv(out / "main_figure_changes.csv", index=False)
    counts.to_csv(out / "sample_counts.csv", index=False)
    pd.DataFrame(findings).to_csv(out / "findings.csv", index=False)
    draw_figure(levels, changes, out, len(data), covid_n)
    summary = (f"Figure 1 replication\nUnions: {len(data):,}\nRespondents: {len(multiplicity):,}\n"
               f"Respondents with multiple unions: {(multiplicity > 1).sum():,}\n"
               f"Unions from these respondents: {multiplicity[multiplicity > 1].sum():,}\n"
               f"2020–21 unions: {covid_n:,}\nPoisson bootstrap replicates: {BOOTSTRAPS}\nBase seed: {SEED}\n"
               f"Python {platform.python_version()}; NumPy {np.__version__}; pandas {pd.__version__}; Matplotlib {matplotlib.__version__}\n")
    (out / "run_summary.txt").write_text(summary, encoding="utf-8")
    (out / "input_files.sha256").write_text("\n".join(checksums) + "\n", encoding="utf-8")
    print(summary)
    print(f"Outputs: {out.resolve()}")


if __name__ == "__main__":
    main()
