#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------
# matplotlib setup
# ------------------------------------------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18

# ------------------------------------------------------------
# utility functions
# ------------------------------------------------------------
def read_simulation_data(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path).fillna("")
    if "actual_hinge_cnt" not in df.columns:
        df["actual_hinge_cnt"] = df["k"]
    df["primary_key"] = (
        df["p_pdb_id"].astype(str).str[:9]
        + "_hinge_"
        + df["actual_hinge_cnt"].astype(str)
        + "_sigma0.5"
    )
    return df

def read_dyndom_data(file_path: str, k: int | None = None) -> pd.DataFrame:
    df = pd.read_csv(file_path).fillna("")
    df["primary_key"] = df["p_pdb_id"] + "_" + df["q_pdb_id"]
    if k is not None:
        df["k"] = k
    return df

# ------------------------------------------------------------
# main pipeline
# ------------------------------------------------------------
def main() -> None:
    os.makedirs("figures", exist_ok=True)

    # ======================
    # 1. Simulation dataset
    # ======================
    sim_shibuya_df = read_simulation_data("rmsdh_result/simulation_shibuya_combined.csv")
    sim_sh_df       = read_simulation_data("rmsdh_result/simulation_sh_combined.csv")

    sim_monge = pd.concat([
        read_simulation_data(f"delta_g_simulation_{k}_sigma0.5.csv")
        for k in range(2, 6)
    ])

    def _merge(df):
        df = df.merge(sim_monge[["primary_key", "delta_g"]],
                      on="primary_key", how="inner")
        return df.merge(
            sim_shibuya_df[["primary_key", "RMSDh"]],
            on="primary_key",
            suffixes=["_heuristic", "_exact"],
        )

    sim = _merge(sim_sh_df)
    sim["absolute_error"] = (sim["RMSDh_heuristic"] - sim["RMSDh_exact"]).abs()
    # theoretical error bound: √(2·Δg·log(n)·k / n)
    sim["theoretical_error"] = np.sqrt(
        2 * sim["delta_g"] * np.log(sim["Residue length"]) * sim["k"]
        / sim["Residue length"]
    )

    # ======================
    # 2. DynDom 2024 dataset
    # ======================
    """
    dyndom_monge = read_dyndom_data("delta_g_dyndom.csv")

    sh_dfs      = [read_dyndom_data(f"rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_{k}.csv", k=k) for k in range(2,6)]
    shibuya_dfs = [read_dyndom_data(f"rmsdh_result/rmsdhk_dyndom_data_{k}.csv",           k=k) for k in range(2,6)]

    dyndom = pd.concat(sh_dfs,      ignore_index=True)
    dyndom_exact = pd.concat(shibuya_dfs, ignore_index=True)
    dyndom["absolute_error"] = (dyndom["RMSDh"] - dyndom_exact["RMSDh"]).abs()
    # attach global Δg and compute theoretical bound
    dyndom["delta_g"] = dyndom_monge["delta_g"]
    dyndom["theoretical_error"] = np.sqrt(
        2 * dyndom["delta_g"] * np.log(dyndom["Residue length"]) * dyndom["k"]
        / dyndom["Residue length"]
    )
    """
        # ------------------------------------------------------------
    # 0. kごとの統計量を計算してLaTeXテーブルを出力
    # ------------------------------------------------------------

    # ratio = absolute_error / theoretical_error
    sim['ratio'] = sim['absolute_error'] / sim['theoretical_error']
    # difference = theoretical_error - absolute_error
    sim['difference'] = sim['theoretical_error'] - sim['absolute_error']

    # kごとに min, average, max を集計
    grouped = sim.groupby('k').agg(
        ratio_min      = ('ratio',      'min'),
        ratio_average  = ('ratio',      'mean'),
        ratio_max      = ('ratio',      'max'),
        diff_min       = ('difference', 'min'),
        diff_average   = ('difference', 'mean'),
        diff_max       = ('difference', 'max'),
    ).rename(index={2:'k=2', 3:'k=3', 4:'k=4', 5:'k=5'})

    # 【1】 全指標をまとめたテーブル
    print(grouped.to_latex(
        float_format="%.3f",
        caption=("Error statistics per $k$ — "
                 "$\\frac{\\mathrm{absolute\\_error}}{\\mathrm{theoretical\\_error}}$ and "
                 "$\\mathrm{theoretical\\_error}-\\mathrm{absolute\\_error}$"),
        label="tab:error_stats_per_k"
    ))

    # 【2】 個別に分けたテーブル（ratio）
    ratio_table = grouped[['ratio_min','ratio_average','ratio_max']].rename(
        columns={'ratio_min':'min','ratio_average':'average','ratio_max':'max'}
    )
    print(ratio_table.to_latex(
        float_format="%.3f",
        caption="Statistics of $\\mathrm{absolute\\_error} / \\mathrm{theoretical\\_error}$ per $k$",
        label="tab:ratio_stats_per_k"
    ))

    # 【3】 個別に分けたテーブル（difference）
    diff_table = grouped[['diff_min','diff_average','diff_max']].rename(
        columns={'diff_min':'min','diff_average':'average','diff_max':'max'}
    )
    print(diff_table.to_latex(
        float_format="%.3f",
        caption="Statistics of $\\mathrm{theoretical\\_error} - \\mathrm{absolute\\_error}$ per $k$",
        label="tab:diff_stats_per_k"
    ))
    # ------------------------------------------------------------
    # 0.1 diff（theoretical_error - absolute_error）のヒストグラムをkごとにプロット
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # 0.1 diff（theoretical_error - absolute_error）のヒストグラムを一行4列でプロット
    # ------------------------------------------------------------
    ks = [2, 3, 4, 5]
    fig_hist, axes_hist = plt.subplots(1, 4, figsize=(24, 6), sharex=True, sharey=True)
    for ax, k in zip(axes_hist, ks):
        data = sim.loc[sim['k'] == k, 'difference']
        ax.hist(data, bins=30, edgecolor='black', alpha=0.7)
        ax.set_title(f"$k = {k}$")
        if k == 2:
            ax.set_ylabel("Count")
        ax.grid(True)
    plt.tight_layout()
    hist_png = os.path.join("figures", "diff_histograms_k2-5_row.png")
    hist_svg = os.path.join("figures", "diff_histograms_k2-5_row.svg")
    fig_hist.savefig(hist_png, dpi=300)
    fig_hist.savefig(hist_svg)
    plt.close(fig_hist)
    print(f"Saved diff histograms for k=2,3,4,5 (1×4): {hist_png}, {hist_svg}")



    # ======================
    # 3. 2×4 subplot grid
    # ======================

    ks = [2, 3, 4, 5]
    fig, axes = plt.subplots(2, 4, figsize=(24, 12), sharex=True, sharey=True)

    for j, k in enumerate(ks):
        # DynDom (row 0)
        """
        ax = axes[0, j]
        df_dyn = dyndom[dyndom["k"] == k]
        ax.scatter(df_dyn["theoretical_error"], df_dyn["absolute_error"], alpha=0.7)
        ax.set_title(f"DynDom, k = {k}")
        ax.grid(True)
        if j == 0:
            ax.set_ylabel("Absolute error")
        """
        # Simulation (row 1)
        ax2 = axes[1, j]
        df_sim = sim[sim["k"] == k]
        ax2.scatter(df_sim["theoretical_error"], df_sim["absolute_error"], alpha=0.7)
        ax2.set_xlim(0, df_sim["theoretical_error"].max())
        ax2.set_ylim(0, df_sim["theoretical_error"].max() / 500)
        ax2.set_title(f"Simulation, k = {k}")
        ax2.grid(True)
        if j == 0:
            ax2.set_ylabel("Absolute error")
        ax2.set_xlabel("Theoretical error bound")

    plt.tight_layout()
    png_path = os.path.join("figures", "SH_error_scatter.png")
    svg_path = os.path.join("figures", "SH_error_scatter.svg")
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path)
    plt.close(fig)
    print(f"Saved 2×4 scatter grid: {png_path}, {svg_path}")

if __name__ == "__main__":
    main()
