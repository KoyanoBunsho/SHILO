#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the correlation between Monge rate and
    • AE  (absolute error = |RMSDhheuristic – RMSDhexact|)
    • CT  (computation time)
    • #iterations
for SH, SH + LO, SH + ILO on both the simulation dataset and the DynDom 2024 dataset,
and output the results as a LaTeX table.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# matplotlib setup
# ------------------------------------------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18


# ------------------------------------------------------------
# utility functions
# ------------------------------------------------------------
def safe_corr(df: pd.DataFrame, col1: str, col2: str) -> str:
    """
    Return the Pearson correlation coefficient between df[col1] and df[col2],
    rounded to the 3rd decimal place.  If it cannot be computed, return "-".
    """
    try:
        r = df[[col1, col2]].corr().loc[col1, col2]
        return "-" if pd.isna(r) else f"{r:.3f}"
    except Exception:
        return "-"


def read_simulation_data(file_path: str) -> pd.DataFrame:
    """
    Read a simulation CSV, fill NA with "", create primary_key, and
    (if missing) copy `k` to `actual_hinge_cnt`.
    """
    df = pd.read_csv(file_path).fillna("")
    if "actual_hinge_cnt" not in df.columns:
        df["actual_hinge_cnt"] = df["k"]
    df["primary_key"] = (
        df["p_pdb_id"].apply(lambda x: str(x)[:9])
        + "_hinge_"
        + df["actual_hinge_cnt"].astype(str)
        + "_sigma0.5"
    )
    return df


def read_dyndom_data(file_path: str, k: int | None = None) -> pd.DataFrame:
    """
    Read a DynDom CSV, create primary_key from p_pdb_id + q_pdb_id,
    and (optionally) attach a given k value.
    """
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
    sim_sh_ilo_df   = read_simulation_data("rmsdh_result/simulation_sh_ilo_combined.csv")
    sim_sh_lo_df    = read_simulation_data("rmsdh_result/simulation_sh_lo_combined.csv")
    sim_sh_df       = read_simulation_data("rmsdh_result/simulation_sh_combined.csv")

    # Monge rates for k = 2–5 concatenated
    sim_monge = pd.concat(
        [read_simulation_data(f"delta_g_simulation_{k}_sigma0.5.csv") for k in range(2, 6)]
    )

    # --- merge Monge rates & exact RMSDh (Shibuya) ---
    def _merge(df: pd.DataFrame) -> pd.DataFrame:
        df = df.merge(sim_monge[["primary_key", "delta_g"]], on="primary_key", how="inner")
        return df.merge(
            sim_shibuya_df[["primary_key", "RMSDh"]],
            on="primary_key",
            suffixes=["_heuristic", "_exact"],
        )

    sim_sh_df    = _merge(sim_sh_df)
    sim_sh_lo_df = _merge(sim_sh_lo_df)
    sim_sh_ilo_df = _merge(sim_sh_ilo_df)

    # --- absolute error ---
    for df in (sim_sh_df, sim_sh_lo_df, sim_sh_ilo_df):
        df["absolute_error"] = (df["RMSDh_heuristic"] - df["RMSDh_exact"]).abs()

    # --- correlation results per-k ---
    sim_results: list[dict] = []
    for k in (2, 3, 4, 5):
        df_sh    = sim_sh_df[   sim_sh_df["actual_hinge_cnt"] == k]
        df_shlo  = sim_sh_lo_df[sim_sh_lo_df["actual_hinge_cnt"] == k]
        df_shilo = sim_sh_ilo_df[sim_sh_ilo_df["actual_hinge_cnt"] == k]

        sim_results.extend(
            [
                {
                    "k": k,
                    "Method": "SH",
                    "AE": safe_corr(df_sh,    "absolute_error", "delta_g"),
                    "CT": safe_corr(df_sh,    "exec_time (s)",  "delta_g"),
                    "#iterations": "-",
                },
                {
                    "k": k,
                    "Method": "SH + LO",
                    "AE": safe_corr(df_shlo,  "absolute_error", "delta_g"),
                    "CT": safe_corr(df_shlo,  "exec_time (s)",  "delta_g"),
                    "#iterations": "-",
                },
                {
                    "k": k,
                    "Method": "SH + ILO",
                    "AE": safe_corr(df_shilo, "absolute_error", "delta_g"),
                    "CT": safe_corr(df_shilo, "exec_time (s)",  "delta_g"),
                    "#iterations": safe_corr(df_shilo, "iter_num", "delta_g"),
                },
            ]
        )

    # ======================
    # 2. DynDom 2024 dataset
    # ======================
    dyndom_monge   = read_dyndom_data("delta_g_dyndom.csv")

    # SH & Shibuya (exact) results for k=2–5
    sh_dfs = [
        read_dyndom_data(f"rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_{k}.csv", k=k)
        for k in range(2, 6)
    ]
    shibuya_dfs = [
        read_dyndom_data(f"rmsdh_result/rmsdhk_dyndom_data_{k}.csv", k=k)
        for k in range(2, 6)
    ]
    sh_lo_dfs = [
        read_dyndom_data(
            f"rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_{k}_postpro.csv", k=k
        )
        for k in range(2, 6)
    ]
    sh_ilo_dfs = [
        read_dyndom_data(
            f"rmsdh_result/fast_rmsdhk_dyndom_{k}_postpro_loop.csv", k=k
        )
        for k in range(2, 6)
    ]

    dyndom_sh       = pd.concat(sh_dfs,       ignore_index=True)
    dyndom_shibuya  = pd.concat(shibuya_dfs,  ignore_index=True)
    dyndom_sh_lo    = pd.concat(sh_lo_dfs,    ignore_index=True)
    dyndom_sh_ilo   = pd.concat(sh_ilo_dfs,   ignore_index=True)

    # --- absolute error ---
    dyndom_sh["absolute_error"]     = (dyndom_sh["RMSDh"]     - dyndom_shibuya["RMSDh"]).abs()
    dyndom_sh_lo["absolute_error"]  = (dyndom_sh_lo["RMSDh"]  - dyndom_shibuya["RMSDh"]).abs()
    dyndom_sh_ilo["absolute_error"] = (dyndom_sh_ilo["RMSDh"] - dyndom_shibuya["RMSDh"]).abs()

    # --- attach delta_g column (common for all primary keys) ---
    for df in (dyndom_sh, dyndom_sh_lo, dyndom_sh_ilo):
        df["delta_g"] = dyndom_monge["delta_g"]

    # --- correlation results per-k ---
    dyn_results: list[dict] = []
    for k in (2, 3, 4, 5):
        df_sh    = dyndom_sh[    dyndom_sh["k"] == k]
        df_shlo  = dyndom_sh_lo[ dyndom_sh_lo["k"] == k]
        df_shilo = dyndom_sh_ilo[dyndom_sh_ilo["k"] == k]

        dyn_results.extend(
            [
                {
                    "k": k,
                    "Method": "SH",
                    "AE": safe_corr(df_sh,    "absolute_error", "delta_g"),
                    "CT": safe_corr(df_sh,    "exec_time (s)",  "delta_g"),
                    "#iterations": "-",
                },
                {
                    "k": k,
                    "Method": "SH + LO",
                    "AE": safe_corr(df_shlo,  "absolute_error", "delta_g"),
                    "CT": safe_corr(df_shlo,  "exec_time (s)",  "delta_g"),
                    "#iterations": "-",
                },
                {
                    "k": k,
                    "Method": "SH + ILO",
                    "AE": safe_corr(df_shilo, "absolute_error", "delta_g"),
                    "CT": safe_corr(df_shilo, "exec_time (s)",  "delta_g"),
                    "#iterations": safe_corr(df_shilo, "iter_num", "delta_g"),
                },
            ]
        )

    # ======================
    # 3. LaTeX table
    # ======================
    header = r"""\begin{table*}[b]
    \centering
    \caption{Correlation coefficients between the Monge rate and AE (absolute error), CT (computation time), and \#iterations for SH (SMAWK-based heuristic DP for $RMSDh^{(k)}$), SH + LO (local optimisation post-processing for $RMSDh^{(k)}$), and SH + ILO (iterative LO), computed on both the simulation dataset and DynDom 2024 dataset.}\label{tab:corr_monge}
    \begin{tabular}{ccrrrrrr}
    \hline
    \multirow{2}{*}{$k$} & \multirow{2}{*}{Method} & \multicolumn{3}{c}{Simulation dataset} & \multicolumn{3}{c}{DynDom 2024 dataset} \\
    \cline{3-8}
     & & AE & CT & \#iterations & AE & CT & \#iterations \\
    \hline
"""
    footer = r"""\hline
    \end{tabular}
    \vspace{1ex}
\end{table*}
"""
    # --- assemble rows (three per k) ---
    body = ""
    for k in (2, 3, 4, 5):
        sim_group = [r for r in sim_results if r["k"] == k]
        dyn_group = [r for r in dyn_results if r["k"] == k]

        body += (
            r"\multirow{3}{*}{$" + str(k) + r"$} & "
            + f"{sim_group[0]['Method']} & {sim_group[0]['AE']} & {sim_group[0]['CT']} & {sim_group[0]['#iterations']} & "
            + f"{dyn_group[0]['AE']} & {dyn_group[0]['CT']} & {dyn_group[0]['#iterations']} \\\n"
            " & "
            + f"{sim_group[1]['Method']} & {sim_group[1]['AE']} & {sim_group[1]['CT']} & {sim_group[1]['#iterations']} & "
            + f"{dyn_group[1]['AE']} & {dyn_group[1]['CT']} & {dyn_group[1]['#iterations']} \\\n"
            " & "
            + f"{sim_group[2]['Method']} & {sim_group[2]['AE']} & {sim_group[2]['CT']} & {sim_group[2]['#iterations']} & "
            + f"{dyn_group[2]['AE']} & {dyn_group[2]['CT']} & {dyn_group[2]['#iterations']} \\\n"
            r"\hline" + "\n"
        )

    tex_str = header + body + footer

    tex_path = os.path.join("figures", "corr_delta_g.tex")
    with open(tex_path, "w", encoding="utf-8") as fh:
        fh.write(tex_str)
    print(f"LaTeX table saved to {tex_path}")


# ------------------------------------------------------------
if __name__ == "__main__":
    main()
