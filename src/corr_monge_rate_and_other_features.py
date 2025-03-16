import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18


def safe_corr(df, col1, col2):
    """dfの2列間の相関係数を小数第3位で返す（計算不能の場合は"-"を返す）"""
    try:
        corr_val = df[[col1, col2]].corr().loc[col1, col2]
        if pd.isna(corr_val):
            return "-"
        else:
            return f"{corr_val:.3f}"
    except Exception:
        return "-"


def read_simulation_data(file_path):
    df = pd.read_csv(file_path).fillna("")
    # 実際のヒンジ数がない場合は k 列からコピー
    if "actual_hinge_cnt" not in df.columns:
        df["actual_hinge_cnt"] = df["k"]
    df["primary_key"] = (
        df["p_pdb_id"].apply(lambda x: str(x)[:9])
        + "_hinge_"
        + df["actual_hinge_cnt"].astype(str)
        + "_sigma0.5"
    )
    return df


def read_dyndom_data(file_path, k=None):
    """
    dyndomデータ用のCSVを読み込み、p_pdb_id と q_pdb_id を結合して primary_key を生成する。
    また、k の情報がないため、引数で与えた k を付与する。
    """
    df = pd.read_csv(file_path).fillna("")
    df["primary_key"] = df["p_pdb_id"] + "_" + df["q_pdb_id"]
    if k is not None:
        df["k"] = k
    return df


def main():
    os.makedirs("figures", exist_ok=True)

    # -------------------------------
    # シミュレーションデータの処理
    # -------------------------------
    simulation_shibuya_df = read_simulation_data(
        "rmsdh_result/simulation_shibuya_combined.csv"
    )
    simulation_sh_ilo_df = read_simulation_data(
        "rmsdh_result/simulation_sh_ilo_combined.csv"
    )
    simulation_sh_lo_df = read_simulation_data(
        "rmsdh_result/simulation_sh_lo_combined.csv"
    )
    simulation_sh_df = read_simulation_data("rmsdh_result/simulation_sh_combined.csv")

    # simulation_monge は k=2〜5 のデータを結合
    simulation_monge = pd.DataFrame()
    for k in range(2, 6):
        df = read_simulation_data(f"monge_rate_simulation_{k}_sigma0.5.csv")
        simulation_monge = pd.concat([simulation_monge, df])

    # 各手法について、Monge のデータと、SHibuya の exact 値 (RMSDh) を結合
    simulation_sh_df = simulation_sh_df.merge(
        simulation_monge[["primary_key", "monge_rate_1"]], on="primary_key", how="inner"
    )
    simulation_sh_df = simulation_sh_df.merge(
        simulation_shibuya_df[["primary_key", "RMSDh"]],
        on="primary_key",
        suffixes=["_heuristic", "_exact"],
    )

    simulation_sh_lo_df = simulation_sh_lo_df.merge(
        simulation_monge[["primary_key", "monge_rate_1"]], on="primary_key", how="inner"
    )
    simulation_sh_lo_df = simulation_sh_lo_df.merge(
        simulation_shibuya_df[["primary_key", "RMSDh"]],
        on="primary_key",
        suffixes=["_heuristic", "_exact"],
    )

    simulation_sh_ilo_df = simulation_sh_ilo_df.merge(
        simulation_monge[["primary_key", "monge_rate_1"]], on="primary_key", how="inner"
    )
    simulation_sh_ilo_df = simulation_sh_ilo_df.merge(
        simulation_shibuya_df[["primary_key", "RMSDh"]],
        on="primary_key",
        suffixes=["_heuristic", "_exact"],
    )

    # approximation_ratio の計算（simulation dataset）
    for df in [simulation_sh_df, simulation_sh_lo_df, simulation_sh_ilo_df]:
        df["approximation_ratio"] = df["RMSDh_heuristic"] / df["RMSDh_exact"]

    # simulation dataset の各 k ごとの相関結果
    sim_results = []
    for k in [2, 3, 4, 5]:
        df_sh = simulation_sh_df[simulation_sh_df["actual_hinge_cnt"] == k]
        df_shlo = simulation_sh_lo_df[simulation_sh_lo_df["actual_hinge_cnt"] == k]
        df_shilo = simulation_sh_ilo_df[simulation_sh_ilo_df["actual_hinge_cnt"] == k]

        res_sh = {
            "k": k,
            "Method": "SH",
            "AR": safe_corr(df_sh, "approximation_ratio", "monge_rate_1"),
            "CT": safe_corr(df_sh, "exec_time (s)", "monge_rate_1"),
            "#iterations": "-",
        }
        res_shlo = {
            "k": k,
            "Method": "SH + LO",
            "AR": safe_corr(df_shlo, "approximation_ratio", "monge_rate_1"),
            "CT": safe_corr(df_shlo, "exec_time (s)", "monge_rate_1"),
            "#iterations": "-",
        }
        res_shilo = {
            "k": k,
            "Method": "SH + ILO",
            "AR": safe_corr(df_shilo, "approximation_ratio", "monge_rate_1"),
            "CT": safe_corr(df_shilo, "exec_time (s)", "monge_rate_1"),
            "#iterations": safe_corr(df_shilo, "iter_num", "monge_rate_1"),
        }
        sim_results.extend([res_sh, res_shlo, res_shilo])

    # -------------------------------
    # DynDom 2024 dataset の処理
    # -------------------------------
    dyndom_monge = read_dyndom_data("monge_rate_dyndom.csv")
    sh_2_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_2.csv", k=2
    )
    shibuya_2_dyndom = read_dyndom_data("rmsdh_result/rmsdhk_dyndom_data_2.csv", k=2)
    sh_3_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_3.csv", k=3
    )
    shibuya_3_dyndom = read_dyndom_data("rmsdh_result/rmsdhk_dyndom_data_3.csv", k=3)
    sh_4_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_4.csv", k=4
    )
    sh_5_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_5.csv", k=5
    )
    shibuya_4_dyndom = read_dyndom_data("rmsdh_result/rmsdhk_dyndom_data_4.csv", k=4)
    shibuya_5_dyndom = read_dyndom_data("rmsdh_result/rmsdhk_dyndom_data_5.csv", k=5)

    sh_lo_2_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_2_postpro.csv", k=2
    )
    sh_lo_3_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_3_postpro.csv", k=3
    )
    sh_lo_4_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_4_postpro.csv", k=4
    )
    sh_lo_5_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_5_postpro.csv", k=5
    )

    sh_ilo_2_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdhk_dyndom_2_postpro_loop.csv", k=2
    )
    sh_ilo_3_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdhk_dyndom_3_postpro_loop.csv", k=3
    )
    sh_ilo_4_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdhk_dyndom_4_postpro_loop.csv", k=4
    )
    sh_ilo_5_dyndom = read_dyndom_data(
        "rmsdh_result/fast_rmsdhk_dyndom_5_postpro_loop.csv", k=5
    )

    dyndom_sh = pd.concat([sh_2_dyndom, sh_3_dyndom, sh_4_dyndom, sh_5_dyndom])
    dyndom_shibuya = pd.concat(
        [shibuya_2_dyndom, shibuya_3_dyndom, shibuya_4_dyndom, shibuya_5_dyndom]
    )
    dyndom_sh["approximation_ratio"] = dyndom_sh["RMSDh"] / dyndom_shibuya["RMSDh"]

    dyndom_sh_lo = pd.concat(
        [sh_lo_2_dyndom, sh_lo_3_dyndom, sh_lo_4_dyndom, sh_lo_5_dyndom]
    )
    dyndom_sh_lo["approximation_ratio"] = (
        dyndom_sh_lo["RMSDh"] / dyndom_shibuya["RMSDh"]
    )

    dyndom_sh_ilo = pd.concat(
        [sh_ilo_2_dyndom, sh_ilo_3_dyndom, sh_ilo_4_dyndom, sh_ilo_5_dyndom]
    )
    dyndom_sh_ilo["approximation_ratio"] = (
        dyndom_sh_ilo["RMSDh"] / dyndom_shibuya["RMSDh"]
    )
    dyn_results = []
    for k in [2, 3, 4, 5]:
        df_sh = dyndom_sh[dyndom_sh["k"] == k]
        df_shlo = dyndom_sh_lo[dyndom_sh_lo["k"] == k]
        df_shilo = dyndom_sh_ilo[dyndom_sh_ilo["k"] == k]
        df_sh["monge_rate_1"] = dyndom_monge["monge_rate_1"]
        df_shlo["monge_rate_1"] = dyndom_monge["monge_rate_1"]
        df_shilo["monge_rate_1"] = dyndom_monge["monge_rate_1"]

        res_sh = {
            "k": k,
            "Method": "SH",
            "AR": safe_corr(df_sh, "approximation_ratio", "monge_rate_1"),
            "CT": safe_corr(df_sh, "exec_time (s)", "monge_rate_1"),
            "#iterations": "-",
        }
        res_shlo = {
            "k": k,
            "Method": "SH + LO",
            "AR": safe_corr(df_shlo, "approximation_ratio", "monge_rate_1"),
            "CT": safe_corr(df_shlo, "exec_time (s)", "monge_rate_1"),
            "#iterations": "-",
        }
        res_shilo = {
            "k": k,
            "Method": "SH + ILO",
            "AR": safe_corr(df_shilo, "approximation_ratio", "monge_rate_1"),
            "CT": safe_corr(df_shilo, "exec_time (s)", "monge_rate_1"),
            "#iterations": safe_corr(df_shilo, "iter_num", "monge_rate_1"),
        }
        dyn_results.extend([res_sh, res_shlo, res_shilo])

    # -------------------------------
    # シミュレーションと DynDom の結果を左右に並べた LaTeX テーブル作成
    # -------------------------------
    header = r"""\begin{table*}[b]
    \centering
    \caption{Correlation coefficients between the Monge rate and AR (approximation ratio), CT (computation time), and \#iterations for SH (SMAWK-based heuristic DP for $RMSDh^{(k)}$), SH + LO (local optimization postprocessing for $RMSDh^{(k)}$), and SH + ILO (iterative LO), computed on both the simulation dataset and DynDom 2024 dataset.}\label{tab:corr_monge}
    \begin{tabular}{ccrrrrrr}
    \hline
    \multirow{2}{*}{$k$} & \multirow{2}{*}{Method} & \multicolumn{3}{c}{Simulation dataset} & \multicolumn{3}{c}{DynDom 2024 dataset} \\
    \cline{3-8}
     & & AR & CT & \#iterations & AR & CT & \#iterations \\
    \hline
    """
    footer = r"""\hline
    \end{tabular}
    \vspace{1ex}
\end{table*}
    """
    table_rows = ""
    # 各 k ごとに 3 行ずつ出力（順番は SH, SH+LO, SH+ILO）
    for k in [2, 3, 4, 5]:
        # 抽出
        sim_group = [r for r in sim_results if r["k"] == k]
        dyn_group = [r for r in dyn_results if r["k"] == k]
        first_row = (
            r"\multirow{3}{*}{$"
            + f"{k}"
            + r"$} & "
            + f"{sim_group[0]['Method']} & {sim_group[0]['AR']} & {sim_group[0]['CT']} & {sim_group[0]['#iterations']} & "
            + f"{dyn_group[0]['AR']} & {dyn_group[0]['CT']} & {dyn_group[0]['#iterations']} \\\\"
            + "\n"
        )
        second_row = (
            " & "
            + f"{sim_group[1]['Method']} & {sim_group[1]['AR']} & {sim_group[1]['CT']} & {sim_group[1]['#iterations']} & "
            + f"{dyn_group[1]['AR']} & {dyn_group[1]['CT']} & {dyn_group[1]['#iterations']} \\\\"
            + "\n"
        )
        third_row = (
            " & "
            + f"{sim_group[2]['Method']} & {sim_group[2]['AR']} & {sim_group[2]['CT']} & {sim_group[2]['#iterations']} & "
            + f"{dyn_group[2]['AR']} & {dyn_group[2]['CT']} & {dyn_group[2]['#iterations']} \\\\"
            + "\n"
        )
        table_rows += first_row + second_row + third_row + r"\hline" + "\n"

    tex_str = header + table_rows + footer

    # LaTeX テーブルをファイルに出力
    tex_filepath = os.path.join("figures", "corr_monge.tex")
    with open(tex_filepath, "w") as f:
        f.write(tex_str)
    print(f"LaTeX file saved as {tex_filepath}")


if __name__ == "__main__":
    main()
