import pandas as pd
import os


def main():
    os.makedirs("figures", exist_ok=True)
    simulation_sh_ilo_df = read_simulation_data(
        "rmsdh_result/simulation_sh_ilo_combined.csv"
    )
    simulation_sh_lo_df = read_simulation_data(
        "rmsdh_result/simulation_sh_lo_combined.csv"
    )
    simulation_sh_df = read_simulation_data("rmsdh_result/simulation_sh_combined.csv")
    simulation_r_ilo_df = read_simulation_data(
        "rmsdh_result/simulation_r_ilo_combined.csv"
    )
    simulation_r_lo_df = read_simulation_data(
        "rmsdh_result/simulation_r_lo_combined.csv"
    )
    simulation_shibuya_df = read_simulation_data(
        "rmsdh_result/simulation_shibuya_combined.csv"
    )
    simulation_dyndom_df = read_fatcat_dyndom_data(
        "dyndom_simulation_hinge_count_result.csv",
        "dyndom_simulation_result/dyndom_simulation_execution_time_improved.csv",
        simulation_shibuya_df,
    )
    simulation_fatcat_df = read_fatcat_dyndom_data(
        "fatcat_simulation_hinge_count_result.csv",
        "fatcat_simulation_execution_time_improved.csv",
        simulation_shibuya_df,
    )
    assert (
        len(simulation_dyndom_df)
        == len(simulation_dyndom_df)
        == len(simulation_r_ilo_df)
        == len(simulation_sh_df)
        == len(simulation_r_lo_df)
        == len(simulation_sh_lo_df)
        == len(simulation_sh_ilo_df)
    )
    df_dict = {}
    df_dict["R + LO"] = simulation_r_lo_df
    df_dict["R + ILO"] = simulation_r_ilo_df
    df_dict["SH"] = simulation_sh_df
    df_dict["SH + LO"] = simulation_sh_lo_df
    df_dict["SH + ILO"] = simulation_sh_ilo_df
    df_dict["Shibuya's method"] = simulation_shibuya_df
    df_dict["FATCAT"] = simulation_fatcat_df
    df_dict["DynDom"] = simulation_dyndom_df

    ks = [2, 3, 4, 5]
    results = []

    for method, df in df_dict.items():
        df_filtered = df[df["actual_hinge_cnt"].isin(ks)]
        for k in ks:
            if method in ["R + LO", "R + ILO"]:
                col_name = f"{k}_computation_time"
                if col_name in df_filtered.columns:
                    sub_df = df_filtered[df_filtered["actual_hinge_cnt"] == k]
                    total_time = sub_df[col_name].sum() / 1000  # ms -> s
                    avg_time = sub_df[col_name].mean() / 1000
                else:
                    total_time = None
                    avg_time = None
            else:
                if "exec_time (s)" in df_filtered.columns:
                    sub_df = df_filtered[df_filtered["actual_hinge_cnt"] == k]
                    total_time = sub_df["exec_time (s)"].sum()
                    avg_time = sub_df["exec_time (s)"].mean()
                else:
                    total_time = None
                    avg_time = None
            results.append((k, method, total_time, avg_time))

    # 結果を k ごとにグループ化
    results_by_k = {k: [] for k in ks}
    for k, method, total, avg in results:
        results_by_k[k].append((method, total, avg))

    # LaTeX の表形式コードを生成（multirowを使用）
    latex_lines = []
    latex_lines.append(r"\begin{table*}[t]")
    latex_lines.append(r"\centering")
    latex_lines.append(r"\begin{tabular}{llrr}")
    latex_lines.append(r"\hline")
    latex_lines.append(r"$k$ & Method & Total (s) & Average (s) \\")
    latex_lines.append(r"\hline")

    for k in ks:
        method_results = results_by_k[k]
        n_methods = len(method_results)
        for i, (method, total, avg) in enumerate(method_results):
            total_str = f"{total:.3f}" if total is not None else "--"
            avg_str = f"{avg:.3f}" if avg is not None else "--"
            if i == 0:
                # multirow: k値をまとめて表示
                latex_lines.append(
                    r"\multirow{"
                    + f"{n_methods}"
                    + r"}{*}{"
                    + f"{k}"
                    + r"} & "
                    + f"{method} & {total_str} & {avg_str} \\\\"
                )
            else:
                latex_lines.append(" & " + f"{method} & {total_str} & {avg_str} \\\\")
        latex_lines.append(r"\hline")

    latex_lines.append(r"\end{tabular}")
    latex_lines.append(
        r"\caption{Simulation Computation Times: Total and Average (s) for each method grouped by $k$.}"
    )
    latex_lines.append(r"\label{tab:simulation_comp_time}")
    latex_lines.append(r"\end{table*}")

    tex_file = os.path.join("figures", "simulation_comp_time_table.tex")
    with open(tex_file, "w", encoding="utf-8") as f:
        f.write("\n".join(latex_lines))

    print(f"LaTeXコードを {tex_file} に出力しました。")


def read_simulation_data(file_path):
    df = pd.read_csv(file_path).fillna("")
    df["actual_hinge_cnt"] = df["k"]
    df["primary_key"] = (
        df["p_pdb_id"].apply(lambda x: str(x)[:9])
        + "_hinge_"
        + df["actual_hinge_cnt"].astype(str)
        + "_sigma0.5"
    )
    return df


def read_fatcat_dyndom_data(file_path1, file_path2, simulation_shibuya_df):
    df1 = pd.read_csv(file_path1).fillna("")
    df1["primary_key"] = (
        df1["q_pdb"].apply(lambda x: str(x).split("/")[1]).str.replace(".pdb", "")
    )
    df1 = df1.rename(columns={"detected_hinge_count": "hinge_cnt"})
    df2 = pd.read_csv(file_path2)
    if "fatcat" in file_path2:
        df2["primary_key"] = df2["pdb"].str.replace(".pdb", "")
    else:
        df2["primary_key"] = df2["pdb_pair"].apply(lambda x: x[19:-8])
    df1 = df1.merge(df2, on=["primary_key"]).rename(
        columns={"execution_time": "exec_time (s)"}
    )
    df1 = df1.merge(
        simulation_shibuya_df[["primary_key", "actual_hinge_cnt"]], on=["primary_key"]
    )
    return df1


if __name__ == "__main__":
    main()
