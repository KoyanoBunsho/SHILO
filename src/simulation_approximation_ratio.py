import pandas as pd
import matplotlib.pyplot as plt
import os
from tqdm import tqdm


def main():
    # 各シミュレーション結果の読み込み
    simulation_sh_ilo_df = read_simulation_data(
        "rmsdh_result/simulation_sh_ilo_combined.csv"
    )
    simulation_sh_lo_df = read_simulation_data(
        "rmsdh_result/simulation_sh_lo_combined.csv"
    )
    simulation_sh_df = read_simulation_data("rmsdh_result/simulation_sh_combined.csv")
    simulation_r_ilo_df = read_simulation_data(
        "rmsdh_result/simulation_r_ilo_combined.csv.gz"
    )
    simulation_r_lo_df = read_simulation_data(
        "rmsdh_result/simulation_r_lo_combined.csv.gz"
    )
    simulation_shibuya_df = read_simulation_data(
        "rmsdh_result/simulation_shibuya_combined.csv"
    )
    df_dict = {
        "R + LO": simulation_r_lo_df,
        "R + ILO": simulation_r_ilo_df,
        "SH": simulation_sh_df,
        "SH + LO": simulation_sh_lo_df,
        "SH + ILO": simulation_sh_ilo_df,
    }
    results = []
    for method, df in df_dict.items():
        for k in [2, 3, 4]:
            df_k = df[df["actual_hinge_cnt"] == k]
            shibuya_k = simulation_shibuya_df[
                simulation_shibuya_df["actual_hinge_cnt"] == k
            ]
            merge_df = pd.merge(
                df_k, shibuya_k, on=["p_pdb_id"], suffixes=["_heuristic", "_exact"]
            )
            if (method == "R + LO") or (method == "R + ILO"):
                rmsdh_cols = [f"{i}_RMSDhk" for i in range(100)]
                ratio_series = merge_df[rmsdh_cols].div(merge_df["RMSDh_exact"], axis=0)
                avg_ratio = ratio_series.mean().mean()
                max_ratio = ratio_series.max().mean()
                min_ratio = ratio_series.min().mean()
            else:
                ratio_series = merge_df["RMSDh_heuristic"] / merge_df["RMSDh_exact"]
                avg_ratio = ratio_series.mean()
                max_ratio = ratio_series.max()
                min_ratio = ratio_series.min()

            results.append(
                {
                    "Method": method,
                    "k": k,
                    "Avg": format_stats(avg_ratio),
                    "Max": format_stats(max_ratio),
                    "Min": format_stats(min_ratio),
                }
            )

    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values(by=["k", "Method"]).reset_index(drop=True)
    latex_table = result_df.to_latex(index=False, escape=False)
    print(latex_table)


def read_simulation_data(file_path):
    df = pd.read_csv(file_path).fillna("")
    """
    df["primary_key"] = (
        df["p_pdb_id"].apply(lambda x: str(x)[:9])
        + "_hinge_"
        + df["k"].astype(str)
        + "_sigma"
        + df["sigma"].astype(str)
    )
    """
    df["actual_hinge_cnt"] = df["k"]
    return df


def format_stats(stats):
    return f"{stats:.4f}"


if __name__ == "__main__":
    main()
