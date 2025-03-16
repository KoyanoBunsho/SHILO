import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18


def main():
    os.makedirs("figures", exist_ok=True)
    ok_pdb = pd.read_csv("../notebooks/ok_pdb.csv")
    ok_pdb["p_pdb_id"] = ok_pdb["p_pdb_id"].str.lower()
    ok_pdb["q_pdb_id"] = ok_pdb["q_pdb_id"].str.lower()
    ok_keys = set(
        zip(
            ok_pdb["p_pdb_id"],
            ok_pdb["q_pdb_id"],
        )
    )

    # 解析対象の k 値
    k_values = [2, 3, 4, 5]

    # ---------------------------
    # SH+ILO Combined のデータ
    # ---------------------------
    simulation_sh_ilo_df = read_simulation_data(
        "rmsdh_result/simulation_sh_ilo_combined.csv"
    )
    shilo_avgs = []
    shilo_maxs = []
    for k in k_values:
        shilo_subset = simulation_sh_ilo_df[simulation_sh_ilo_df["k"] == k]
        s_avg, s_max = compute_metrics(shilo_subset)
        shilo_avgs.append(s_avg)
        shilo_maxs.append(s_max)
        print(f"SH+ILO Combined: k = {k}: Average: {s_avg}, Maximum: {s_max}")

    # ---------------------------
    # Shibuya のデータ
    # ---------------------------
    shibuya_dfs = {
        2: pd.read_csv("rmsdh_result/fast_rmsdh_hingek_cnt_2_postpro_loop.csv"),
        3: pd.read_csv("rmsdh_result/fast_rmsdh_hingek_cnt_3_postpro_loop.csv"),
        4: pd.read_csv("rmsdh_result/fast_rmsdh_hingek_cnt_4_postpro_loop.csv"),
        5: pd.read_csv("rmsdh_result/fast_rmsdh_hingek_cnt_5_postpro_loop.csv"),
    }
    shibuya_avgs = []
    shibuya_maxs = []
    for k in k_values:
        df = shibuya_dfs[k]
        df["actual_hinge_cnt"] = k
        avg_val, max_val = compute_metrics(df)
        shibuya_avgs.append(avg_val)
        shibuya_maxs.append(max_val)
        print(f"Shibuya: k = {k}: Average: {avg_val}, Maximum: {max_val}")

    # ---------------------------
    # PAR のデータ
    # ---------------------------
    par_dfs = {
        2: filter_by_ok_keys(pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_2_pospro_loop.csv"), ok_keys),
        3: filter_by_ok_keys(pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_3_pospro_loop.csv"), ok_keys),
        4: filter_by_ok_keys(pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_4_pospro_loop.csv"), ok_keys),
        5: filter_by_ok_keys(pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_5_pospro_loop.csv"), ok_keys),
    }
    par_avgs = []
    par_maxs = []
    for k in k_values:
        df = par_dfs[k]
        df["actual_hinge_cnt"] = k
        avg_val, max_val = compute_metrics(df)
        par_avgs.append(avg_val)
        par_maxs.append(max_val)
        print(f"PAR: k = {k}: Average: {avg_val}, Maximum: {max_val}")

    # ---------------------------
    # Dyndom のデータ
    # ---------------------------
    dyndom_dfs = {
        2: pd.read_csv("rmsdh_result/fast_rmsdhk_dyndom_2_postpro_loop.csv"),
        3: pd.read_csv("rmsdh_result/fast_rmsdhk_dyndom_3_postpro_loop.csv"),
        4: pd.read_csv("rmsdh_result/fast_rmsdhk_dyndom_4_postpro_loop.csv"),
        5: pd.read_csv("rmsdh_result/fast_rmsdhk_dyndom_5_postpro_loop.csv"),
    }
    dyndom_avgs = []
    dyndom_maxs = []
    for k in k_values:
        df = dyndom_dfs[k]
        df["actual_hinge_cnt"] = k
        avg_val, max_val = compute_metrics(df)
        dyndom_avgs.append(avg_val)
        dyndom_maxs.append(max_val)
        print(f"Dyndom: k = {k}: Average: {avg_val}, Maximum: {max_val}")

    # ---------------------------
    # 1行4列のサブプロットでプロット
    # ---------------------------
    fig, axes = plt.subplots(1, 4, figsize=(20, 6), sharey=True)

    # 各サブプロットごとにタイトルとデータを設定
    datasets = [
        ("Simulation", shilo_avgs, shilo_maxs),
        ("Shibuya", shibuya_avgs, shibuya_maxs),
        ("PAR", par_avgs, par_maxs),
        ("DynDom", dyndom_avgs, dyndom_maxs),
    ]
    caption_dict = {"Simulation": "(a)","Shibuya": "(b)","PAR": "(c)","DynDom": "(d)"}
    for ax, (title, avg_list, max_list) in zip(axes, datasets):
        # 平均値を丸マーカー、最大値を四角マーカーでプロット
        ax.plot(k_values, avg_list, marker="o", linestyle="-", label="Average")
        ax.plot(k_values, max_list, marker="s", linestyle="-", label="Maximum")
        ax.set_xlabel(f"{caption_dict[title]}", fontsize=24)
        ax.set_xticks(k_values)
        ax.grid(True)
        ax.legend()

    axes[0].set_ylabel("#Iterations", fontsize=24)

    plt.tight_layout()
    plt.savefig("figures/all_sh_ilo_iteration_metrics.svg", format="svg")
    plt.savefig("figures/all_sh_ilo_iteration_metrics.png")
    plt.close()


def compute_metrics(df):
    # "iter_num" 列を数値に変換し、変換できない値は 0 に置換
    iteration_values = pd.to_numeric(df["iter_num"], errors="coerce").fillna(0)
    avg_val = iteration_values.mean()
    max_val = iteration_values.max()
    return avg_val, max_val


def read_simulation_data(file_path):
    df = pd.read_csv(file_path)
    # k 列を実際の hinge 数として利用
    df["actual_hinge_cnt"] = df["k"]
    return df


def filter_by_ok_keys(df, ok_keys):
    return df[
        df.apply(
            lambda row: (row["p_pdb_id"].lower(), row["q_pdb_id"].lower()) in ok_keys,
            axis=1,
        )
    ]


if __name__ == "__main__":
    main()
