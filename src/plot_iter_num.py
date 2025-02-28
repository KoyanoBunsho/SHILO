import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"


def main():
    os.makedirs("figures", exist_ok=True)
    # TODO: parデータをok_pdbで限定する
    # shilo データの読み込み
    shilo_2_shibuya = pd.read_csv(
        "rmsdh_result/fast_rmsdh_hingek_cnt_2_postpro_loop.csv"
    )
    shilo_3_shibuya = pd.read_csv(
        "rmsdh_result/fast_rmsdh_hingek_cnt_3_postpro_loop.csv"
    )
    shilo_4_shibuya = pd.read_csv(
        "rmsdh_result/fast_rmsdh_hingek_cnt_4_postpro_loop.csv"
    )

    shilo_2_par = pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_2_pospro_loop.csv")
    shilo_3_par = pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_3_pospro_loop.csv")
    shilo_4_par = pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_4_pospro_loop.csv")

    shilo_2_dyndom = pd.read_csv("rmsdh_result/fast_rmsdhk_dyndom_2_postpro_loop.csv")
    shilo_3_dyndom = pd.read_csv("rmsdh_result/fast_rmsdhk_dyndom_3_postpro_loop.csv")
    shilo_4_dyndom = pd.read_csv("rmsdh_result/fast_rmsdhk_dyndom_4_postpro_loop.csv")

    # rilo データの読み込み
    rilo_2_shibuya = pd.read_csv("rmsdh_result/ablation_study_loop2.csv")
    rilo_3_shibuya = pd.read_csv("rmsdh_result/ablation_study_loop3.csv")
    rilo_4_shibuya = pd.read_csv("rmsdh_result/ablation_study_loop4.csv")

    rilo_2_par = pd.read_csv("rmsdh_result/ablation_study_loop_par2.csv")
    rilo_3_par = pd.read_csv("rmsdh_result/ablation_study_loop_par3.csv")
    rilo_4_par = pd.read_csv("rmsdh_result/ablation_study_loop_par4.csv")

    rilo_2_dyndom = pd.read_csv("rmsdh_result/ablation_study_loop_dyndom2.csv")
    rilo_3_dyndom = pd.read_csv("rmsdh_result/ablation_study_loop_dyndom3.csv")
    rilo_4_dyndom = pd.read_csv("rmsdh_result/ablation_study_loop_dyndom4.csv")

    # データセットのリスト化 (k=2,3,4 の各 Shibuya, Par, Dyndom)
    rilo_datasets = [
        [rilo_2_shibuya, rilo_2_par, rilo_2_dyndom],
        [rilo_3_shibuya, rilo_3_par, rilo_3_dyndom],
        [rilo_4_shibuya, rilo_4_par, rilo_4_dyndom],
    ]

    shilo_datasets = [
        [shilo_2_shibuya, shilo_2_par, shilo_2_dyndom],
        [shilo_3_shibuya, shilo_3_par, shilo_3_dyndom],
        [shilo_4_shibuya, shilo_4_par, shilo_4_dyndom],
    ]

    titles = ["k=2", "k=3", "k=4"]
    col_titles = ["Shibuya 2008", "PAR 2020", "DynDom 2024"]

    iter_num_cols_list = [f"{i}_iter_num" for i in range(100)]

    _, axes = plt.subplots(3, 3, figsize=(15, 10))

    for i, (rilo_list, shilo_list, row_title) in enumerate(
        zip(rilo_datasets, shilo_datasets, titles)
    ):
        for j, (rilo_df, shilo_df, col_title) in enumerate(
            zip(rilo_list, shilo_list, col_titles)
        ):
            iter_num_cols = [
                col for col in iter_num_cols_list if col in rilo_df.columns
            ]
            rilo_counts = (
                rilo_df[iter_num_cols]
                .apply(pd.value_counts)
                .fillna(0)
                .sum(axis=1)
                .reset_index()
            )
            rilo_counts.columns = ["Value", "Count"]
            shilo_counts = shilo_df["iter_num"].value_counts().reset_index()
            shilo_counts.columns = ["Value", "Count"]
            total_rilo = rilo_counts["Count"].sum()
            total_shilo = shilo_counts["Count"].sum()

            if total_rilo > 0:
                rilo_counts["Percentage"] = (rilo_counts["Count"] / total_rilo) * 100
            else:
                rilo_counts["Percentage"] = 0

            if total_shilo > 0:
                shilo_counts["Percentage"] = (shilo_counts["Count"] / total_shilo) * 100
            else:
                shilo_counts["Percentage"] = 0

            axes[i, j].bar(
                rilo_counts["Value"],
                rilo_counts["Percentage"],
                label="R + ILO",
                alpha=0.7,
            )
            axes[i, j].bar(
                shilo_counts["Value"],
                shilo_counts["Percentage"],
                label="SH + ILO",
                alpha=0.5,
            )
            axes[i, j].set_title(f"{row_title} - {col_title}")
            axes[i, j].set_xlabel("Value")
            axes[i, j].set_ylabel("Percentage (%)")
            axes[i, j].set_xlim(1, 7)
            axes[i, j].legend()

    plt.tight_layout()
    plt.savefig("figures/iteration_percentage.svg", format="svg")
    plt.savefig("figures/iteration_percentage.png")
    plt.close()


if __name__ == "__main__":
    main()
