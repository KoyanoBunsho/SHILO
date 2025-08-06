import pandas as pd
import matplotlib.pyplot as plt

# フォント設定
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18

def main():
    #─── 静的手法３種の読み込み & プロット ───#
    methods = ["shibuya", "par", "dyndom"]
    _, axes = plt.subplots(1, 4, figsize=(20, 5), sharey=True)

    for i, (ax, method) in enumerate(zip(axes[:3], methods)):
        delta_g_df = pd.read_csv(f"delta_g_{method}.csv")
        avg_delta_df = pd.read_csv(f"average_delta_g_{method}.csv")
        # df に 'delta_g' と 'average_delta_g' 列があることを想定
        df = pd.concat([delta_g_df, avg_delta_df], axis=1)
        ax.scatter(df["delta_g"], df["average_delta_g"], s=1)
        ax.set_xlim(0, df["delta_g"].max())
        ax.set_ylim(0, df["average_delta_g"].max())
        ax.set_xlabel(r"$\Delta_G$"+f"\n({chr(ord('a') + i)})")
        if i == 0:
            ax.set_ylabel(r"$\overline{\Delta_G}$")
        ax.grid(True)

    #─── シミュレーションの読み込み & マージ ───#
    sim_nums = [2, 3, 4, 5]
    sim_dfs = []
    for num in sim_nums:
        df_sim_delta = pd.read_csv(f"delta_g_simulation_{num}_sigma0.5.csv")
        df_sim_avg_delta = pd.read_csv(f"average_delta_g_simulation_{num}_sigma0.5.csv")
        df_sim = df_sim_delta.merge(df_sim_avg_delta, on="p_pdb_id")
        df_sim["sim_id"] = f"sim_{num}"
        sim_dfs.append(df_sim)
    df_sim_all = pd.concat(sim_dfs, ignore_index=True)

    #─── シミュレーションまとめプロット ───#
    ax = axes[3]
    ax.scatter(df_sim_all["delta_g"], df_sim_all["average_delta_g"], s=1)
    ax.set_xlim(0, df_sim_all["delta_g"].max())
    ax.set_ylim(0, df_sim_all["average_delta_g"].max())
    ax.set_xlabel(r"$\Delta_G$"+f"\n(d)")
    ax.grid(True)

    plt.tight_layout()
    plt.savefig("comparison_average_delta_g.png", bbox_inches="tight")
    plt.savefig("comparison_average_delta_g.svg", bbox_inches="tight", format="svg")
    plt.close()

if __name__ == "__main__":
    main()
