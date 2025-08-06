import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# フォント設定
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 28

def main():
    methods = ["shibuya", "par", "dyndom"]
    _, axes = plt.subplots(2, 2, figsize=(10, 10)) 
    axes_flat = axes.ravel()

    for i, (ax, method) in enumerate(zip(axes_flat[:3], methods)):
        delta_g_df = pd.read_csv(f"delta_g_{method}.csv")
        avg_delta_df = pd.read_csv(f"average_delta_g_{method}.csv")
        df = pd.concat([delta_g_df, avg_delta_df], axis=1)

        # 散布図
        ax.scatter(df["delta_g"], df["average_delta_g"], s=1)
        ax.set_xlim(0, df["delta_g"].max())
        ax.set_ylim(0, df["average_delta_g"].max())
        ax.set_xlabel(r"$\Delta_G$" + f"\n({chr(ord('a') + i)})")
        ax.set_ylabel(r"$\overline{\Delta_G}$")
        ax.grid(True)

        # 追加する直線
        x_vals = np.linspace(0, df["delta_g"].max(), 100)
        for factor in [1, 0.1, 0.01, 0.001]:
            ax.plot(x_vals, x_vals * factor, linestyle="--")

        # 相関係数の計算と出力
        corr = df["delta_g"].corr(df["average_delta_g"])
        print(f"Correlation for {method}: {corr:.4f}")
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
    ax = axes_flat[3]
    ax.scatter(df_sim_all["delta_g"], df_sim_all["average_delta_g"], s=1)
    ax.set_xlim(0, df_sim_all["delta_g"].max())
    ax.set_ylim(0, df_sim_all["average_delta_g"].max())
    ax.set_xlabel(r"$\Delta_G$" + "\n(d)")
    ax.set_ylabel(r"$\overline{\Delta_G}$")
    ax.grid(True)

    # 追加する直線（シミュレーション）
    x_vals_sim = np.linspace(0, df_sim_all["delta_g"].max(), 100)
    for factor in [1, 0.1, 0.01, 0.001]:
        ax.plot(x_vals_sim, x_vals_sim * factor, linestyle="--")

    # 相関係数の計算と出力（シミュレーション）
    corr_sim = df_sim_all["delta_g"].corr(df_sim_all["average_delta_g"])
    print(f"Correlation for simulation: {corr_sim:.4f}")

    # レイアウト調整・保存
    plt.tight_layout()
    plt.savefig("comparison_average_delta_g.png", bbox_inches="tight")
    plt.savefig("comparison_average_delta_g.svg", bbox_inches="tight", format="svg")
    plt.close()

if __name__ == "__main__":
    main()
