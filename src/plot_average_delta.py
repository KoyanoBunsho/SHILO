import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# フォント設定
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18

def main():
    methods = ["shibuya", "par", "dyndom"]
    _, axes = plt.subplots(2, 2, figsize=(10, 10)) 
    axes_flat = axes.ravel()
    linestyles = ["-", "--", ":", "-."]  # 実線、破線、点線、‐．点線
    colors     = ["red", 'green',  'magenta',    'black']
    labels   = [r"$y=x$", r'$y=\frac{x}{10}$', r'$y=\frac{x}{100}$', r'$y=\frac{x}{1000}$']
    for i, (ax, method) in enumerate(zip(axes_flat[:3], methods)):
        delta_g_df = pd.read_csv(f"delta_g_{method}.csv")
        avg_delta_df = pd.read_csv(f"average_delta_g_{method}.csv")
        df = pd.concat([delta_g_df, avg_delta_df], axis=1)
        print(method)
        print((df["average_delta_g"] / df["delta_g"]).mean())

        # 散布図
        xticks = np.linspace(0, df["delta_g"].max(), 4)        # 0 〜 max を 11 点に分割
        yticks = np.linspace(0, df["average_delta_g"].max(), 4)
        round_base = 100
        xticks = np.round(xticks / round_base) * round_base
        if i == 0:
            yticks = np.round(yticks / 10) * 10
        else:
            yticks = np.round(yticks / round_base) * round_base
        ax.scatter(df["delta_g"], df["average_delta_g"], s=1)
        ax.set_xlim(0, df["delta_g"].max())
        ax.set_ylim(0, df["average_delta_g"].max())
        ax.set_xlabel(r"$\Delta_G$" + f"\n({chr(ord('a') + i)})")
        ax.set_ylabel(r"$\tilde{\Delta_G}$")
        ax.set_xticks(xticks)
        if i == 2:
            ax.set_yticks([0, 1000, 2000, 3000, 4000])
        else:
            ax.set_yticks(yticks)
        ax.grid(True)

        # 追加する直線
        x_vals = np.linspace(0, df["delta_g"].max(), 100)
        #for factor, ls, color, label in zip([1, 0.1, 0.01, 0.001], linestyles, colors, labels):
            #ax.plot(x_vals, x_vals * factor, linestyle=ls, label=label, color=color)

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
    print("Simulation")
    print((df_sim_all["average_delta_g"] / df_sim_all["delta_g"]).mean())
    ax.scatter(df_sim_all["delta_g"], df_sim_all["average_delta_g"], s=1)
    ax.set_xlim(0, df_sim_all["delta_g"].max())
    ax.set_ylim(0, df_sim_all["average_delta_g"].max())
    ax.set_xlabel(r"$\Delta_G$" + "\n(d)")
    ax.set_ylabel(r"$\tilde{\Delta_G}$")
    xticks = np.linspace(0, df_sim_all["delta_g"].max(), 4)
    yticks = np.linspace(0, df_sim_all["average_delta_g"].max(), 4)
    xticks = np.linspace(0, df["delta_g"].max(), 4)        # 0 〜 max を 11 点に分割
    yticks = np.linspace(0, df["average_delta_g"].max(), 4)
    round_base = 132000
    xticks = np.round(xticks / round_base) * round_base
    yticks = np.round(yticks / 1000) * 1000
    ax.set_xticks(xticks)
    ax.set_yticks([0, 1000, 2000, 3000, 4000])
    ax.grid(True)

    # 追加する直線（シミュレーション）
    x_vals_sim = np.linspace(0, df_sim_all["delta_g"].max(), 100)
    #for factor, ls, color, label in zip([1, 0.1, 0.01, 0.001], linestyles, colors, labels):
        #ax.plot(x_vals_sim, x_vals_sim * factor, linestyle=ls, color=color, label=label)
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
