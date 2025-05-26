import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['mathtext.fontset'] = 'stix'

def main():
    delta_g_simulation_2_df = pd.read_csv("delta_g_simulation_2_sigma0.5.csv")
    delta_g_simulation_3_df = pd.read_csv("delta_g_simulation_3_sigma0.5.csv")
    delta_g_simulation_4_df = pd.read_csv("delta_g_simulation_4_sigma0.5.csv")
    delta_g_simulation_5_df = pd.read_csv("delta_g_simulation_5_sigma0.5.csv")

    delta_g_shibuya_df = pd.read_csv("delta_g_shibuya.csv")
    delta_g_par_df     = pd.read_csv("delta_g_par.csv")
    delta_g_dyndom_df  = pd.read_csv("delta_g_dyndom.csv")

    delta_g_simulation_df = pd.concat([
        delta_g_simulation_2_df,
        delta_g_simulation_3_df,
        delta_g_simulation_4_df,
        delta_g_simulation_5_df
    ], ignore_index=True)

    dataframes = [
        (delta_g_simulation_df, "(a)"),
        (delta_g_shibuya_df,   "(b)"),
        (delta_g_par_df,       "(c)"),
        (delta_g_dyndom_df,    "(d)")
    ]

    fig, axes = plt.subplots(
        nrows=1, ncols=4,
        figsize=(16, 4),
        tight_layout=True
    )

    for i, (ax, (df, title)) in enumerate(zip(axes, dataframes)):
        ax.hist(df["delta_g"], bins=np.arange(0, 10000, 1000),
                color="black", edgecolor="white")
        ax.set_xlabel(title)
        if i == 0:
            ax.set_ylabel("#Data")

    plt.savefig("delta_g_histograms.svg", format="svg")
    plt.savefig("delta_g_histograms.png")
    plt.close(fig)

    stats = []
    for df, title in dataframes:
        mn = df["delta_g"].min()
        mean = df["delta_g"].mean()
        mx = df["delta_g"].max()
        stats.append((title, mn, mean, mx))

    tex_path = "delta_g_stats.tex"
    with open(tex_path, "w") as tex:
        tex.write(r"\begin{table}[ht]" + "\n")
        tex.write(r"  \centering" + "\n")
        tex.write(r"  \begin{tabular}{lrrr}" + "\n")
        tex.write(r"    \hline" + "\n")
        tex.write(r"    Data set & Min & Mean & Max \\" + "\n")
        tex.write(r"    \hline" + "\n")
        for title, mn, mean, mx in stats:
            tex.write(f"    {title} & {mn:.2f} & {mean:.2f} & {mx:.2f} \\\\\n")
        tex.write(r"    \hline" + "\n")
        tex.write(r"  \end{tabular}" + "\n")
        tex.write(r"  \caption{各データセットにおける $\Delta g$ の最小値、平均値、最大値}" + "\n")
        tex.write(r"  \label{tab:delta_g_stats}" + "\n")
        tex.write(r"\end{table}" + "\n")

    print(f"LaTeX テーブルを '{tex_path}' に書き出しました。")
    # ── ここまで既存の stats 出力 ──

    # 1000 刻みの度数分布を集計し，10000 より大きい値もカウントして .tex 形式で出力
    # -----------------------------------------
    # ビンのエッジを 0,1000,…,10000 とし，
    # np.histogram で度数を計算
    bin_edges = np.arange(0, 10001, 1000)  # [0,1000,2000,…,10000]
    counts_list = []
    for df, title in dataframes:
        # 各ビン内の度数
        counts, _ = np.histogram(df["delta_g"], bins=bin_edges)
        # 10000 より大きい値を個別にカウント
        over_count = int((df["delta_g"] > 10000).sum())
        counts_list.append((title, counts, over_count))

    counts_tex_path = "delta_g_counts.tex"
    with open(counts_tex_path, "w") as tex:
        tex.write(r"\begin{table}[ht]" + "\n")
        tex.write(r"  \centering" + "\n")
        # カラム数：Data set + 10 ビン + >10000
        tex.write(r"  \begin{tabular}{l" + "r" * (len(bin_edges)-1 + 1) + "}" + "\n")
        tex.write(r"    \hline" + "\n")
        # ヘッダー行：ビン区間のラベル
        bin_labels = [f"{int(bin_edges[i])}\\text{'–'}{int(bin_edges[i+1])}" for i in range(len(bin_edges)-1)]
        header = "Data set & " + " & ".join(bin_labels) + " & >10000 \\\\"
        tex.write(f"    {header}\n")
        tex.write(r"    \hline" + "\n")
        # 各データセットごとの度数を出力
        for title, counts, over in counts_list:
            counts_str = " & ".join(str(c) for c in counts)
            tex.write(f"    {title} & {counts_str} & {over} \\\\\n")
        tex.write(r"    \hline" + "\n")
        tex.write(r"  \end{tabular}" + "\n")
        tex.write(r"  \caption{The distribution of $\Delta_G$ values for all graphs $G$ considered in the (a) Simulation, (b) Shibuya 2008, (c) PAR 2020, and (d) DynDom2024 datasets.}" + "\n")
        tex.write(r"  \label{tab:delta_g_hist_kn}" + "\n")
        tex.write(r"\end{table}" + "\n")

    print(f"LaTeX カウントテーブルを '{counts_tex_path}' に書き出しました。")
    # -----------------------------------------

if __name__ == "__main__":
    main()
