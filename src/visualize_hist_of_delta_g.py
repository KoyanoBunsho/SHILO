import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# フォント設定
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['font.size'] = 24

def main():
    # CSV 読み込み
    delta_g_simulation_df = pd.read_csv("rmsdh_result/simulation_delta_g_output_simulation_file.csv")

    # 描画対象が1つだけなので、直接取り出し
    df = delta_g_simulation_df
    title = "Simulation"

    # サブプロット1つ
    fig, ax = plt.subplots(
        nrows=1, ncols=1,
        figsize=(4, 4),
        tight_layout=True
    )

    # ヒストグラム描画
    ax.hist(
        df["delta_g"],
        bins=np.arange(0, 10000, 1000),
        color="black",
        edgecolor="white"
    )
    ax.set_title(title, fontsize=12)
    ax.set_ylabel("#Data", fontsize=10)

    # SVG & PNGで保存
    plt.savefig("delta_g_histograms_simulation.svg", format="svg")
    plt.savefig("delta_g_histograms_simulation.png", dpi=300)
    plt.close(fig)

    # ――― ここから TeX 出力の追記 ―――
    # 1000 刻みのビンエッジを 0,1000,…,10000 とし
    # 各ビンの度数と 10000 を超える値をカウント
    bin_edges = np.arange(0, 10001, 2000)  # [0,1000,2000,…,10000]
    counts, _ = np.histogram(df["delta_g"], bins=bin_edges)
    over_count = int((df["delta_g"] > 10000).sum())

    tex_path = "delta_g_counts_simulation.tex"
    with open(tex_path, "w") as tex:
        tex.write(r"\begin{table}[ht]" + "\n")
        tex.write(r"  \centering" + "\n")
        # Data set + 10ビン + >10000 列
        tex.write(r"  \begin{tabular}{l" + "r" * (len(bin_edges)-1 + 1) + "}" + "\n")
        tex.write(r"    \hline" + "\n")
        # ヘッダー：区間ラベル
        bin_labels = [f"[{int(bin_edges[i])}, {int(bin_edges[i+1])})" for i in range(len(bin_edges)-1)]
        header = "Dataset & " + " & ".join(bin_labels) + " & [10000, $\infty$)] \\\\"
        tex.write(f"    {header}\n")
        tex.write(r"    \hline" + "\n")
        # シミュレーションデータの度数を出力
        counts_str = " & ".join(str(c) for c in counts)
        tex.write(f"    {title} & {counts_str} & {over_count} \\\\\n")
        tex.write(r"    \hline" + "\n")
        tex.write(r"  \end{tabular}" + "\n")
        tex.write(r"  \caption{The distribution of $\Delta_G$ values for all graphs $G$ considered in the simulation dataset.}" + "\n")
        tex.write(r"  \label{tab:delta_g_simulation}" + "\n")
        tex.write(r"\end{table}" + "\n")

    print(f"LaTeX カウントテーブルを '{tex_path}' に書き出しました。")
    # ――― ここまで ―――

if __name__ == "__main__":
    main()
