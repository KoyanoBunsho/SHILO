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

if __name__ == "__main__":
    main()
