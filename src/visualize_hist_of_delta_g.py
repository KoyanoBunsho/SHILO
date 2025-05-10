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
        (delta_g_simulation_df, "Simulation"),
        (delta_g_shibuya_df,   "Shibuya 2008"),
        (delta_g_par_df,       "PAR 2020"),
        (delta_g_dyndom_df,    "DynDom 2024")
    ]

    fig, axes = plt.subplots(
        nrows=1, ncols=4,
        figsize=(16, 4),  # お好みで調整
        tight_layout=True
    )

    for i, (ax, (df, title)) in enumerate(zip(axes, dataframes)):
        ax.hist(df["delta_g"], bins=np.arange(0, 10000, 1000), color="black", edgecolor="white")
        ax.set_title(title)
        if i == 0:
            ax.set_ylabel("#Data")

    plt.savefig("delta_g_histograms.svg", format="svg")
    plt.savefig("delta_g_histograms.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
