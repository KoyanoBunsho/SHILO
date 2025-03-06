import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18


def main():
    os.makedirs("figures", exist_ok=True)
    ok_pdb = pd.read_csv("../notebooks/ok_pdb.csv")
    small_monge = pd.read_csv("monge_rate_shibuya.csv")
    more_monge = pd.read_csv("monge_rate_par.csv")
    more_monge = more_monge[
        more_monge["p_pdb_id"].isin(ok_pdb["p_pdb_id"])
    ].reset_index()
    dyndom_monge = pd.read_csv("monge_rate_dyndom.csv")
    _, axes = plt.subplots(1, 4, figsize=(20, 5))
    simulation_df = pd.DataFrame()
    for idx, k in enumerate(range(2, 6)):
        df = pd.read_csv(f"monge_rate_simulation_{k}_sigma0.5.csv")
        simulation_df = pd.concat([simulation_df, df])
    print(len(simulation_df))
    print(simulation_df["Residue length"].mean())
    print(simulation_df["monge_rate_1"].mean())
    print(simulation_df["monge_rate_1"].std(ddof=0))
    print(simulation_df["monge_rate_1"].max())
    print(simulation_df["monge_rate_1"].min())
    _, axes = plt.subplots(nrows=1, ncols=4, figsize=(18, 3))
    simulation_df["monge_rate_1"].hist(
        ax=axes[0], bins=[0.4, 0.5, 0.6, 0.7, 0.8, 0.900, 1.00], color="Black"
    )
    axes[0].set_xlabel("(a)")
    axes[0].set_yticks([0, 5000, 10000])
    axes[0].set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    axes[0].set_xlim(0.4, 1.0)
    axes[0].set_ylim(0, 12000)
    small_monge["monge_rate_1"].hist(
        ax=axes[1], bins=[0.4, 0.5, 0.6, 0.7, 0.8, 0.900, 1.00], color="Black"
    )
    axes[1].set_xlabel("(b)")
    axes[1].set_yticks([0, 2, 4])
    axes[1].set_xlim(0.4, 1.0)
    axes[1].set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    axes[1].set_ylim(0, 5)
    more_monge["monge_rate_1"].hist(
        ax=axes[2], bins=[0.4, 0.5, 0.6, 0.7, 0.8, 0.900, 1.00], color="Black"
    )
    axes[2].set_xlabel("(c)")
    axes[2].set_xlim(0.4, 1.0)
    axes[2].set_yticks([0, 20, 40])
    axes[2].set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    axes[2].set_ylim(0, 50)
    dyndom_monge["monge_rate_1"].hist(
        ax=axes[3], bins=[0.4, 0.5, 0.6, 0.7, 0.8, 0.900, 1.00], color="Black"
    )
    axes[3].set_xlabel("(d)")
    axes[3].set_xlim(0.4, 1.0)
    axes[3].set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    axes[3].set_yticks([0, 500, 1000])
    axes[3].set_ylim(0, 1250)
    plt.savefig(
        "figures/distribution_of_monge_rate.svg", bbox_inches="tight", format="svg"
    )
    plt.tight_layout()
    plt.savefig("figures/simulation_monge_rate.png", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
