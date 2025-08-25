import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 24


def main():
    os.makedirs("figures", exist_ok=True)

    # valid PDB pairs
    ok_pdb = pd.read_csv("../notebooks/ok_pdb.csv")
    ok_pdb["p_pdb_id"] = ok_pdb["p_pdb_id"].str.lower()
    ok_pdb["q_pdb_id"] = ok_pdb["q_pdb_id"].str.lower()
    ok_keys = set(zip(ok_pdb["p_pdb_id"], ok_pdb["q_pdb_id"]))

    k_values = [2, 3, 4, 5]

    # --------------------------------
    # SH+ILO (Simulation) dataset
    # --------------------------------
    exact_dfs = []
    exact_df = pd.read_csv("rmsdh_result/simulation_sh_combined.csv")
    shi_df = pd.read_csv("rmsdh_result/simulation_shibuya_combined.csv")
    for df in (exact_df, shi_df):
        df["p_pdb_id"] = df["p_pdb_id"].str.lower()

    for k in k_values:
        # load delta_g
        delta = pd.read_csv(f"delta_g_simulation_{k}_sigma0.5.csv")
        delta["p_pdb_id"] = delta["p_pdb_id"].str.lower()
        # theoretical error
        delta["theoretical_error"] = np.sqrt(
            2
            * delta["delta_g"]
            * np.log2(delta["Residue length"])
            * k
            / delta["Residue length"]
        )
        # subset RMSDh values
        sim_k = exact_df[exact_df["k"] == k][["p_pdb_id", "RMSDh"]].rename(
            columns={"RMSDh": "RMSDh_exact"}
        )
        shi_k = shi_df[shi_df["k"] == k][["p_pdb_id", "RMSDh"]].rename(
            columns={"RMSDh": "RMSDh_shi"}
        )
        # merge and compute actual error
        df_m = delta.merge(sim_k, on=["p_pdb_id"]).merge(
            shi_k, on=["p_pdb_id"]
        )
        df_m["actual_error"] = df_m["RMSDh_exact"] - df_m["RMSDh_shi"]
        df_m["k"] = k
        exact_dfs.append(df_m[["theoretical_error", "actual_error", "k"]])
    sim_error_df = pd.concat(exact_dfs, ignore_index=True)

    # --------------------------------
    # Shibuya dataset
    # --------------------------------
    shibuya_exact = {
        k: pd.read_csv(f"rmsdh_result/rmsdh_hingek_cnt_{k}.csv") for k in k_values
    }
    shibuya_fast = {
        k: pd.read_csv(f"rmsdh_result/fast_rmsdh_hingek_cnt_{k}.csv") for k in k_values
    }
    for d in (shibuya_exact, shibuya_fast):
        for k, df in d.items():
            df["p_pdb_id"] = df["p_pdb_id"].str.lower()
            df["q_pdb_id"] = df["q_pdb_id"].str.lower()

    delta_sh = pd.read_csv("delta_g_shibuya.csv")
    delta_sh["p_pdb_id"] = delta_sh["p_pdb_id"].str.lower()
    delta_sh["q_pdb_id"] = delta_sh["q_pdb_id"].str.lower()

    shi_dfs = []
    for k in k_values:
        d0 = delta_sh.copy()
        d0["theoretical_error"] = np.sqrt(
            2 * d0["delta_g"] * np.log2(d0["Residue length"]) * k / d0["Residue length"]
        )
        ex = shibuya_exact[k][["p_pdb_id", "q_pdb_id", "RMSDh"]].rename(
            columns={"RMSDh": "RMSDh_exact"}
        )
        fa = shibuya_fast[k][["p_pdb_id", "q_pdb_id", "RMSDh"]].rename(
            columns={"RMSDh": "RMSDh_fast"}
        )
        df_m = d0.merge(ex, on=["p_pdb_id", "q_pdb_id"]).merge(
            fa, on=["p_pdb_id", "q_pdb_id"]
        )
        df_m["actual_error"] = df_m["RMSDh_fast"] - df_m["RMSDh_exact"]
        df_m["k"] = k
        shi_dfs.append(df_m[["theoretical_error", "actual_error", "k"]])
    shi_error_df = pd.concat(shi_dfs, ignore_index=True)

    # --------------------------------
    # PAR dataset
    # --------------------------------
    par_exact = {
        k: filter_by_ok_keys(
            pd.read_csv(f"rmsdh_result/rmsdhk_more_data_{k}.csv"), ok_keys
        )
        for k in k_values
    }
    par_fast = {
        k: filter_by_ok_keys(
            pd.read_csv(f"rmsdh_result/fast_rmsdhk_more_data_{k}.csv"), ok_keys
        )
        for k in k_values
    }
    for d in (par_exact, par_fast):
        for k, df in d.items():
            df["p_pdb_id"] = df["p_pdb_id"].str.lower()
            df["q_pdb_id"] = df["q_pdb_id"].str.lower()

    delta_par = pd.read_csv("delta_g_par.csv")
    delta_par["p_pdb_id"] = delta_par["p_pdb_id"].str.lower()
    delta_par["q_pdb_id"] = delta_par["q_pdb_id"].str.lower()
    # filter by ok_keys
    delta_par = delta_par[
        delta_par.apply(lambda r: (r["p_pdb_id"], r["q_pdb_id"]) in ok_keys, axis=1)
    ]

    par_dfs = []
    for k in k_values:
        d0 = delta_par.copy()
        d0["theoretical_error"] = np.sqrt(
            2 * d0["delta_g"] * np.log2(d0["Residue length"]) * k / d0["Residue length"]
        )
        ex = par_exact[k][["p_pdb_id", "q_pdb_id", "RMSDh"]].rename(
            columns={"RMSDh": "RMSDh_exact"}
        )
        fa = par_fast[k][["p_pdb_id", "q_pdb_id", "RMSDh"]].rename(
            columns={"RMSDh": "RMSDh_fast"}
        )
        df_m = d0.merge(ex, on=["p_pdb_id", "q_pdb_id"]).merge(
            fa, on=["p_pdb_id", "q_pdb_id"]
        )
        df_m["actual_error"] = df_m["RMSDh_fast"] - df_m["RMSDh_exact"]
        df_m["k"] = k
        print(len(df_m))
        par_dfs.append(df_m[["theoretical_error", "actual_error", "k"]])
    par_error_df = pd.concat(par_dfs, ignore_index=True)

    # --------------------------------
    # DynDom dataset
    # --------------------------------
    dyn_exact = {
        k: pd.read_csv(f"rmsdh_result/rmsdhk_dyndom_data_{k}.csv") for k in k_values
    }
    dyn_fast = {
        k: pd.read_csv(f"rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_{k}.csv") for k in k_values
    }
    for d in (dyn_exact, dyn_fast):
        for k, df in d.items():
            df["p_pdb_id"] = df["p_pdb_id"].str.lower()
            df["q_pdb_id"] = df["q_pdb_id"].str.lower()

    delta_dy = pd.read_csv("delta_g_dyndom.csv")
    delta_dy["p_pdb_id"] = delta_dy["p_pdb_id"].str.lower()
    delta_dy["q_pdb_id"] = delta_dy["q_pdb_id"].str.lower()

    dy_dfs = []
    for k in k_values:
        d0 = delta_dy.copy()
        d0["theoretical_error"] = np.sqrt(
            2 * d0["delta_g"] * np.log2(d0["Residue length"]) * k / d0["Residue length"]
        )
        ex = dyn_exact[k][["p_pdb_id", "q_pdb_id", "RMSDh"]].rename(
            columns={"RMSDh": "RMSDh_exact"}
        )
        fa = dyn_fast[k][["p_pdb_id", "q_pdb_id", "RMSDh"]].rename(
            columns={"RMSDh": "RMSDh_fast"}
        )
        df_m = pd.concat([d0, ex, fa], axis=1, join="inner")
        df_m["actual_error"] = df_m["RMSDh_fast"] - df_m["RMSDh_exact"]
        df_m["k"] = k
        print(len(df_m))
        dy_dfs.append(df_m[["theoretical_error", "actual_error", "k"]])
    dyndom_error_df = pd.concat(dy_dfs, ignore_index=True)

    # ---------------------------
    # 4行4列の散布図プロット
    # ---------------------------
    datasets = [
        ("Simulation", sim_error_df),
        ("Shibuya 2008", shi_error_df),
        ("PAR 2020",    par_error_df),
        ("DynDom 2024", dyndom_error_df),
    ]

    # Figureとaxesを4x4で作成
    _, axes = plt.subplots(
        nrows=4, ncols=4,
        figsize=(20, 20)
    )


    # 各行・各列に散布図をプロット
    for i, (name, df_all) in enumerate(datasets):
        for j, k in enumerate(k_values):
            ax = axes[i, j]
            df = df_all[df_all["k"] == k]
            ax.scatter(df["theoretical_error"], df["actual_error"],
                       s=10, alpha=0.7)
            ax.set_xlim(0, 200)
            ax.set_xticks(np.arange(0, 200, 50))
            ax.set_ylim(0, df["actual_error"].max())
            ax.grid(True)
            # 左端の列には行ラベル（データセット名）を y 軸ラベルとして表示
            if j == 0:
                ax.set_ylabel(f"{name}\n"+r"Actual error ($\AA$)")
            # 最下段の行には x 軸ラベルを表示
            if j == 0:
                ax.set_xlabel(r"Theoretical error bound ($\AA$)"+"\n(a)")
            elif j == 1:
                ax.set_xlabel(r"Theoretical error bound ($\AA$)"+"\n(b)")
            elif j == 2:
                ax.set_xlabel(r"Theoretical error bound ($\AA$)"+"\n(c)")
            else:
                ax.set_xlabel(r"Theoretical error bound ($\AA$)"+"\n(d)")

    plt.tight_layout()
    plt.savefig("figures/all_sh_error_scatter_4x4.svg", format="svg")
    plt.savefig("figures/all_sh_error_scatter_4x4.png")
    plt.close()



def filter_by_ok_keys(df, ok_keys):
    return df[
        df.apply(
            lambda row: (row["p_pdb_id"].lower(), row["q_pdb_id"].lower()) in ok_keys,
            axis=1,
        )
    ]


if __name__ == "__main__":
    main()
