import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18


def preprocess_df(df):
    df["p_pdb_id"] = df["p_pdb"].apply(lambda x: str(x).lower().split("_")[0])
    df["p_chain_id"] = df["p_pdb"].apply(lambda x: str(x).split("_")[1])
    df["q_pdb_id"] = df["q_pdb"].apply(lambda x: str(x).lower().split("_")[0])
    df["q_chain_id"] = df["q_pdb"].apply(lambda x: str(x).split("_")[1])
    return df


def format_stats(stats):
    return f"{stats:.3f}"


def main():
    os.makedirs("figures", exist_ok=True)
    ok_pdb = pd.read_csv("../notebooks/ok_pdb.csv")
    ok_pdb["p_pdb_id"] = ok_pdb["p_pdb_id"].str.lower()
    ok_pdb["q_pdb_id"] = ok_pdb["q_pdb_id"].str.lower()
    ok_keys = set(zip(ok_pdb["p_pdb_id"], ok_pdb["q_pdb_id"]))
    dyndom_comp_time_for_par_df = preprocess_df(
        pd.read_csv("all_pdb/dyndom_execution_time_par_improved.csv")
    )
    dyndom_comp_time_for_dyn_df = pd.read_csv("all_pdb/dyndom_execution_time_dyn_improved.csv")
    fatcat_comp_time_for_par_df = preprocess_df(
        pd.read_csv("all_pdb/fatcat_execution_time_par_improved.csv")
    )
    fatcat_comp_time_for_dyn_df = pd.read_csv("all_pdb/fatcat_execution_time_dyn_improved.csv")

    dyn_data_for_experiment_df = pd.read_csv("rmsdh_result/dyn_data_for_experiment.csv")
    dyn_data_for_experiment_df["p_q_pdb"] = (
        dyn_data_for_experiment_df["p_pdb_id"] + "_" + dyn_data_for_experiment_df["q_pdb_id"]
    )
    dyndom_comp_time_for_dyn_df["p_q_pdb"] = (
        dyndom_comp_time_for_dyn_df["p_pdb"] + "_" + dyndom_comp_time_for_dyn_df["q_pdb"]
    )
    fatcat_comp_time_for_dyn_df["p_q_pdb"] = (
        fatcat_comp_time_for_dyn_df["p_pdb"] + "_" + fatcat_comp_time_for_dyn_df["q_pdb"]
    )
    dyndom_comp_time_for_dyn_df = dyndom_comp_time_for_dyn_df[
        dyndom_comp_time_for_dyn_df["p_q_pdb"].isin(dyn_data_for_experiment_df["p_q_pdb"].to_list())
    ]
    fatcat_comp_time_for_dyn_df = fatcat_comp_time_for_dyn_df[
        fatcat_comp_time_for_dyn_df["p_q_pdb"].isin(dyn_data_for_experiment_df["p_q_pdb"].to_list())
    ]
    dyndom_comp_time_for_par_df = dyndom_comp_time_for_par_df[
        dyndom_comp_time_for_par_df.apply(
            lambda row: (row["p_pdb_id"], row["q_pdb_id"]) in ok_keys, axis=1
        )
    ]
    fatcat_comp_time_for_par_df = fatcat_comp_time_for_par_df[
        fatcat_comp_time_for_par_df.apply(
            lambda row: (row["p_pdb_id"], row["q_pdb_id"]) in ok_keys, axis=1
        )
    ]
    assert len(fatcat_comp_time_for_dyn_df) == len(dyndom_comp_time_for_dyn_df) == len(dyn_data_for_experiment_df)
    assert len(fatcat_comp_time_for_par_df) == len(dyndom_comp_time_for_par_df) == len(ok_pdb)
    comp_time_cols_list = [f"{i}_computation_time" for i in range(100)]
    rlo_par_dict, rlo_dyndom_dict = {}, {}
    rilo_par_dict, rilo_dyndom_dict = {}, {}
    sh_par_dict, sh_dyndom_dict = {}, {}
    sh_lo_par_dict, sh_lo_dyndom_dict = {}, {}
    shilo_par_dict, shilo_dyndom_dict = {}, {}

    k_max = 5
    for k in range(2, k_max + 1):
        df = pd.read_csv(f"rmsdh_result/ablation_study_par_{k}.csv")
        rlo_par_dict[k] = df[df.apply(lambda row: (row["p_pdb_id"], row["q_pdb_id"]) in ok_keys, axis=1)]
        df = pd.read_csv(f"rmsdh_result/ablation_study_loop_par{k}.csv")
        rilo_par_dict[k] = df[df.apply(lambda row: (row["p_pdb_id"], row["q_pdb_id"]) in ok_keys, axis=1)]
        df = pd.read_csv(f"rmsdh_result/fast_rmsdhk_more_data_{k}.csv")
        sh_par_dict[k] = df[df.apply(lambda row: (row["p_pdb_id"], row["q_pdb_id"]) in ok_keys, axis=1)]
        df = pd.read_csv(f"rmsdh_result/fast_rmsdhk_more_data_{k}_pospro.csv")
        sh_lo_par_dict[k] = df[df.apply(lambda row: (row["p_pdb_id"], row["q_pdb_id"]) in ok_keys, axis=1)]
        df = pd.read_csv(f"rmsdh_result/fast_rmsdhk_more_data_{k}_pospro_loop.csv")
        shilo_par_dict[k] = df[df.apply(lambda row: (row["p_pdb_id"], row["q_pdb_id"]) in ok_keys, axis=1)]
        rlo_dyndom_dict[k] = pd.read_csv(f"rmsdh_result/ablation_study_dyndom_{k}.csv")
        rilo_dyndom_dict[k] = pd.read_csv(f"rmsdh_result/ablation_study_loop_dyndom{k}.csv")
        sh_dyndom_dict[k] = pd.read_csv(f"rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_{k}.csv")
        sh_lo_dyndom_dict[k] = pd.read_csv(f"rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_{k}_postpro.csv")
        shilo_dyndom_dict[k] = pd.read_csv(f"rmsdh_result/fast_rmsdhk_dyndom_{k}_postpro_loop.csv")
    comp_time_par = {"R+LO": [], "R+ILO": [], "SH": [], "SH+LO": [], "SH+ILO": []}
    comp_time_dyndom = {"R+LO": [], "R+ILO": [], "SH": [], "SH+LO": [], "SH+ILO": []}
    for k in range(2, k_max + 1):
        comp_time_par["R+LO"].append(rlo_par_dict[k][comp_time_cols_list].mean().mean() / 1000)
        comp_time_par["R+ILO"].append(rilo_par_dict[k][comp_time_cols_list].mean().mean() / 1000)
        comp_time_par["SH"].append(sh_par_dict[k]["exec_time (s)"].mean())
        comp_time_par["SH+LO"].append(sh_lo_par_dict[k]["exec_time (s)"].mean())
        comp_time_par["SH+ILO"].append(shilo_par_dict[k]["exec_time (s)"].mean())
        comp_time_dyndom["R+LO"].append(rlo_dyndom_dict[k][comp_time_cols_list].mean().mean() / 1000)
        comp_time_dyndom["R+ILO"].append(rilo_dyndom_dict[k][comp_time_cols_list].mean().mean() / 1000)
        comp_time_dyndom["SH"].append(sh_dyndom_dict[k]["exec_time (s)"].mean())
        comp_time_dyndom["SH+LO"].append(sh_lo_dyndom_dict[k]["exec_time (s)"].mean())
        comp_time_dyndom["SH+ILO"].append(shilo_dyndom_dict[k]["exec_time (s)"].mean())

    marker_styles = {
        "R+LO": "o",
        "R+ILO": "s",
        "SH": "D",
        "SH+LO": "^",
        "SH+ILO": "v",
    }
    _, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6), sharey=True)
    k_values = list(range(2, k_max + 1))
    for key, values in comp_time_par.items():
        axes[0].plot(
            k_values, values,
            marker=marker_styles[key], linestyle="-", label=key
        )
    axes[0].set_xlabel("(a)", fontsize=28)
    axes[0].set_xticks(k_values)
    axes[0].set_ylabel("Computation Time (s)", fontsize=28)
    axes[0].legend()
    axes[0].grid(True)
    for key, values in comp_time_dyndom.items():
        axes[1].plot(
            k_values, values,
            marker=marker_styles[key], linestyle="-", label=key
        )
    axes[1].set_xlabel("(b)", fontsize=28)
    axes[1].set_xticks(k_values)
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()
    plt.savefig("figures/computation_time_comparison.svg", format="svg")
    plt.savefig("figures/computation_time_comparison.png")
    plt.close()


if __name__ == "__main__":
    main()
