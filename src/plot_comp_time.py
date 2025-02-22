import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18


def main():
    os.makedirs("figures", exist_ok=True)
    comp_time_cols_list = [f"{i}_computation_time" for i in range(100)]
    rlo_par_dict, rlo_dyndom_dict = {}, {}
    rilo_par_dict, rilo_dyndom_dict = {}, {}
    sh_par_dict, sh_dyndom_dict = {}, {}
    sh_lo_par_dict, sh_lo_dyndom_dict = {}, {}
    shilo_par_dict, shilo_dyndom_dict = {}, {}
    k_max = 10
    for k in range(2, k_max + 1):
        rlo_par_dict[k] = pd.read_csv(f"rmsdh_result/ablation_study_{k}_paper.csv")
        rlo_dyndom_dict[k] = pd.read_csv(f"rmsdh_result/ablation_study_dyndom_{k}.csv")

        rilo_par_dict[k] = pd.read_csv(f"rmsdh_result/ablation_study_loop_par{k}.csv")
        rilo_dyndom_dict[k] = pd.read_csv(
            f"rmsdh_result/ablation_study_loop_dyndom{k}.csv"
        )

        sh_par_dict[k] = pd.read_csv(f"rmsdh_result/fast_rmsdhk_more_data_{k}.csv")
        sh_dyndom_dict[k] = pd.read_csv(
            f"rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_{k}.csv"
        )

        sh_lo_par_dict[k] = pd.read_csv(
            f"rmsdh_result/fast_rmsdhk_more_data_{k}_pospro.csv"
        )
        sh_lo_dyndom_dict[k] = pd.read_csv(
            f"rmsdh_result/fast_rmsdh_hingek_cnt_dyndom_{k}_postpro.csv"
        )

        shilo_par_dict[k] = pd.read_csv(
            f"rmsdh_result/fast_rmsdhk_more_data_{k}_pospro_loop.csv"
        )
        shilo_dyndom_dict[k] = pd.read_csv(
            f"rmsdh_result/fast_rmsdhk_dyndom_{k}_postpro_loop.csv"
        )

    comp_time_par = {"R+LO": [], "R+ILO": [], "SH": [], "SH+LO": [], "SH+ILO": []}
    comp_time_dyndom = {"R+LO": [], "R+ILO": [], "SH": [], "SH+LO": [], "SH+ILO": []}
    for k in range(2, k_max + 1):
        comp_time_par["R+LO"].append(
            rlo_par_dict[k][comp_time_cols_list].sum().mean() / 1000
        )
        comp_time_par["R+ILO"].append(
            rilo_par_dict[k][comp_time_cols_list].sum().mean() / 1000
        )
        comp_time_par["SH"].append(sh_par_dict[k]["exec_time (s)"].sum())
        comp_time_par["SH+LO"].append(sh_lo_par_dict[k]["exec_time (s)"].sum())
        comp_time_par["SH+ILO"].append(shilo_par_dict[k]["exec_time (s)"].sum())

        comp_time_dyndom["R+LO"].append(
            rlo_dyndom_dict[k][comp_time_cols_list].sum().mean() / 1000
        )
        comp_time_dyndom["R+ILO"].append(
            rilo_dyndom_dict[k][comp_time_cols_list].sum().mean() / 1000
        )
        comp_time_dyndom["SH"].append(sh_dyndom_dict[k]["exec_time (s)"].sum())
        comp_time_dyndom["SH+LO"].append(sh_lo_dyndom_dict[k]["exec_time (s)"].sum())
        comp_time_dyndom["SH+ILO"].append(shilo_dyndom_dict[k]["exec_time (s)"].sum())
    k_values = range(2, k_max + 1)
    _, axes = plt.subplots(1, 2, figsize=(14, 6))
    for key, values in comp_time_par.items():
        axes[0].plot(k_values, values, marker="o", linestyle="-", label=key)
    axes[0].set_xlabel("k")
    axes[0].set_xticks(range(2, k_max + 1))
    axes[0].set_ylim(0)
    axes[0].set_ylabel("Computation Time (s)")
    axes[0].set_title("PAR 2020")
    axes[0].legend()
    axes[0].grid(True)
    for key, values in comp_time_dyndom.items():
        axes[1].plot(k_values, values, marker="o", linestyle="-", label=key)
    axes[1].set_xlabel("k")
    axes[1].set_xticks(range(2, k_max + 1))
    axes[1].set_ylabel("Computation Time (s)")
    axes[1].set_ylim(0)
    axes[1].set_title("DynDom 2024")
    axes[1].legend()
    axes[1].grid(True)
    plt.tight_layout()
    plt.savefig("figures/computation_time_comparison.svg", format="svg")
    plt.savefig("figures/computation_time_comparison.png")
    plt.close()


if __name__ == "__main__":
    main()
