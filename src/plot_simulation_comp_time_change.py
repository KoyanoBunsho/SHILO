import pandas as pd
import matplotlib.pyplot as plt
import os
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18
def read_simulation_data(file_path):
    df = pd.read_csv(file_path).fillna("")
    df["actual_hinge_cnt"] = df["k"]
    df["primary_key"] = (
        df["p_pdb_id"].apply(lambda x: str(x)[:9])
        + "_hinge_"
        + df["actual_hinge_cnt"].astype(str)
        + "_sigma0.5"
    )
    return df

def main():
    os.makedirs("figures", exist_ok=True)
    
    # simulation_sh_ilo のデータ読み込み
    simulation_sh_ilo_df = read_simulation_data("rmsdh_result/simulation_sh_ilo_combined.csv")
    
    # k の値ごとに計算時間（秒）の平均値を算出
    ks = [2, 3, 4, 5]
    avg_times = []
    for k in ks:
        sub_df = simulation_sh_ilo_df[simulation_sh_ilo_df["actual_hinge_cnt"] == k]
        if "exec_time (s)" in sub_df.columns:
            avg_time = sub_df["exec_time (s)"].sum()
        else:
            avg_time = None
        avg_times.append(avg_time)
        print(f"k = {k}: Average Computation Time = {avg_time:.3f} s" if avg_time is not None else f"k = {k}: No data")
    
    # プロットの作成
    plt.figure(figsize=(8, 6))
    plt.plot(ks, avg_times, marker='o', linestyle='-')
    plt.ylabel("Computation Time (s)")
    plt.xticks(ks)
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig("figures/simulation_sh_ilo_computation_time.svg", format="svg")
    plt.savefig("figures/simulation_sh_ilo_computation_time.png")
    plt.close()

if __name__ == "__main__":
    main()
