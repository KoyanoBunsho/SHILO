import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18

def main():
    os.makedirs("figures", exist_ok=True)
    
    # ok_pdb の読み込みと前処理（pdb ID を小文字に統一）
    ok_pdb = pd.read_csv("../notebooks/ok_pdb.csv")
    ok_pdb["p_pdb_id"] = ok_pdb["p_pdb_id"].str.lower()
    ok_pdb["q_pdb_id"] = ok_pdb["q_pdb_id"].str.lower()
    ok_keys = set(zip(ok_pdb["p_pdb_id"], ok_pdb["q_pdb_id"]))
    
    # k の範囲（例として 2～10）
    k_values = range(2, 11)
    
    shilo_times = []   # SH+ILO の計算時間（合計値）
    dyndom_times = []  # DynDom の計算時間（合計値）
    
    for k in k_values:
        # SH+ILO のデータ読み込み
        shilo_file = f"rmsdh_result/fast_rmsdhk_more_data_{k}_pospro_loop.csv"
        try:
            df_shilo = pd.read_csv(shilo_file)
        except FileNotFoundError:
            print(f"File not found: {shilo_file}")
            shilo_times.append(None)
        else:
            # ok_keys に登録されたペアのみ抽出
            df_shilo = df_shilo[df_shilo.apply(lambda row: (row["p_pdb_id"], row["q_pdb_id"]) in ok_keys, axis=1)]
            if "exec_time (s)" in df_shilo.columns:
                total_time_shilo = df_shilo["exec_time (s)"].sum()
            else:
                total_time_shilo = None
            shilo_times.append(total_time_shilo)
            if total_time_shilo is not None:
                print(f"k = {k}: SH+ILO Computation Time = {total_time_shilo:.3f} s")
            else:
                print(f"k = {k}: SH+ILO Computation Time = No data")
        
        # DynDom のデータ読み込み
        dyndom_file = f"rmsdh_result/fast_rmsdhk_dyndom_{k}_postpro_loop.csv"
        try:
            df_dyndom = pd.read_csv(dyndom_file)
        except FileNotFoundError:
            print(f"File not found: {dyndom_file}")
            dyndom_times.append(None)
        else:
            if "exec_time (s)" in df_dyndom.columns:
                total_time_dyndom = df_dyndom["exec_time (s)"].sum()
            else:
                total_time_dyndom = None
            dyndom_times.append(total_time_dyndom)
            if total_time_dyndom is not None:
                print(f"k = {k}: DynDom Computation Time = {total_time_dyndom:.3f} s")
            else:
                print(f"k = {k}: DynDom Computation Time = No data")
    
    # 1行2列のサブプロット作成
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 左側：SH+ILO のプロット
    axes[0].plot(list(k_values), shilo_times, marker='o', linestyle='-')
    axes[0].set_ylabel("Computation Time (s)")
    axes[0].set_xticks(list(k_values))
    axes[0].grid(True)
    
    # 右側：DynDom のプロット
    axes[1].plot(list(k_values), dyndom_times, marker='s', linestyle='-')
    axes[1].set_ylabel("Computation Time (s)")
    axes[1].set_xticks(list(k_values))
    axes[1].grid(True)
    
    plt.tight_layout()
    plt.savefig("figures/computation_time_shilo_dyndom.svg", format="svg")
    plt.savefig("figures/computation_time_shilo_dyndom.png")
    plt.close()

if __name__ == "__main__":
    main()
