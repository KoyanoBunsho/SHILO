import pandas as pd

def main():
    csv_files = [
        "all_pdb/dyndom_execution_time_dyn_improved.csv",
        "all_pdb/fatcat_execution_time_dyn_improved.csv",
        "all_pdb/dyndom_execution_time_par_improved.csv",
        "all_pdb/fatcat_execution_time_par_improved.csv"
    ]
    # なぜかFATCATで計算できないケースがある
    all_execution_times = []
    print("【各CSVファイルごとの execution_time の平均値と合計値】")
    for file in csv_files:
        df = pd.read_csv(file)
        avg = df["execution_time"].mean()
        total = df["execution_time"].sum()
        print(f"{file} -> 平均値: {avg}, 合計値: {total}")
        all_execution_times.append(df["execution_time"])
    combined_execution_times = pd.concat(all_execution_times, ignore_index=True)
    overall_avg = combined_execution_times.mean()
    overall_total = combined_execution_times.sum()    
    print("\n【全CSVファイルの execution_time の総合計結果】")
    print(f"全体の平均値: {overall_avg}")
    print(f"全体の合計値: {overall_total}")

if __name__ == "__main__":
    main()
