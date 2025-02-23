import pandas as pd
import os


def process_files(prefix):
    dfs = []
    for hinge_num in range(2, 11):
        for sigma in [0.5, 1.5]:
            filename = f"rmsdh_result/{prefix}_{hinge_num}_{sigma}.csv"
            if os.path.exists(filename):
                df = pd.read_csv(filename)
                dfs.append(df)
            else:
                print(f"Warning: {filename} が見つかりませんでした。")
    if dfs:
        return pd.concat(dfs, axis=0, ignore_index=True)
    else:
        return None


def main():
    prefixes = [
        "simulation_sh_lo",
        "simulation_sh",
        "simulation_shibuya",
        "simulation_r_lo",
        "simulation_r_ilo",
    ]
    combined_results = {}

    for prefix in prefixes:
        print(f"Processing files with prefix: {prefix}")
        df_combined = process_files(prefix)
        if df_combined is not None:
            combined_results[prefix] = df_combined
            print(f"{prefix} の結合結果：{len(df_combined)} 行")
            output_filename = f"rmsdh_result/{prefix}_combined.csv"
            df_combined.to_csv(output_filename, index=False)
        else:
            print(f"{prefix} のファイルは存在しませんでした。")
    for prefix, df in combined_results.items():
        print(f"\n{prefix} の先頭5行:")
        print(df.head())


if __name__ == "__main__":
    main()
