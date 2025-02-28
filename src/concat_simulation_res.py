import csv
import os


def process_files(prefix):
    combined_rows = []
    header = None
    for hinge_num in range(2, 11):
        for sigma in [0.5]:
            filename = f"rmsdh_result/{prefix}_{hinge_num}_{sigma}.csv"
            if os.path.exists(filename):
                with open(filename, newline="", encoding="utf-8") as csvfile:
                    reader = csv.reader(csvfile)
                    try:
                        file_header = next(reader)
                    except StopIteration:
                        continue
                    if header is None:
                        header = file_header
                        combined_rows.append(header)
                    for row in reader:
                        combined_rows.append(row)
            else:
                print(f"Warning: {filename} が見つかりませんでした。")
    if header is not None:
        return combined_rows
    else:
        return None


def main():
    prefixes = [
        "simulation_sh_ilo",
        "simulation_sh_lo",
        "simulation_sh",
        "simulation_shibuya",
        "simulation_r_lo",
        "simulation_r_ilo",
    ]
    combined_results = {}
    for prefix in prefixes:
        print(f"Processing files with prefix: {prefix}")
        combined_data = process_files(prefix)
        if combined_data is not None:
            combined_results[prefix] = combined_data
            row_count = len(combined_data) - 1
            print(f"{prefix} の結合結果：{row_count} 行")
            output_filename = f"rmsdh_result/{prefix}_combined.csv"
            with open(
                output_filename, mode="w", newline="", encoding="utf-8"
            ) as csvfile:
                writer = csv.writer(csvfile)
                writer.writerows(combined_data)
        else:
            print(f"{prefix} のファイルは存在しませんでした。")
    for prefix, data in combined_results.items():
        print(f"\n{prefix} の先頭5行:")
        for row in data[:6]:
            print(row)


if __name__ == "__main__":
    main()
