import argparse
import pandas as pd
import subprocess
import os
import time
from joblib import Parallel, delayed


def main(csv_file, output_csv):
    if not os.path.exists(csv_file):
        print("Error: CSV file does not exist.")
        return
    pdb_pairs = load_pdb_pairs(csv_file)
    os.chdir("all_pdb")
    num_cores = os.cpu_count()
    results = Parallel(n_jobs=num_cores)(
        delayed(process_pair)(pair) for _, pair in pdb_pairs.iterrows()
    )
    execution_time_list = [
        {"p_pdb": pair["p_pdb"], "q_pdb": pair["q_pdb"], "execution_time": duration}
        for pair, duration in zip(pdb_pairs.to_dict("records"), results)
    ]
    pd.DataFrame(execution_time_list).to_csv(output_csv, index=False)
    print(f"Output file created successfully: {output_csv}")


def convert_pdb_id(pdb_str):
    pdb, chain = pdb_str.split("_")
    return f"{pdb.lower()}_{chain}"


def process_pair(pair):
    make_dyndom_command_file(pair)
    return run_dyndom(pair)


def load_pdb_pairs(csv_file):
    return pd.read_csv(csv_file)


def run_dyndom(pair):
    pdb1 = convert_pdb_id(pair["p_pdb"])
    pdb2 = convert_pdb_id(pair["q_pdb"])
    output_filename = f"{pdb1}_{pdb2}"
    command = ["./DynDom_speed", f"{output_filename}.w5.command"]
    start_time = time.time()
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        output = result.stdout
        duration = None
        for line in output.splitlines():
            if "Time for CLUSTER2:" in line:
                duration = float(line.split(": ")[1])
                print(duration)
        if duration is None:
            duration = time.time() - start_time
        return duration
    except subprocess.CalledProcessError as e:
        print(f"Error in executing command: {command}\nError: {str(e)}")
        return time.time() - start_time


def make_dyndom_command_file(pair):
    pdb1 = convert_pdb_id(pair["p_pdb"])
    pdb2 = convert_pdb_id(pair["q_pdb"])
    title = f"{pdb1}_{pdb2}.w5"
    filename1 = f"{pdb1}.pdb"
    filename2 = f"{pdb2}.pdb"
    chain1id = pdb1.split("_")[-1]
    chain2id = pdb2.split("_")[-1]
    with open(f"{title}.command", "w") as file:
        file.write(f"title={title}\n")
        file.write("\n")
        file.write(f"filename1={filename1}\n")
        file.write(f"chain1id={chain1id}\n")
        file.write(f"filename2={filename2}\n")
        file.write(f"chain2id={chain2id}\n")
        file.write("window=5\n")
        file.write("domain=20\n")
        file.write("ratio=1.0\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="DynDom の実行時間を計測し，結果を CSV に保存するスクリプト"
    )
    parser.add_argument(
        "--csv_file",
        type=str,
        default="pdb_with_hinges.csv",
        help="入力 CSV ファイル名 (default: pdb_with_hinges.csv)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="dyndom_execution_time_dyn_improved.csv",
        help="出力 CSV ファイル名 (default: dyndom_execution_time_dyn_improved.csv)",
    )
    args = parser.parse_args()
    csv_file = args.csv_file
    output_csv = args.output
    main(csv_file=csv_file, output_csv=output_csv)
