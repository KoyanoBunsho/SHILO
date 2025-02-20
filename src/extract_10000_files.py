import os
import gzip
import random
import pandas as pd
from glob import glob


def main():
    pdb_files = glob("all_pdb/*.ent.gz")
    sampled_files = random.sample(pdb_files, min(10000, len(pdb_files)))
    csv_filename = "selected_files.csv"
    df = pd.DataFrame(sampled_files, columns=["FilePath"])
    df.to_csv(csv_filename, index=False, encoding="utf-8", header=None)
    print(f"サンプリング結果を {csv_filename} に保存しました。")


if __name__ == "__main__":
    main()
