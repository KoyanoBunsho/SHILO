import argparse
import pandas as pd
from biopandas.pdb import PandasPdb


def main(csv_path):
    df = pd.read_csv(csv_path)
    filtered_df = df[df["execution_time"].isna()]
    ppdb = PandasPdb()
    for idx, row in filtered_df.iterrows():
        p_pdb = row["p_pdb"]
        q_pdb = row["q_pdb"]
        pdb_id = p_pdb.split("_")[0].lower()
        try:
            ppdb.fetch_pdb(pdb_id).to_pdb(f"pdb{pdb_id}.ent.gz")
            print(f"pdb{pdb_id}.ent.gz のダウンロードに成功しました。")
        except Exception as e:
            print(f"pdb{pdb_id}.ent.gz のダウンロードに失敗しました: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CSVからexecution_timeが空の行を抽出し、PDBファイルをダウンロードします。"
    )
    parser.add_argument(
        "csv_path",
        help="CSVファイルのパス (例: all_pdb/fatcat_execution_time_par_improved.csv)",
    )
    args = parser.parse_args()
    main(args.csv_path)
