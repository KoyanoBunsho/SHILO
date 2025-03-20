import os
from glob import glob

import pandas as pd
from biopandas.pdb import PandasPdb
from joblib import Parallel, delayed
from tqdm import tqdm


def main():
    pdb_directory = "simulation_data"
    output_directory = "coord_csv_simulation"
    os.makedirs(output_directory, exist_ok=True)
    pdb_files = glob(os.path.join(pdb_directory, "*.pdb"))
    n_jobs = os.cpu_count() or 1
    print(f"Processing in {n_jobs} parallel jobs")

    Parallel(n_jobs=n_jobs)(
        delayed(process_pdb_file)(pdb_file, output_directory) for pdb_file in tqdm(pdb_files, total=len(pdb_files))
    )


def process_pdb_file(pdb_file, output_directory):
    pdb_id = os.path.basename(pdb_file).split(".pdb")[0]
    ppdb = PandasPdb().read_pdb(pdb_file)
    ca_atoms = ppdb.df["ATOM"][ppdb.df["ATOM"]["atom_name"] == "CA"].drop_duplicates(subset=["residue_number"])
    coordinates = ca_atoms[["x_coord", "y_coord", "z_coord"]].copy()

    if not coordinates.empty:
        output_file_path = os.path.join(
            output_directory, f"{pdb_id}_CA_coordinates.csv"
        )
        coordinates.T.to_csv(output_file_path, index=False, header=None)


if __name__ == "__main__":
    main()
