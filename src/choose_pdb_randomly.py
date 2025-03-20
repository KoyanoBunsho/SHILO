from glob import glob
import random
import csv
import os
from biopandas.pdb import PandasPdb
from joblib import Parallel, delayed
import os
from tqdm import tqdm


def main():
    pdb_files = glob("all_pdb/*.ent.gz")
    random.seed(62)
    n_jobs = os.cpu_count()
    results = Parallel(n_jobs=n_jobs)(
        delayed(has_ca_atom)(pdb_file) for pdb_file in tqdm(pdb_files)
    )
    ca_files = [result for result in results if result is not None]

    selected_files = random.sample(ca_files, min(10000, len(ca_files)))
    with open("selected_files.csv", "w", newline="") as csvfile:
        filewriter = csv.writer(csvfile)
        for filename in selected_files:
            filewriter.writerow([filename])


def has_ca_atom(pdb_file):
    try:
        ppdb = PandasPdb().read_pdb(pdb_file)
        if (ppdb.df["ATOM"]["atom_name"] == "CA").any():
            return pdb_file
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
    return None


if __name__ == "__main__":
    main()
