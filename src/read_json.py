import json
import csv
import itertools


def main():
    with open("res.json", "r", encoding="utf-8") as f:
        data = json.load(f)

    with open("output.csv", "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PED_id", "ensemble_id_p", "ensemble_id_q"])
        for entry in data.get("result", []):
            ped_id = entry.get("entry_id")
            ensembles = entry.get("ensembles", [])
            ensemble_ids = [ens.get("ensemble_id") for ens in ensembles]
            if len(ensemble_ids) < 2:
                continue
            for id_p, id_q in itertools.combinations(ensemble_ids, 2):
                writer.writerow([ped_id, id_p, id_q])


if __name__ == "__main__":
    main()
