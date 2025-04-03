import json
import csv
import itertools

# res.json を読み込み
with open("res.json", "r", encoding="utf-8") as f:
    data = json.load(f)

# CSVファイルとして出力 (output.csv)
with open("output.csv", "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    # ヘッダー行を書き込む
    writer.writerow(["PED_id", "ensemble_id_p", "ensemble_id_q"])

    # 各PEDエントリについて処理
    for entry in data.get("result", []):
        ped_id = entry.get("entry_id")
        ensembles = entry.get("ensembles", [])
        # ensemble id のリストを抽出
        ensemble_ids = [ens.get("ensemble_id") for ens in ensembles]
        # 2個以上のensembleがある場合のみ処理
        if len(ensemble_ids) < 2:
            continue

        # 組み合わせ（順序の入れ替えは同一とみなす）を生成
        for id_p, id_q in itertools.combinations(ensemble_ids, 2):
            writer.writerow([ped_id, id_p, id_q])
