import pandas as pd
import os

expert_annotated_list = [
    "35",
    "25 : 55",
    "98 : 110",
    "13 : 81",
    "89 : 264",
    "135",
    "150 : 318",
    "86 : 177",
    "91 : 251",
    "91 : 193",
    "104 : 236",
    "36 : 63 : 102",
]

shibuya_method_list = [
    "35 : 52",
    "45 : 57",
    "98 : 110",
    "12 : 76",
    "84 : 253",
    "84 : 134",
    "42 : 141",
    "85 : 179",
    "92 : 251",
    "91 : 192",
    "103 : 235",
    "35 : 71",
]
flexprot_list = [
    "39",
    "",
    "",
    "",
    "89",
    "",
    "",
    "88 : 181",
    "85 : 245",
    "84 : 177",
    "101",
    "35 : 67",
]
fatcat_list = ["", "", "", "", "91", "", "", "89 : 183", "92 : 252", "81", "99", ""]
nahal_list = [
    "61",
    "48 : 59",
    "98 : 111",
    "14 : 83",
    "87 : 252",
    "72 : 136",
    "47 : 148",
    "79 : 179",
    "90 : 253",
    "90 : 162",
    "103 : 237",
    "36 : 69",
]


def main():
    os.makedirs("figures", exist_ok=True)
    small_fast_rmsdhk_2 = pd.read_csv("rmsdh_result/fast_rmsdh_hingek_cnt_2.csv")
    small_fast_rmsdhk_2_postpro = pd.read_csv(
        "rmsdh_result/fast_rmsdh_hingek_cnt_2_postpro.csv"
    )
    small_fast_rmsdhk_2_postpro_loop = pd.read_csv(
        "rmsdh_result/fast_rmsdh_hingek_cnt_2_postpro_loop.csv"
    )
    lsp_2_shibuya = pd.read_csv(
        "rmsdh_result/ablation_study_2_paper.csv"
    )
    lsp_2_shibuya_loop = pd.read_csv(
        "rmsdh_result/ablation_study_loop2.csv"
    )
    calc_ans_shibuya(small_fast_rmsdhk_2, small_fast_rmsdhk_2_postpro,small_fast_rmsdhk_2_postpro_loop)
    eval_lsp_shibuya(lsp_2_shibuya)
    eval_lsp_shibuya(lsp_2_shibuya_loop)


def calc_ans(detect_method_list, d=3):
    ans_list = []
    for i in range(len(expert_annotated_list)):
        exp = expert_annotated_list[i]
        detect = detect_method_list[i]
        true_hinge_indices = exp.split(" : ")
        detected_hinge_indices = detect.split(" : ")
        TP = 0
        FP = 0
        FN = 0
        if detected_hinge_indices != [""]:
            detected_ranges = [
                (int(label) - d, int(label) + d) for label in detected_hinge_indices
            ]
        else:
            detected_ranges = []
        for true in true_hinge_indices:
            if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
                TP += 1
            else:
                FN += 1
        for lower, upper in detected_ranges:
            if not any(lower <= int(true) <= upper for true in true_hinge_indices):
                FP += 1
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (
            (2 * precision * recall) / (precision + recall)
            if (precision + recall) != 0
            else 0
        )
        ans_list.append(
            {"precision": precision, "recall": recall, "f_measure": f_measure}
        )
    return ans_list


def calc_ans_par(exp, detect, d=3):
    true_hinge_indices = exp.split(" : ")
    detected_hinge_indices = detect.split(" : ")
    TP = 0
    FP = 0
    FN = 0
    if detected_hinge_indices != [""]:
        detected_ranges = [
            (int(label) - d, int(label) + d) for label in detected_hinge_indices
        ]
    else:
        detected_ranges = []
    for true in true_hinge_indices:
        if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
            TP += 1
        else:
            FN += 1
    for lower, upper in detected_ranges:
        if not any(lower <= int(true) <= upper for true in true_hinge_indices):
            FP += 1
    return {"TP": TP, "FP": FP, "FN": FN}


def calc_approximation_ratio_df(df, heuristic_df, more_monge_2):
    df = pd.DataFrame(
        {
            "tmr": more_monge_2["monotonicity_rate_1"].to_list(),
            "RMSDhkD": (heuristic_df["RMSDh"] - df["RMSDh"]).to_list(),
            "approximation_ratio": (heuristic_df["RMSDh"] / df["RMSDh"]).to_list(),
            "RMSDhk": df["RMSDh"].to_list(),
            "heuristic_RMSDhk": heuristic_df["RMSDh"].to_list(),
        },
        index=more_monge_2.index,
    )
    return df


def calc_ans_tp(detect_method_list, d=3):
    ans_list = []
    for i in range(len(expert_annotated_list)):
        exp = expert_annotated_list[i]
        detect = detect_method_list[i]
        true_hinge_indices = exp.split(" : ")
        detected_hinge_indices = detect.split(" : ")
        TP = 0
        FP = 0
        FN = 0
        if detected_hinge_indices != [""]:
            detected_ranges = [
                (int(label) - d, int(label) + d) for label in detected_hinge_indices
            ]
        else:
            detected_ranges = []
        for true in true_hinge_indices:
            if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
                TP += 1
            else:
                FN += 1
        for lower, upper in detected_ranges:
            if not any(lower <= int(true) <= upper for true in true_hinge_indices):
                FP += 1
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (
            (2 * precision * recall) / (precision + recall)
            if (precision + recall) != 0
            else 0
        )
        ans_list.append({"TP": TP, "FP": FP, "FN": FN})
    return ans_list


def calc_ans_shibuya(small_fast_rmsdhk_2, small_fast_rmsdhk_2_postpro, small_fast_rmsdhk_2_postpro_loop, d=3):
    sh_list = small_fast_rmsdhk_2["hinge_index"].tolist()
    sh_lo_list = small_fast_rmsdhk_2_postpro["hinge_index"].tolist()
    sh_ilo_list = small_fast_rmsdhk_2_postpro_loop["hinge_index"].tolist()
    method_f_measure_list = []
    for method, method_list in {
        "SH": sh_list,
        "SH + LO": sh_lo_list,
        "SH + ILO": sh_ilo_list,
        "shibuya": shibuya_method_list,
        "flexprot": flexprot_list,
        "fatcat": fatcat_list,
        "nahal": nahal_list,
    }.items():
        TP = 0
        FP = 0
        FN = 0
        ans_df = pd.DataFrame(calc_ans_tp(method_list, d=d))
        TP = ans_df["TP"].sum()
        FP = ans_df["FP"].sum()
        FN = ans_df["FN"].sum()
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (
            (2 * precision * recall) / (precision + recall)
            if (precision + recall) != 0
            else 0
        )
        method_f_measure_list.append(
            {
                "method": method,
                "d": d,
                "F-measure": f_measure,
                "Precision": precision,
                "Recall": recall,
            }
        )
        print(method, "%.3f" % f_measure, "&", "%.3f" % precision, "&", "%.3f" % recall)

def calc_ans_par(exp, detect, d=0):
    true_hinge_indices = exp.split(" : ")
    detected_hinge_indices = detect.split(" : ")
    TP = 0
    FP = 0
    FN = 0
    if detected_hinge_indices != [""]:
        detected_ranges = [
            (int(label) - d, int(label) + d) for label in detected_hinge_indices
        ]
    else:
        detected_ranges = []
    for true in true_hinge_indices:
        if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
            TP += 1
        else:
            FN += 1
    for lower, upper in detected_ranges:
        if not any(lower <= int(true) <= upper for true in true_hinge_indices):
            FP += 1
    return {"TP": TP, "FP": FP, "FN": FN}


def calc_acc_df(df, heuristic_df):
    f_measure_dict = {}
    for d in [0, 3]:
        acc = []
        for i in range(len(df)):
            acc.append(
                calc_ans_par(
                    df.loc[i]["hinge_index"], heuristic_df.loc[i]["hinge_index"], d
                )
            )
        acc_df = pd.DataFrame(acc)
        TP = acc_df["TP"].sum()
        FP = acc_df["FP"].sum()
        FN = acc_df["FN"].sum()
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (
            (2 * precision * recall) / (precision + recall)
            if (precision + recall) != 0
            else 0
        )
        f_measure_dict[f"F-measure_distance_{d}"] = f_measure
        f_measure_dict[f"Precision_distance_{d}"] = precision
        f_measure_dict[f"Recall_distance_{d}"] = recall
    return f_measure_dict


def eval_lsp_shibuya(lsp_2_shibuya):
    ans_f_measure = 0
    ans_precision = 0
    ans_recall = 0
    for i in range(100):
        acc_dict_2_lsp = calc_acc_df(pd.DataFrame({"hinge_index": expert_annotated_list}), lsp_2_shibuya[f"{i}"].reset_index().rename(columns={f"{i}": "hinge_index"}))
        ans_f_measure += acc_dict_2_lsp["F-measure_distance_3"]
        ans_recall += acc_dict_2_lsp["Recall_distance_3"]
        ans_precision += acc_dict_2_lsp["Precision_distance_3"]
    print("%.3f" % (ans_f_measure / 100), "&",
          "%.3f" % (ans_precision / 100), "&",
          "%.3f" % (ans_recall / 100))

if __name__ == "__main__":
    main()
