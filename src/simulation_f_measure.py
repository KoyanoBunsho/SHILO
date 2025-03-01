import pandas as pd
import os
from tqdm import tqdm

def main():
    os.makedirs("figures", exist_ok=True)    
    simulation_sh_ilo_df = read_simulation_data("rmsdh_result/simulation_sh_ilo_combined.csv")
    simulation_sh_lo_df = read_simulation_data("rmsdh_result/simulation_sh_lo_combined.csv")
    simulation_sh_df = read_simulation_data("rmsdh_result/simulation_sh_combined.csv")
    simulation_r_ilo_df = read_simulation_data("rmsdh_result/simulation_r_ilo_combined.csv")
    simulation_r_lo_df = read_simulation_data("rmsdh_result/simulation_r_lo_combined.csv")
    simulation_shibuya_df = read_simulation_data("rmsdh_result/simulation_shibuya_combined.csv")
    # TODO: DynDomとFATCATの結果をまとめる
    df_dict = {}
    df_dict["R + LO"] = simulation_r_lo_df
    df_dict["R + ILO"] = simulation_r_ilo_df
    df_dict["SH"] = simulation_sh_df
    df_dict["SH + LO"] = simulation_sh_lo_df
    df_dict["SH + ILO"] = simulation_sh_ilo_df
    df_dict["Shibuya's method"] = simulation_shibuya_df
    latex_table = generate_latex_table_accuracy(df_dict)
    with open("figures/latex_table.tex", "w", encoding="utf-8") as f:
        f.write(latex_table)
    print("LaTeX table code written to figures/simulation_f_measure.tex")

def read_simulation_data(file_path):
    df = pd.read_csv(file_path).fillna("")
    df["actual_hinge_cnt"] = df["k"]
    return df

def calc_ans_dyndom(exp, detect, d=0):
    true_hinge_indices = exp.split(" : ")
    detected_hinge_indices = detect.split(" : ")
    TP = 0
    FP = 0
    FN = 0
    if detected_hinge_indices != [""]:
        detected_ranges = [(int(label) - d, int(label) + d) for label in detected_hinge_indices]
    else:
        detected_ranges = []
    for true in true_hinge_indices:
        if true != "":
            if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
                TP += 1
            else:
                FN += 1
        else:
            if detected_hinge_indices == [""]:
                TP += 1
    for lower, upper in detected_ranges:
        if true != "":
            if not any(lower <= int(true) <= upper for true in true_hinge_indices):
                FP += 1
        else:
            FP += 1
    return {"TP": TP, "FP": FP, "FN": FN}

def calc_acc_df(df, heuristic_df):
    df["hinge_index"] = df["hinge_index"].fillna("")
    heuristic_df["hinge_index"] = heuristic_df["hinge_index"].fillna("")
    f_measure_dict = {}
    for d in [3]:
        acc = []
        for i in range(len(df)):
            acc.append(
                calc_ans_dyndom(
                    df.loc[i]["hinge_index"],
                    heuristic_df.loc[i]["hinge_index"],
                    d
                )
            )
        acc_df = pd.DataFrame(acc)
        TP = acc_df["TP"].sum()
        FP = acc_df["FP"].sum()
        FN = acc_df["FN"].sum()
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (2 * precision * recall) / (precision + recall) if (precision + recall) != 0 else 0
        f_measure_dict[f"F-measure_distance_{d}"] = f_measure
        f_measure_dict[f"Precision_distance_{d}"] = precision
        f_measure_dict[f"Recall_distance_{d}"] = recall
        f_measure_dict[f"TP_{d}"] = TP
        f_measure_dict[f"FP_{d}"] = FP
        f_measure_dict[f"FN_{d}"] = FN
    return f_measure_dict


def calc_acc_df_multi(true_df, pred_df):
    metrics_list = []
    hinge_cols = [f"{i}_hinge_index" for i in range(100)]
    if not pred_df.columns.str.contains("0_hinge_index").any():
        return calc_acc_df(true_df, pred_df)
    for col in tqdm(hinge_cols, desc="Calculating multi-hinge metrics"):
        temp_pred = pred_df[[col]].copy().rename(columns={col: "hinge_index"})
        metrics = calc_acc_df(true_df.reset_index(drop=True), temp_pred.reset_index(drop=True))
        metrics_list.append(metrics)
    avg_metrics = {}
    keys = metrics_list[0].keys()
    for key in keys:
        avg_metrics[key] = sum(m[key] for m in metrics_list) / len(metrics_list)
    return avg_metrics


def generate_latex_table_accuracy(df_dict):
    results = []
    results_precision = []
    results_recall = []
    method_order = ["R + LO", "R + ILO", "SH", "SH + LO", "SH + ILO", "Shibuya's method"]
    for method in method_order:
        acc_df = df_dict[method]
        for cnt in range(2, 6):
            tmp_acc = acc_df[acc_df["actual_hinge_cnt"] == cnt]
            tmp_true = tmp_acc["actual_hinge_indices"].reset_index(drop=True).to_frame()
            tmp_true = tmp_true.rename(columns={"actual_hinge_indices": "hinge_index"})
            tmp_pred = tmp_acc.reset_index(drop=True)
            res_acc = calc_acc_df_multi(tmp_true, tmp_pred)
            f_measure = res_acc.get("F-measure_distance_3", 0)
            precision = res_acc.get("Precision_distance_3", 0)
            recall = res_acc.get("Recall_distance_3", 0)
            results.append({"Method": method, "k": cnt, "F-measure": f_measure})
            results_precision.append({"Method": method, "k": cnt, "Precision": precision})
            results_recall.append({"Method": method, "k": cnt, "Recall": recall})
    df_f = pd.DataFrame(results)
    df_p = pd.DataFrame(results_precision)
    df_r = pd.DataFrame(results_recall)
    df_merged = pd.merge(df_f, df_p, on=["Method", "k"])
    df_merged = pd.merge(df_merged, df_r, on=["Method", "k"])
    method_priority = {m: i for i, m in enumerate(method_order)}
    df_merged["method_order"] = df_merged["Method"].map(method_priority)
    df_merged = df_merged.sort_values(by=["k", "method_order"])
    lines = []
    lines.append(r"\begin{table*}[t]")
    lines.append(r"    \centering")
    lines.append(r"    \caption{Evaluation results (F-measure, Precision, and Recall) on the Simulation dataset for different assumed number of hinges $k$.}")
    lines.append(r"    \begin{tabular}{ccccc}")
    lines.append(r"        \hline")
    lines.append(r"        \multirow{2}{*}{$k$} & \multirow{2}{*}{Method} & \multicolumn{3}{c}{Simulation dataset} \\")
    lines.append(r"        \cline{3-5} \\")
    lines.append(r"         &  & F-measure & Precision & Recall \\")
    lines.append(r"        \hline")
    for k_value, group in df_merged.groupby("k"):
        group = group.reset_index(drop=True)
        n_rows = len(group)
        for idx, row in group.iterrows():
            if idx == 0:
                line = f"        \\multirow{{{n_rows}}}{{*}}{{{k_value}}} & {row['Method']} & {row['F-measure']:.3f} & {row['Precision']:.3f} & {row['Recall']:.3f} \\\\"
            else:
                line = f"         & {row['Method']} & {row['F-measure']:.3f} & {row['Precision']:.3f} & {row['Recall']:.3f} \\\\"
            lines.append(line)
        lines.append(r"        \hline")    
    lines.append(r"    \end{tabular}")
    lines.append(r"    \label{tab:simulation_dataset_results}")
    lines.append(r"\end{table*}")
    return "\n".join(lines)

if __name__ == "__main__":
    main()
