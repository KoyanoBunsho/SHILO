#include "rmsdh_new.h"
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

std::vector<int> selectRandomHinges(int n, int k) {
  std::vector<int> hinges;
  for (int i = 1; i <= n; ++i) {
    hinges.push_back(i);
  }
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(hinges.begin(), hinges.end(), g);
  hinges.resize(k);
  std::sort(hinges.begin(), hinges.end());
  return hinges;
}

int main(int argc, char **argv) {
  if (argc <= 1) {
    std::cout << "Command line error" << std::endl;
    return 0;
  }
  int hinge_num = atoi(argv[1]);
  std::string coord_path = "coord_csv/";
  std::ofstream myfile;
  std::vector<std::vector<std::string>> pdb_chain_data;
  read_csv(pdb_chain_data, "pdb_with_hinges.csv");
  std::vector<std::pair<std::string, std::string>> pdb_pair_vec;
  for (const auto &pdb_chain : pdb_chain_data) {
    pdb_pair_vec.push_back(std::make_pair(pdb_chain[0], pdb_chain[1]));
  }
  myfile.open("rmsdh_result/ablation_study_par_" + std::to_string(hinge_num) +
              ".csv");
  const int iter_num = 100;
  myfile << "p_pdb_id,q_pdb_id,Residue length,";
  for (int i = 0; i < iter_num; i++) {
    if (i <= iter_num - 2) {
      myfile << i << "," << std::to_string(i) + "_RMSDhk,"
             << std::to_string(i) + "_computation_time,";
    } else {
      myfile << i << "," << std::to_string(i) + "_RMSDhk,"
             << std::to_string(i) + "_computation_time" << std::endl;
    }
  }
  std::chrono::duration<double, std::milli> exec_time_ms;
  for (int i = 0; i < (int)pdb_pair_vec.size(); i++) {
    std::string p_pdb_id, q_pdb_id;
    std::string p_pdb_chain_id = pdb_pair_vec[i].first;
    std::string q_pdb_chain_id = pdb_pair_vec[i].second;
    std::transform(p_pdb_chain_id.begin(), p_pdb_chain_id.begin() + 4,
                   std::back_inserter(p_pdb_id), ::tolower);
    std::string p_chain_id = p_pdb_chain_id.substr(5, 1);
    std::transform(q_pdb_chain_id.begin(), q_pdb_chain_id.begin() + 4,
                   std::back_inserter(q_pdb_id), ::tolower);
    std::string q_chain_id = q_pdb_chain_id.substr(5, 1);
    Eigen::MatrixXd p, q;
    p = openMatrixData(coord_path + "coord_" + p_pdb_id + "_" + p_chain_id +
                       "_" + q_pdb_id + "_" + q_chain_id + ".csv");
    q = openMatrixData(coord_path + "coord_" + q_pdb_id + "_" + q_chain_id +
                       "_" + p_pdb_id + "_" + p_chain_id + ".csv");
    int total_residue_length = p.cols();
    myfile << p_pdb_id << "," << q_pdb_id << "," << total_residue_length << ",";
    if (p.cols() != q.cols()) {
      std::cout << "p length: " << p.cols() << " q length: " << q.cols()
                << std::endl;
      std::cout << "The residue length is different" << std::endl;
      continue;
    }
    for (int iter = 0; iter < iter_num; iter++) {
      auto start = std::chrono::high_resolution_clock::now();
      std::vector<int> random_hinges =
          selectRandomHinges(total_residue_length, hinge_num);
      ProteinRMSDhinge rmsdh_calculator(p, q, 100);
      RMSDhHingeCnt res =
          rmsdh_calculator.RMSDhkPostProcessing(random_hinges, hinge_num);
      auto end = std::chrono::high_resolution_clock::now();
      exec_time_ms = end - start;
      RMSDhHingeCnt rmsdhk = rmsdh_calculator.calcRMSDhKAfterHingeUpdate(
          res.hinge_index_vec, hinge_num);
      std::string hinge_index = "";
      for (int i = 0; i < res.hinge_index_vec.size(); i++) {
        if (i < res.hinge_index_vec.size() - 1) {
          hinge_index += (std::to_string(res.hinge_index_vec[i]) + " : ");
        } else {
          hinge_index += (std::to_string(res.hinge_index_vec[i]));
        }
      }
      if (iter < iter_num - 1) {
        myfile << hinge_index << "," << rmsdhk.rmsdh_result << ","
               << exec_time_ms.count() << ",";
      } else {
        myfile << hinge_index << "," << rmsdhk.rmsdh_result << ","
               << exec_time_ms.count();
      }
    }
    myfile << std::endl;
  }
  return 0;
}
