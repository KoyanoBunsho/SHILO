#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")

#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh_new.h"
#include <algorithm>
#include <chrono>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <omp.h>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::experimental::filesystem;

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <hinge_num> <sigma>" << std::endl;
    return 1;
  }
  int hinge_num = std::stoi(argv[1]);
  std::string sigma = argv[2];
  std::string save_name = "rmsdh_result/simulation_r_lo_" +
                          std::to_string(hinge_num) + "_" + sigma + ".csv";
  std::ofstream myfile(save_name);
  const int iter_num = 100;

  myfile << "p_pdb_id,Residue length,k,actual_hinge_indices,";
  for (int i = 0; i < iter_num; i++) {
    if (i < iter_num - 1)
      myfile << std::to_string(i) + "_hinge_index,"
             << std::to_string(i) + "_RMSDhk,"
             << std::to_string(i) + "_computation_time,";
    else
      myfile << std::to_string(i) + "_hinge_index,"
             << std::to_string(i) + "_RMSDhk,"
             << std::to_string(i) + "_computation_time" << std::endl;
  }
  std::vector<std::tuple<std::string, std::string, std::string>> file_triples;
  for (const auto &entry : fs::directory_iterator("simulation_data")) {
    std::string filename = entry.path().filename().string();
    std::smatch match;
    std::regex pattern(R"(pdb([\w\d]+)_([A-Z])_original\.pdb)");

    if (std::regex_match(filename, match, pattern)) {
      std::string pdb_id = match[1].str();
      std::string chain_id = match[2].str();
      std::string p_path = "coord_csv_simulation/pdb" + pdb_id + "_" +
                           chain_id + "_original_CA_coordinates.csv";
      std::string q_path = "coord_csv_simulation/pdb" + pdb_id + "_" +
                           chain_id + "_hinge_" + std::to_string(hinge_num) +
                           "_sigma" + sigma + "_CA_coordinates.csv";
      std::string hinge_path =
          "simulation_data_info/pdb" + pdb_id + "_" + chain_id + "_hinge_" +
          std::to_string(hinge_num) + "_sigma" + sigma + ".csv";
      if (fs::exists(q_path) && fs::exists(hinge_path)) {
        file_triples.push_back(std::make_tuple(p_path, q_path, hinge_path));
      }
    }
  }

  std::vector<std::string> results(file_triples.size());

#pragma omp parallel for
  for (size_t i = 0; i < file_triples.size(); i++) {
    auto triple = file_triples[i];
    std::string p_path = std::get<0>(triple);
    std::string q_path = std::get<1>(triple);
    std::string hinge_path = std::get<2>(triple);

    std::string p_pdb_id = p_path.substr(p_path.find_last_of("/") + 1);

    std::ostringstream oss;

    Eigen::MatrixXd p = openMatrixData(p_path);
    Eigen::MatrixXd q = openMatrixData(q_path);
    int total_residue_length = p.cols();

    if (p.cols() != q.cols()) {
      std::cerr << "p length: " << p.cols() << " q length: " << q.cols()
                << std::endl;
      std::cerr << "The residue length is different" << std::endl;
      results[i] = "";
      continue;
    }
    if (p.cols() == 0) {
      std::cerr << "No data" << std::endl;
      results[i] = "";
      continue;
    }
    std::string hingeIndices = extractHingeIndices(hinge_path);

    oss << p_pdb_id << "," << total_residue_length << "," << hinge_num << ","
        << hingeIndices << ",";
    for (int iter = 0; iter < iter_num; iter++) {
      auto start = std::chrono::high_resolution_clock::now();
      std::vector<int> random_hinges =
          selectRandomHinges(total_residue_length, hinge_num);
      ProteinRMSDhinge rmsdh_calculator(p, q, 100);
      RMSDhHingeCnt res =
          rmsdh_calculator.RMSDhkPostProcessing(random_hinges, hinge_num);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> exec_time_ms = end - start;
      RMSDhHingeCnt rmsdhk = rmsdh_calculator.calcRMSDhKAfterHingeUpdate(
          res.hinge_index_vec, hinge_num);
      std::string hinge_index = "";
      for (size_t j = 0; j < res.hinge_index_vec.size(); j++) {
        hinge_index += std::to_string(res.hinge_index_vec[j]);
        if (j != res.hinge_index_vec.size() - 1) {
          hinge_index += " : ";
        }
      }

      if (iter < iter_num - 1)
        oss << hinge_index << "," << rmsdhk.rmsdh_result << ","
            << exec_time_ms.count() << ",";
      else
        oss << hinge_index << "," << rmsdhk.rmsdh_result << ","
            << exec_time_ms.count();
    }
    oss << "\n";
    results[i] = oss.str();
  }
  for (const auto &line : results) {
    if (!line.empty()) {
      myfile << line;
    }
  }
  myfile.close();
  return 0;
}
