#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")

#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh_new.h"
#include <algorithm>
#include <experimental/filesystem>
#include <iostream>
#include <map>
#include <numeric>
#include <omp.h>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>

namespace fs = std::experimental::filesystem;

int main(int argc, char **argv) {
  std::string save_method_name;
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <hinge_num> <method> <sigma>"
              << std::endl;
    return 1;
  }
  int hinge_num = std::stoi(argv[1]);
  std::ofstream myfile;
  std::string simulation_data_path = "coord_csv_simulation/";
  std::string simulation_data_info_path = "simulation_data_info/";
  std::string sigma = argv[2];
  std::string save_name = "rmsdh_result/simulation_r_ilo_" +
                          std::to_string(hinge_num) + "_" + sigma + ".csv";
  myfile.open(save_name);
  const int iter_num = 100;
  myfile << "p_pdb_id,Residue length,hinge_num,actual_hinge_indices,";
  for (int i = 0; i < iter_num; i++) {
    if (i <= iter_num - 2) {
      myfile << i << "," << std::to_string(i) + "_hinge_index"
             << std::to_string(i) + "_RMSDhk,"
             << std::to_string(i) + "_computation_time,";
    } else {
      myfile << i << "," << std::to_string(i) + "_hinge_index"
             << std::to_string(i) + "_RMSDhk,"
             << std::to_string(i) + "_computation_time" << std::endl;
    }
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
#pragma omp parallel for
  for (const auto &triple : file_triples) {
    std::string p_pdb_id =
        std::get<0>(triple).substr(std::get<0>(triple).find_last_of("/") + 1);
    Eigen::MatrixXd p = openMatrixData(std::get<0>(triple));
    Eigen::MatrixXd q = openMatrixData(std::get<1>(triple));
    std::string hinge_file = std::get<2>(triple);
    std::string hingeIndices = extractHingeIndices(hinge_file);
    int total_residue_length = p.cols();
    if (p.cols() != q.cols()) {
      std::cout << "p length: " << p.cols() << " q length: " << q.cols()
                << std::endl;
      std::cout << "The residue length is different" << std::endl;
      continue;
    }
    if (p.cols() == 0) {
      std::cout << "No data" << std::endl;
      continue;
    }
    std::cout << total_residue_length << std::endl;
    std::chrono::duration<double, std::milli> exec_time_ms;
#pragma omp critical
    {
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
        for (int i = 0; i < (int)res.hinge_index_vec.size(); i++) {
          if (i < (int)res.hinge_index_vec.size() - 1) {
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
                 << exec_time_ms.count() << std::endl;
        }
      }
    }
  }
  myfile.close();
}
