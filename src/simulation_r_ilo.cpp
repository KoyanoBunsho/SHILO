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

std::string getToken(const std::string &str, size_t index) {
  size_t start = 0, end = 0, count = 0;
  while ((end = str.find('_', start)) != std::string::npos) {
    if (count == index) {
      return str.substr(start, end - start);
    }
    start = end + 1;
    count++;
  }
  if (count == index) {
    return str.substr(start);
  }
  return "";
}

namespace fs = std::experimental::filesystem;
std::string extractHingeIndices(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "ファイルが開けません: " << filename << std::endl;
    return "";
  }
  std::string line;
  std::getline(file, line);
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> row;
    while (std::getline(ss, cell, ',')) {
      row.push_back(cell);
    }
    if (row.size() > 5) {
      file.close();
      return row[4];
    }
  }
  file.close();
  return "";
}

int main(int argc, char **argv) {
  std::string save_method_name;
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <hinge_num> <method> <sigma>"
              << std::endl;
    return 1;
  }
  int hinge_num = std::stoi(argv[1]);
  std::ofstream myfile;
  std::string simulation_data_path = "simulation_data/";
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
      std::string p_path = entry.path().string();
      std::string q_path = "simulation_data/pdb" + pdb_id + "_" + chain_id +
                           "_hinge_" + std::to_string(hinge_num) + "_sigma" +
                           sigma + ".pdb";
      std::string hinge_path =
          "simulation_data_info/pdb" + pdb_id + "_" + chain_id + "_hinge_" +
          std::to_string(hinge_num) + "_sigma" + sigma + ".csv";

      if (fs::exists(q_path) && fs::exists(hinge_path)) {
        file_triples.push_back(std::make_tuple(p_path, q_path, hinge_path));
      }
    }
  }
  std::chrono::duration<double, std::milli> exec_time_ms;
#pragma omp parallel for
  for (int i = 0; i < (int)file_triples.size(); i++) {
    const auto &triple = file_triples[i];
    std::string p_pdb_id =
        std::get<0>(triple).substr(std::get<0>(triple).find_last_of("/") + 1);
    std::string p_chain_id = getToken(p_pdb_id, 1);
    PDBReader reader1(std::get<0>(triple));
    PDBReader reader2(std::get<1>(triple));
    std::set<int> res_numbers1 = reader1.get_residue_numbers(p_chain_id);
    std::set<int> res_numbers2 = reader2.get_residue_numbers(p_chain_id);
    std::set<int> intersected_res_numbers =
        intersect_residue_numbers(res_numbers1, res_numbers2);
    std::vector<std::tuple<double, double, double>> coordinatesP =
        reader1.get_CA_coordinates(intersected_res_numbers, p_chain_id);
    std::vector<std::tuple<double, double, double>> coordinatesQ =
        reader2.get_CA_coordinates(intersected_res_numbers, p_chain_id);
    Eigen::MatrixXd p, q;
    p = convert_to_matrix(coordinatesP).transpose();
    q = convert_to_matrix(coordinatesQ).transpose();
    std::string hinge_file = std::get<2>(triple);
    std::string hingeIndices = extractHingeIndices(hinge_file);
    int total_residue_length = p.cols();
    myfile << p_pdb_id << "," << total_residue_length << "," << hinge_num << ","
           << hingeIndices << ",";
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
    for (int iter = 0; iter < iter_num; iter++) {
      auto start = std::chrono::high_resolution_clock::now();
      std::vector<int> random_hinges =
          selectRandomHinges(total_residue_length, hinge_num);
      ProteinRMSDhinge rmsdh_calculator(p, q, 100);
      AblationResult res =
          rmsdh_calculator.RMSDhkPostProcessingLoop(random_hinges, hinge_num);
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
#pragma omp critical
      {
        if (iter < iter_num - 1) {
          myfile << hinge_index << "," << rmsdhk.rmsdh_result << ","
                 << exec_time_ms.count() << "," << res.iter_num << ",";
        } else {
          myfile << hinge_index << "," << rmsdhk.rmsdh_result << ","
                 << exec_time_ms.count() << "," << res.iter_num << std::endl;
        }
      }
    }
  }
  myfile.close();
}
