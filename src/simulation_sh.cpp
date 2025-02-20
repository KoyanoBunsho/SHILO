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
  int hinge_num = 2;
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <hinge_num> <method> <sigma>"
              << std::endl;
    return 1;
  }
  std::ofstream myfile;
  std::string simulation_data_path = "simulation_data/";
  std::string simulation_data_info_path = "simulation_data_info/";
  if (std::string(argv[2]) == "sh") {
    save_method_name = "shilo";
  } else if (std::string(argv[2]) == "shlo") {
    save_method_name = "sh_lo";
  } else if (std::string(argv[2]) == "shibuya") {
    save_method_name = "shibuya";
  } else {
    std::cerr << "Invalid method name: " << argv[2] << std::endl;
    return 1;
  }
  double sigma = std::stod(argv[3]);
  std::string save_name = "rmsdh_result/simulation_" + save_method_name + "_" +
                          std::to_string(hinge_num) + "_" +
                          std::to_string(sigma) + ".csv";
  myfile.open(save_name);
  myfile << "p_pdb_id,Residue "
            "length,RMSD,RMSDh,k,hinge_cnt,actual_hinge_"
            "indices,hinge_index,sigma,exec_time (s)"
         << std::endl;
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
                           "_sigma" + std::to_string(sigma) + ".pdb";
      std::string hinge_path =
          "simulation_data_info/pdb" + pdb_id + "_" + chain_id + "_hinge_" +
          std::to_string(hinge_num) + "_sigma" + std::to_string(sigma) + ".csv";

      if (fs::exists(q_path) && fs::exists(hinge_path)) {
        file_triples.push_back(std::make_tuple(p_path, q_path, hinge_path));
      }
    }
  }
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
    std::vector<double> default_weights;
    for (int i = 0; i < total_residue_length; i++) {
      default_weights.push_back(1.0);
    }
    ConformationPair PQ_pair = MoveToOrigin(p, q, default_weights);
    double rmsd_result = CalcRMSD(PQ_pair.P, PQ_pair.Q, default_weights);
    auto start = std::chrono::high_resolution_clock::now();
    ProteinRMSDhinge rmsdh_hinge(PQ_pair.P, PQ_pair.Q, hinge_num);
    RMSDhHingeCnt rmsdh_hinge_cnt_result = {};
    if (save_method_name == "sh") {
      rmsdh_hinge_cnt_result = rmsdh_hinge.CalcFastRMSDhK();
    } else if (save_method_name == "sh_lo") {
      rmsdh_hinge_cnt_result = rmsdh_hinge.CalcFastRMSDhK(true);
    } else if (save_method_name == "shibuya") {
      rmsdh_hinge_cnt_result = rmsdh_hinge.CalcRMSDhK();
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_time_ms = end - start;
    double exec_time_s = exec_time_ms.count() / 1000.0;
    std::cout << exec_time_s << " s" << std::endl;
    double rmsdh_final_result = rmsdh_hinge_cnt_result.rmsdh_result;
    std::vector<int> hinge_index_vec = rmsdh_hinge_cnt_result.hinge_index_vec;
    std::string hinge_index = "";
    for (int i = hinge_index_vec.size() - 1; i >= 0; i--) {
      if (i > 0) {
        hinge_index += (std::to_string(hinge_index_vec[i]) + " : ");
      } else {
        hinge_index += (std::to_string(hinge_index_vec[i]));
      }
    }
#pragma omp critical
    {
      myfile << p_pdb_id << "," << total_residue_length << "," << rmsd_result
             << "," << rmsdh_final_result << "," << hinge_num << ","
             << hingeIndices << "," << hinge_index << "," << sigma << ","
             << exec_time_s << std::endl;
    }
  }
  myfile.close();
}
