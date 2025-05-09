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

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <hinge_num> <sigma>" << std::endl;
    return 1;
  }
  std::vector<std::vector<std::string>> pdb_chain_data;
  int hinge_num = std::stoi(argv[1]);
  std::string sigma = argv[2];
  std::ofstream myfile;
  std::string simulation_data_path = "simulation_data/";
  std::string simulation_data_info_path = "simulation_data_info/";
  myfile.open("monge_rate_simulation_" + std::to_string(hinge_num) + "_" +
              "sigma" + sigma + ".csv");
  myfile << "p_pdb_id,Residue length,actual_hinge_cnt,sigma,monge_rate_1"
         << std::endl;
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
  for (int i = 0; i < (int)file_triples.size(); i++) {
    const auto &triple = file_triples[i];
    std::string p_pdb_id =
        std::get<0>(triple).substr(std::get<0>(triple).find_last_of("/") + 1);
    Eigen::MatrixXd p = openMatrixData(std::get<0>(triple));
    Eigen::MatrixXd q = openMatrixData(std::get<1>(triple));
    if (p.cols() == 0 || q.cols() == 0) {
      std::cout << "No data" << std::endl;
      continue;
    }
    if (p.cols() != q.cols()) {
      std::cout << "p length: " << p.cols() << " q length: " << q.cols()
                << std::endl;
      std::cout << "The residue length is different" << std::endl;
      continue;
    }
    myfile << p_pdb_id;
    int total_residue_length = p.cols();
    std::vector<double> default_weights;
    for (int i = 0; i < total_residue_length; i++) {
      default_weights.push_back(1.0);
    }
    auto start = std::chrono::high_resolution_clock::now();
    ConformationPair PQ_pair = MoveToOrigin(p, q, default_weights);
    ProteinRMSDhinge rmsdh_calculator(PQ_pair.P, PQ_pair.Q, 100);
    double monge_rate2 = rmsdh_calculator.calcMongeRate2();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_time_ms = end - start;
    myfile << "," << total_residue_length << "," << hinge_num << "," << sigma
           << "," << monge_rate2 << std::endl;
  }
  myfile.close();
}
