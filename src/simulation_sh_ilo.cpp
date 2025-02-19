#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")

#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh_new.h"
#include <experimental/filesystem>
#include <regex>
#include <sstream>

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
  int hinge_num;
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <hinge_num> <method> <sigma>"
              << std::endl;
    return 1;
  }
  std::ofstream myfile;
  std::string simulation_data_path = "simulation_data/";
  std::string simulation_data_info_path = "simulation_data_info/";
  if (std::string(argv[2]) == "sh_ilo") {
    save_method_name = "shilo";
  } else if (std::string(argv[2]) == "sh") {
    save_method_name = "sh";
  } else if (std::string(argv[2]) == "shibuya") {
    save_method_name = "shibuya";
  } else {
    std::cerr << "Invalid method name: " << argv[2] << std::endl;
    return 1;
  }
  double sigma = std::stod(argv[3]);
  // simulation_dataディレクトリの下にpdb{pdb_id}_{chain_id}_original.pdbというオリジナルのファイルとpdb{pdb_id}_{chain_id}_hinge_{hinge_num}_sigma{sigma}.pdbというファイルが存在する
  // simulation_dataディレクトリの下にあるファイル名を参考にして，3つのファイル名のpathを格納したtupleのvectorを作る
  // イメージ: (simulation_data/pdb1cbu_B_original.pdb,
  // simulation_data/pdb1cbu_B_sigma0.5.pdb,
  // simulation_data_info/pdb1cbu_B_hinge_1_sigma0.5.csv)
  // sigmaとhinge_numは標準入力から与えられる
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
  for (const auto &triple : file_triples) {
    std::string p_pdb_id =
        std::get<0>(triple).substr(std::get<0>(triple).find_last_of("/") + 1);
    PDBReader reader1(std::get<0>(triple));
    PDBReader reader2(std::get<1>(triple));
    std::set<int> res_numbers1 = reader1.get_residue_numbers(chain_id1);
    std::set<int> res_numbers2 = reader2.get_residue_numbers(chain_id2);
    std::set<int> intersected_res_numbers =
        intersect_residue_numbers(res_numbers1, res_numbers2);
    std::vector<std::tuple<double, double, double>> coordinatesP =
        reader1.get_CA_coordinates(intersected_res_numbers, chain_id1);
    std::vector<std::tuple<double, double, double>> coordinatesQ =
        reader2.get_CA_coordinates(intersected_res_numbers, chain_id2);
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
    std::cout << p_pdb_id << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    ProteinRMSDhinge rmsdh_hinge(PQ_pair.P, PQ_pair.Q, hinge_num);
    RMSDhHingeCnt rmsdh_hinge_cnt_result;
    if (save_method_name == "sh_ilo") {
      rmsdh_hinge_cnt_result = rmsdh_hinge.CalcFastRMSDhKLoop();
    } else if (save_method_name == "sh") {
      rmsdh_hinge_cnt_result = rmsdh_hinge.CalcFastRMSDhK();
    } else if (save_method_name == "shibuya") {
      rmsdh_hinge_cnt_result = rmsdh_hinge.CalcRMSDhK();
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_time_ms = end - start;
    double exec_time_s = exec_time_ms.count() / 1000.0;
    std::cout << exec_time_s << " s" << std::endl;
    double rmsdh_result = rmsdh_hinge_cnt_result.rmsdh_result;
    int hinge_cnt = rmsdh_hinge_cnt_result.hinge_cnt;
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
      myfile << p_pdb_id << ",";
      myfile << total_residue_length << ",";
      myfile << rmsd_result << ",";
      myfile << rmsdh_result << "," << hinge_num << "," << hingeIndices << ","
             << hinge_index << "," << sigma << "," << exec_time_s << std::endl;
    }
  }
  myfile.close();
}
