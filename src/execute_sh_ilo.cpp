#if defined(__GNUC__) && !defined(__clang__) && defined(__x86_64__)
#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")
#endif

#include "rmsdh_new.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

std::map<int, int> map_hinge_to_residue(const std::set<int> &residue_numbers) {
  std::map<int, int> hinge_to_residue_map;
  std::vector<int> hinge_indices(residue_numbers.size());
  std::iota(hinge_indices.begin(), hinge_indices.end(), 1);
  std::vector<int> sorted_residues(residue_numbers.begin(),
                                   residue_numbers.end());
  if (sorted_residues.size() != hinge_indices.size()) {
    std::cerr << "Error: Residue numbers and hinge indices size mismatch."
              << std::endl;
    return hinge_to_residue_map;
  }
  for (size_t i = 0; i < hinge_indices.size(); ++i) {
    hinge_to_residue_map[hinge_indices[i]] = sorted_residues[i];
  }

  return hinge_to_residue_map;
}

int main(int argc, char *argv[]) try {
  if (argc != 6) {
    std::cout << "Usage: " << argv[0]
              << " <pdb filename 1> <pdb filename 2> <chain_id 1> <chain id 2> "
                 "<hinge num>\n";
    return 1;
  }
  std::string pdb_filename1 = argv[1];
  std::string pdb_filename2 = argv[2];
  std::string chain_id1 = argv[3];
  std::string chain_id2 = argv[4];
  const std::string hinge_argument = argv[5];
  if (hinge_argument.empty() ||
      hinge_argument.find_first_not_of("0123456789") != std::string::npos) {
    throw std::invalid_argument("Hinge count must be a positive integer");
  }
  int hinge_num = std::stoi(hinge_argument);
  PDBReader reader1(pdb_filename1);
  PDBReader reader2(pdb_filename2);
  std::set<int> res_numbers1 = reader1.get_residue_numbers(chain_id1);
  std::set<int> res_numbers2 = reader2.get_residue_numbers(chain_id2);
  std::set<int> intersected_res_numbers =
      intersect_residue_numbers(res_numbers1, res_numbers2);
  std::vector<std::tuple<double, double, double>> coordinatesP =
      reader1.get_CA_coordinates(intersected_res_numbers, chain_id1);
  std::vector<std::tuple<double, double, double>> coordinatesQ =
      reader2.get_CA_coordinates(intersected_res_numbers, chain_id2);
  if (coordinatesP.empty() || coordinatesP.size() != coordinatesQ.size()) {
    throw std::invalid_argument("Structures must have matching CA coordinates");
  }
  if (hinge_num < 1 || static_cast<size_t>(hinge_num) >= coordinatesP.size()) {
    throw std::invalid_argument("Hinge count must be positive and less than the residue count");
  }
  Eigen::MatrixXd p, q;
  p = convert_to_matrix(coordinatesP).transpose();
  q = convert_to_matrix(coordinatesQ).transpose();
  int total_residue_length = p.cols();
  std::cout << "Residue number: " << total_residue_length << std::endl;
  std::vector<double> default_weights;
  for (int i = 0; i < total_residue_length; i++) {
    default_weights.push_back(1.0);
  }
  ConformationPair PQ_pair = MoveToOrigin(p, q, default_weights);
  double rmsd_result = CalcRMSD(PQ_pair.P, PQ_pair.Q, default_weights);
  std::cout << "RMSD: " << rmsd_result << std::endl;
  ProteinRMSDhinge rmsdh_hinge(PQ_pair.P, PQ_pair.Q, hinge_num);
  AblationResult rmsdh_hinge_cnt_result = rmsdh_hinge.CalcFastRMSDhKLoop();
  double rmsdh_result = rmsdh_hinge_cnt_result.rmsdh_result;
  std::vector<int> hinge_index_vec = rmsdh_hinge_cnt_result.hinge_index_vec;
  std::cout << "RMSDh: " << rmsdh_result << std::endl;
  std::string hinge_index = "";
  std::map<int, int> hinge_to_residue_map =
      map_hinge_to_residue(intersected_res_numbers);
  for (int i = 0; i < hinge_num; i++) {
    if (i < hinge_num - 1) {
      hinge_index +=
          (std::to_string(hinge_to_residue_map[hinge_index_vec[i]]) + " : ");
    } else {
      hinge_index += (std::to_string(hinge_to_residue_map[hinge_index_vec[i]]));
    }
  }
  std::cout << "Hinge index is: " << hinge_index << std::endl;
  return 0;
} catch (const std::exception &error) {
  std::cerr << "Error: " << error.what() << std::endl;
  return 1;
}
