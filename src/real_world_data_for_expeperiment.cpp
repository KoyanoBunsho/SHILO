#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")

#include "rmsdh_new.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <input_csv_path> <output_csv_path>"
              << std::endl;
    return 1;
  }
  std::string input_coord_path = argv[1];
  std::string input_csv_path = argv[2];
  std::string output_csv_path = argv[3];

  std::vector<std::vector<std::string>> pdb_chain_data;
  read_csv(pdb_chain_data, input_csv_path);

  std::vector<std::pair<std::string, std::string>> pdb_pair_vec;
  for (const auto &pdb_chain : pdb_chain_data) {
    pdb_pair_vec.push_back(std::make_pair(pdb_chain[0], pdb_chain[1]));
  }

  std::ofstream myfile(output_csv_path);
  if (!myfile) {
    std::cerr << "Error: Cannot open output file " << output_csv_path
              << std::endl;
    return 1;
  }
  myfile << "p_pdb_id,q_pdb_id" << std::endl;

  for (size_t i = 0; i < pdb_pair_vec.size(); i++) {
    std::string p_pdb_id, q_pdb_id;
    std::string p_pdb_chain_id = pdb_pair_vec[i].first;
    std::string q_pdb_chain_id = pdb_pair_vec[i].second;

    std::transform(p_pdb_chain_id.begin(), p_pdb_chain_id.begin() + 4,
                   std::back_inserter(p_pdb_id), ::tolower);
    std::transform(q_pdb_chain_id.begin(), q_pdb_chain_id.begin() + 4,
                   std::back_inserter(q_pdb_id), ::tolower);

    std::string p_chain_id = p_pdb_chain_id.substr(5, 1);
    std::string q_chain_id = q_pdb_chain_id.substr(5, 1);

    myfile << p_pdb_id << "_" << p_chain_id << "," << q_pdb_id << "_"
           << q_chain_id << std::endl;
  }
  myfile.close();
  return 0;
}
