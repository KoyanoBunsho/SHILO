#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh_new.h"

int main() {
  std::string coord_path = "coord_csv_dyndom/";
  std::ofstream myfile;
  std::vector<std::vector<std::string>> pdb_chain_data;
  read_csv(pdb_chain_data, "DynDom_database.csv");
  std::vector<std::pair<std::string, std::string>> pdb_pair_vec;
  for (const auto &pdb_chain : pdb_chain_data) {
    pdb_pair_vec.push_back(std::make_pair(pdb_chain[0], pdb_chain[1]));
  }
  myfile.open("monge_rate_dyndom.csv");
  myfile << "p_pdb_id"
         << ","
         << "q_pdb_id"
         << ","
         << "monge_rate,monge_rate_1" << std::endl;
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
    myfile << p_pdb_id << "," << q_pdb_id;
    int total_residue_length = p.cols();
    std::cout << total_residue_length << std::endl;
    std::cout << p_pdb_id << ", " << q_pdb_id << std::endl;
    std::vector<double> default_weights;
    for (int i = 0; i < total_residue_length; i++) {
      default_weights.push_back(1.0);
    }
    auto start = std::chrono::high_resolution_clock::now();
    ConformationPair PQ_pair = MoveToOrigin(p, q, default_weights);
    ProteinRMSDhinge rmsdh_calculator(PQ_pair.P, PQ_pair.Q, 100);
    double monge_rate = rmsdh_calculator.calcMongeRate();
    double monge_rate2 = rmsdh_calculator.calcMongeRate2();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_time_ms = end - start;
    std::cout << exec_time_ms.count() << " ms" << std::endl;
    myfile << "," << monge_rate << "," << monge_rate2 << std::endl;
  }
  myfile.close();
}
