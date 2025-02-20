#ifndef rmsd_io
#define rmsd_io
#include "rmsd_struct.h"
#include <set>

void read_csv(std::vector<std::vector<std::string>> &data,
              std::string csv_path) {
  std::ifstream ifs(csv_path);
  if (!ifs) {
    std::cerr << "Failed to open file." << std::endl;
    return;
  }
  // Skip the column name line
  std::string skipline;
  std::getline(ifs, skipline);

  std::string line;
  while (std::getline(ifs, line)) {
    std::vector<std::string> row;
    size_t pos = 0;
    std::string delimiter = ",";
    while ((pos = line.find(delimiter)) != std::string::npos) {
      std::string token = line.substr(0, pos);
      row.push_back(token);
      line.erase(0, pos + delimiter.length());
    }
    row.push_back(line);
    data.push_back(row);
  }
}

template <typename T>
bool getFileContent(std::string fileName, std::vector<T> &Flexibility) {
  std::ifstream in(fileName.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fileName << std::endl;
    return false;
  }
  T val;
  while (in >> val) {
    Flexibility.push_back(val);
  }
  in.close();
  return true;
}

Eigen::MatrixXd openMatrixData(std::string fileToOpen) {
  std::vector<double> matrixEntries;
  std::ifstream matrixDataFile(fileToOpen);
  if (!matrixDataFile) {
    return Eigen::MatrixXd();
  }
  std::string matrixRowString;
  std::string matrixEntry;
  int matrixRowNumber = 0;
  while (getline(matrixDataFile, matrixRowString)) {
    std::stringstream matrixRowStringStream(matrixRowString);
    while (getline(matrixRowStringStream, matrixEntry, ',')) {
      matrixEntries.push_back(stod(matrixEntry));
    }
    matrixRowNumber++;
  }
  if (matrixRowNumber > 0) {
    return Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        matrixEntries.data(), matrixRowNumber,
        matrixEntries.size() / matrixRowNumber);
  } else {
    return Eigen::MatrixXd();
  }
}
bool getFileStrContent(
    std::string fileName,
    std::vector<std::pair<std::string, std::string>> &pq_pair) {
  std::ifstream in(fileName.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fileName << std::endl;
    return false;
  }
  std::string val_p, val_q;
  while (in >> val_p >> val_q) {
    pq_pair.push_back(make_pair(val_p, val_q));
  }
  in.close();
  return true;
}
class PDBReader {
private:
  std::string pdb_filename;

public:
  PDBReader(const std::string &filename) : pdb_filename(filename) {}

  std::set<int> get_residue_numbers(const std::string &chain_id) const {
    std::ifstream file(pdb_filename);
    std::string line;
    std::set<int> residue_numbers;

    if (file.is_open()) {
      while (getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        std::getline(ss, field, ' ');
        if (field == "ATOM") {
          std::string atom_name;
          std::string chain;
          int res_num;
          ss >> field;
          ss >> atom_name;
          ss >> field;
          ss >> chain;
          ss >> res_num;
          if (chain == chain_id) {
            residue_numbers.insert(res_num);
          }
        }
      }
      file.close();
    } else {
      std::cout << "Unable to open file: " << pdb_filename << std::endl;
    }

    return residue_numbers;
  }

  std::vector<std::tuple<double, double, double>>
  get_CA_coordinates(const std::set<int> &residue_numbers,
                     const std::string &chain_id) const {
    std::ifstream file(pdb_filename);
    std::string line;
    std::vector<std::tuple<double, double, double>> coordinates;

    if (file.is_open()) {
      while (getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        std::getline(ss, field, ' ');
        if (field == "ATOM") {
          std::string atom_name;
          std::string chain;
          int res_num;
          double x, y, z;
          ss >> field;
          ss >> atom_name;
          ss >> field;
          ss >> chain;
          ss >> res_num;
          ss >> x >> y >> z;
          if (atom_name == "CA" && chain == chain_id &&
              residue_numbers.find(res_num) != residue_numbers.end()) {
            coordinates.push_back(std::make_tuple(x, y, z));
          }
        }
      }
      file.close();
    } else {
      std::cout << "Unable to open file: " << pdb_filename << std::endl;
    }

    return coordinates;
  }
};

std::set<int> intersect_residue_numbers(const std::set<int> &set1,
                                        const std::set<int> &set2) {
  std::set<int> intersection;
  std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                        std::inserter(intersection, intersection.begin()));
  return intersection;
}

Eigen::MatrixXd
convert_to_matrix(const std::vector<std::tuple<double, double, double>> &vec) {
  Eigen::MatrixXd matrix(vec.size(), 3);
  for (int i = 0; i < (int)vec.size(); ++i) {
    matrix(i, 0) = std::get<0>(vec[i]);
    matrix(i, 1) = std::get<1>(vec[i]);
    matrix(i, 2) = std::get<2>(vec[i]);
  }
  return matrix;
}

#endif
