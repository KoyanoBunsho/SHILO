#ifndef rmsd_io
#define rmsd_io
#include "rmsd_struct.h"
#include "pdb_reader.h"
#include <algorithm>
#include <sstream>
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

#endif
