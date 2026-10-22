#ifndef SHILO_PDB_READER_H
#define SHILO_PDB_READER_H

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

class PDBReader {
public:
  explicit PDBReader(const std::string &filename) : pdb_filename(filename) {}

  std::set<int> get_residue_numbers(const std::string &chain_id) const {
    std::set<int> residues;
    for (const auto &entry : read_ca_atoms(chain_id)) {
      residues.insert(entry.first);
    }
    return residues;
  }

  std::vector<std::tuple<double, double, double>>
  get_CA_coordinates(const std::set<int> &residue_numbers,
                     const std::string &chain_id) const {
    std::vector<std::tuple<double, double, double>> coordinates;
    const auto atoms = read_ca_atoms(chain_id);
    for (int residue : residue_numbers) {
      const auto atom = atoms.find(residue);
      if (atom == atoms.end()) {
        throw std::runtime_error("Missing CA coordinates in " + pdb_filename);
      }
      coordinates.push_back(atom->second);
    }
    return coordinates;
  }

private:
  std::string pdb_filename;

  std::map<int, std::tuple<double, double, double>>
  read_ca_atoms(const std::string &chain_id) const {
    if (chain_id.empty()) {
      throw std::invalid_argument("Chain ID must not be empty");
    }
    const auto dot = pdb_filename.find_last_of('.');
    const auto extension = dot == std::string::npos ? "" : pdb_filename.substr(dot);
    if (extension != ".pdb" && extension != ".ent" && extension != ".cif") {
      throw std::invalid_argument("Unsupported structure file: " + pdb_filename);
    }
    std::ifstream file(pdb_filename);
    if (!file) {
      throw std::runtime_error("Cannot open structure file: " + pdb_filename);
    }
    std::map<int, std::tuple<double, double, double>> atoms;
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream fields(line);
      std::string record, ignored, atom_name, chain;
      fields >> record;
      if (record != "ATOM") {
        continue;
      }
      int residue = 0;
      double x = 0.0, y = 0.0, z = 0.0;
      if (extension == ".cif") {
        fields >> ignored >> ignored >> atom_name >> ignored >> ignored >> chain
               >> ignored >> residue >> ignored >> x >> y >> z;
      } else {
        fields >> ignored >> atom_name >> ignored >> chain >> residue >> x >> y >> z;
      }
      if (!fields || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
        throw std::runtime_error("Invalid ATOM record in " + pdb_filename);
      }
      if (atom_name == "CA" && chain == chain_id) {
        // Keep the first alternate location and align both files by residue ID.
        atoms.emplace(residue, std::make_tuple(x, y, z));
      }
    }
    if (file.bad()) {
      throw std::runtime_error("Cannot read structure file: " + pdb_filename);
    }
    return atoms;
  }
};

#endif
