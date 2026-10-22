#include "rmsdh_new.h"
#include <functional>
#include <limits>
#include <stdexcept>

void require(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void rejects(const std::function<void()> &operation) {
  try {
    operation();
  } catch (const std::exception &) {
    return;
  }
  throw std::runtime_error("Invalid input was accepted");
}

int main(int argc, char **argv) {
  require(argc == 2, "Expected a temporary directory");
  const std::string directory = argv[1];
  const std::string pdb = directory + "/input.pdb";
  const std::string cif = directory + "/input.cif";
  const auto write_pdb = [&](const std::string &data) {
    std::ofstream(pdb) << data;
  };

  write_pdb("ATOM 1 CA GLY A 2 4 5 6\n"
            "ATOM 2 CA GLY A 1 1 2 3\n"
            "ATOM 3 CA GLY A 1 7 8 9\n"
            "ATOM 4 N GLY A 3 0 0 0\n"
            "ATOM 5 CA GLY B 4 0 0 0\n");
  PDBReader reader(pdb);
  const std::set<int> expected = {1, 2};
  require(reader.get_residue_numbers("A") == expected, "Only CA residues are aligned");
  const auto coordinates = reader.get_CA_coordinates(expected, "A");
  require(coordinates.size() == 2, "Alternate locations must not duplicate residues");
  require(coordinates[0] == std::make_tuple(1.0, 2.0, 3.0), "Residue order must be stable");
  rejects([&] { reader.get_CA_coordinates({3}, "A"); });
  rejects([&] { reader.get_residue_numbers(""); });
  rejects([&] { PDBReader(directory + "/missing.pdb").get_residue_numbers("A"); });
  rejects([&] { PDBReader(directory + "/input.txt").get_residue_numbers("A"); });

  for (const auto &record : {
           "ATOM 1 CA GLY A\n", "ATOM 1 CA GLY A 1 2\n",
           "ATOM 1 CA GLY A invalid 1 2 3\n",
           "ATOM 1 CA GLY A 1 nan 2 3\n",
           "ATOM 1 CA GLY A 1 1e999 2 3\n"}) {
    write_pdb(record);
    rejects([&] { reader.get_residue_numbers("A"); });
    rejects([&] { reader.get_CA_coordinates({1}, "A"); });
  }
  std::ofstream(cif) << "ATOM 1 C CA . GLY 1 1 1 ? 1 2 3\n";
  require(PDBReader(cif).get_residue_numbers("1") == std::set<int>{1},
          "Numeric CIF chains must be supported");
  std::ofstream(cif) << "ATOM 1 C CA . GLY 1 1 1 ? 1\n";
  rejects([&] { PDBReader(cif).get_CA_coordinates({1}, "1"); });

  const Eigen::Matrix3Xd points = Eigen::Matrix3Xd::Random(3, 5);
  const Eigen::Matrix3Xd empty(3, 0);
  const Eigen::Matrix3Xd shorter = points.leftCols(4);
  for (int hinges : {-1, 0, 5, std::numeric_limits<int>::max()}) {
    rejects([&] { ProteinRMSDhinge invalid(points, points, hinges); });
  }
  rejects([&] { ProteinRMSDhinge invalid(empty, empty, 1); });
  rejects([&] { ProteinRMSDhinge invalid(points, shorter, 1); });
  auto nonfinite = points;
  nonfinite(0, 0) = std::numeric_limits<double>::infinity();
  rejects([&] { ProteinRMSDhinge invalid(points, nonfinite, 1); });
  ProteinRMSDhinge valid(points, points, 1);
  require(std::isfinite(valid.CalcRMSDhK().rmsdh_result), "Valid exact RMSDh must be finite");
  require(std::isfinite(valid.CalcFastRMSDhKLoop().rmsdh_result), "Valid fast RMSDh must be finite");
  std::cout << "Input security regression tests passed" << std::endl;
}
