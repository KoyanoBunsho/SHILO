#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

namespace fs = std::experimental::filesystem;

int main() {
  const std::string simulation_data_path = "simulation_data/";
  const std::string simulation_data_info_path = "simulation_data_info/";
  const std::string output_file_name = "output_simulation_file.csv";
  const std::regex original_pattern(".*_original\\.pdb");
  const std::regex hinge_pattern(".*_hinge_.*\\.pdb");
  try {
    std::ofstream output_file(output_file_name);
    output_file << "p_pdb,q_pdb,hinge_file_path\n";
    if (fs::exists(simulation_data_path) &&
        fs::is_directory(simulation_data_path)) {
      for (const auto &entry : fs::directory_iterator(simulation_data_path)) {
        const auto &path = entry.path();
        std::string filename = path.filename().string();

        if (std::regex_match(filename, original_pattern)) {
          std::string base_name =
              filename.substr(0, filename.find("_original"));
          std::string original_file_path = path.string();
          std::regex hinge_file_pattern(base_name + "_hinge_.*\\.pdb");
          for (const auto &hinge_entry :
               fs::directory_iterator(simulation_data_path)) {
            const auto &hinge_path = hinge_entry.path();
            std::string hinge_filename = hinge_path.filename().string();

            if (std::regex_match(hinge_filename, hinge_file_pattern)) {
              std::string hinge_indices_csv =
                  hinge_filename.substr(0, hinge_filename.find(".pdb")) +
                  ".csv";
              std::string hinge_file_csv =
                  hinge_indices_csv.substr(0, hinge_indices_csv.find(".csv")) +
                  ".pdb";
              output_file << original_file_path << ","
                          << "simulation_data/" + hinge_file_csv << ","
                          << simulation_data_info_path + hinge_indices_csv
                          << "\n";
            }
          }
        }
      }
      output_file.close();
      std::cout << "Output file created successfully: " << output_file_name
                << std::endl;
    } else {
      std::cerr << "Directory does not exist or is not a directory: "
                << simulation_data_path << std::endl;
    }
  } catch (const fs::filesystem_error &e) {
    std::cerr << "Filesystem error: " << e.what() << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Standard exception: " << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown error occurred." << std::endl;
  }
  return 0;
}
