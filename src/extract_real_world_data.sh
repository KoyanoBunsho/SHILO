#!/bin/bash

g++ real_world_data_for_expeperiment.cpp -o real_world_data_for_expeperiment -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
./real_world_data_for_expeperiment "coord_csv/" "pdb_with_hinges.csv" "rmsdh_result/par_data_for_experiment.csv"
./real_world_data_for_expeperiment "coord_csv_dyndom/" "DynDom_database.csv" "rmsdh_result/dyn_data_for_experiment.csv"
