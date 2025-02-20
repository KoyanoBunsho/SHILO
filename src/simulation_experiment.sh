#!/bin/bash

g++ simulation_sh_ilo.cpp -o simulation_sh_ilo -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ simulation_sh.cpp -o simulation_sh -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

sigma_val=(0.5 1.0 1.5 2.0 2.5 3.0)
model_type=("shlo" "sh" "shibuya")

for sigma in "${sigma_val[@]}"; do
    for hinge_num in {2..10}; do
        for model in "${model_type[@]}"; do
            ./simulation_sh "$hinge_num" "$model" "$sigma"
        done
    done
done

for sigma in "${sigma_val[@]}"; do
    for hinge_num in {2..10}; do
        for model in "${model_type[@]}"; do
            ./simulation_sh_ilo "$hinge_num" "shilo" "$sigma"
        done
    done
done

