#!/bin/bash

g++ simulation_r_lo.cpp -o simulation_r_lo -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ simulation_r_ilo.cpp -o simulation_r_ilo -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

sigma_val=(0.5 1 1.5)

for sigma in "${sigma_val[@]}"; do
    for hinge_num in {2..10}; do
        ./simulation_r_lo "$hinge_num" "$sigma"
        ./simulation_r_ilo "$hinge_num" "$sigma"
    done
done

