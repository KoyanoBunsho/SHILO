#! /usr/local/bin/nosh
#$ -S /usr/local/bin/nosh
#$ -cwd
#$ -l s_vmem=2.1G
#$ -pe mpi_32 128
g++ calc_monge_rate_simulation.cpp -o calc_monge_rate_simulation -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp


sigma_val=(0.5)

for sigma in "${sigma_val[@]}"; do
    for hinge_num in {2..10}; do
        ./calc_monge_rate_simulation "$hinge_num" "$sigma"
    done
done
