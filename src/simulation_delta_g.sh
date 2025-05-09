#! /usr/local/bin/nosh
#$ -S /usr/local/bin/nosh
#$ -cwd
#$ -l s_vmem=2.1G
#$ -pe mpi_32 128
g++ calc_delta_g.cpp -o calc_delta_g -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ calc_delta_g_par.cpp -o calc_delta_g_par -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ calc_delta_g_dyndom.cpp -o calc_delta_g_dyndom -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
#g++ calc_delta_g_simulation.cpp -o calc_delta_g_simulation -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

./calc_delta_g
./calc_delta_g_par
./calc_delta_g_dyndom
#sigma_val=(0.5)

#for sigma in "${sigma_val[@]}"; do
 #   for hinge_num in {2..5}; do
  #      ./calc_delta_g_simulation "$hinge_num" "$sigma"
   # done
#done
