#!/bin/bash
python make_simulation_data.py
g++ output_simulation_file.cpp -o output_simulation_file -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
./output_simulation_file


g++ simulation_sh_ilo.cpp -o simulation_sh_ilo -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ simulation_sh.cpp -o simulation_sh -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ simulation_r_lo.cpp -o simulation_r_lo -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ simulation_r_ilo.cpp -o simulation_r_ilo -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

sigma_val=(0.5 1.0 1.5)
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
        ./simulation_sh_ilo "$hinge_num" "shilo" "$sigma"
        ./simulation_r_lo "$hinge_num" "$sigma"
        ./simulation_r_ilo "$hinge_num" "$sigma"
    done
done

bash fatcat_dyndom_simulation.sh
bash fatcat_dyndom_result.sh
