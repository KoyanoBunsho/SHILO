#!/bin/bash

g++ specific_rmsdhk_dyndom_data.cpp -o specific_rmsdhk_dyndom_data -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk_dyndom.cpp -o specific_fast_rmsdhk_dyndom_data -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk_postpro_dyndom.cpp -o specific_fast_rmsdhk_postpro_dyndom -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk_postpro_dyndom_loop.cpp -o specific_fast_rmsdhk_postpro_dyndom_loop -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

g++ specific_rmsdhk.cpp -o specific_rmsdhk -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk.cpp -o specific_fast_rmsdhk -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk_postpro.cpp -o specific_fast_rmsdhk_postpro -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk_postpro_loop.cpp -o specific_fast_rmsdhk_postpro_loop -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

g++ specific_rmsdhk_more_data.cpp -o specific_rmsdhk_more_data -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk_more_data.cpp -o specific_fast_rmsdhk_more_data -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk_more_data_postpro.cpp -o specific_fast_rmsdhk_more_data_postpro -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ specific_fast_rmsdhk_more_data_postpro_loop.cpp -o specific_fast_rmsdhk_more_data_postpro_loop -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

for k in {5..10}
do
echo "-----$k-----"
echo "Shibuya 2008 dataset"
echo "-----Shibuya's method-----"
./specific_rmsdhk $k
echo "-----SH-----"
./specific_fast_rmsdhk $k
echo "-----SH + LO-----"
./specific_fast_rmsdhk_postpro $k
echo "-----SH + ILO-----"
./specific_fast_rmsdhk_postpro_loop $k

echo "-----PAR 2020 dataset-----"
echo "-----Shibuya's method-----"
./specific_rmsdhk_more_data $k
echo "-----SH-----"
./specific_fast_rmsdhk_more_data $k
echo "-----SH + LO-----"
./specific_fast_rmsdhk_more_data_postpro $k
echo "-----SH + ILO-----"
./specific_fast_rmsdhk_more_data_postpro_loop $k

echo "-----DynDom 2024 dataset-----"
echo "-----Shibuya's method-----"
./specific_rmsdhk_dyndom_data $k
echo "-----SH-----"
./specific_fast_rmsdhk_dyndom_data $k
echo "-----SH + LO-----"
./specific_fast_rmsdhk_postpro_dyndom $k
echo "-----SH + ILO-----"
./specific_fast_rmsdhk_postpro_dyndom_loop $k
done
