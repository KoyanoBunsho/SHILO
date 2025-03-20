g++ ablation_study_postprocessing.cpp -O3 -o ablation_study_postprocessing
g++ ablation_study_postprocessing_loop.cpp -O3 -o ablation_study_postprocessing_loop
g++ ablation_study_postprocessing_par.cpp -O3 -o ablation_study_postprocessing_par
g++ ablation_study_postprocessing_par_loop.cpp -O3 -o ablation_study_postprocessing_par_loop
g++ ablation_study_postprocessing_dyndom.cpp -O3 -o ablation_study_postprocessing_dyndom
g++ ablation_study_postprocessing_dyndom_loop.cpp -O3 -o ablation_study_postprocessing_dyndom_loop


for k in {5..10}
do
echo "-----$k-----"
echo "-----Shibuya 2008 dataset-----"
echo "-----R + LO-----"
./ablation_study_postprocessing $k
echo "-----R + ILO-----"
./ablation_study_postprocessing_loop $k

echo "-----PAR 2020 dataset-----"
echo "-----R + LO-----"
./ablation_study_postprocessing_par $k
echo "-----R + ILO-----"
./ablation_study_postprocessing_par_loop $k

echo "-----DynDom 2024 dataset-----"
echo "-----R + LO-----"
./ablation_study_postprocessing_dyndom $k
echo "-----R + ILO-----"
./ablation_study_postprocessing_dyndom_loop $k
done
