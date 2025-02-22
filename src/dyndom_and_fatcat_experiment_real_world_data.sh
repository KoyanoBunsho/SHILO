# DynDom experiment
python dyndom_experiment_comp_time.py --csv_file pdb_12_with_hinges_lower.csv --output dyndom_execution_time_shibuya_improved.csv
python dyndom_experiment_comp_time.py --csv_file pdb_with_hinges.csv --output dyndom_execution_time_par_improved.csv
python dyndom_experiment_comp_time.py --csv_file DynDon_database.csv --output dyndom_execution_time_dyn_improved.csv

# FATCAT experiment
python fatcat_experiment_comp_time.py --csv_file pdb_12_with_hinges_lower.csv --output fatcat_execution_time_shibuya_improved.csv
python fatcat_experiment_comp_time.py --csv_file pdb_with_hinges.csv --output fatcat_execution_time_par_improved.csv
python fatcat_experiment_comp_time.py --csv_file DynDon_database.csv --output fatcat_execution_time_dyn_improved.csv
