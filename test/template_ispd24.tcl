set design_name __DESIGN__

read_lef /home/jkim/ContestBenchmarks/ispd2024/Nangate.lef
read_def ispd2024/${design_name}/${design_name}.gp.def

legalize -qp_refine 
