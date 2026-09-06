set design_name __DESIGN_NAME__

read_lef iccad2017/${design_name}/tech.lef
read_lef iccad2017/${design_name}/cells_modified.lef
read_def iccad2017/${design_name}/placed.def

legalize \
  -constraint iccad2017/${design_name}/placement.constraints \
  -coeff_to_init_rho 3.0 \
  -tech_penalty 1.95

write_def iccad2017/${design_name}/${design_name}_test.def
