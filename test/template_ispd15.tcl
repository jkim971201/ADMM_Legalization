set design_name __DESIGN__

read_lef ispd2015/${design_name}/tech.lef
read_lef ispd2015/${design_name}/cells_modified.lef
read_def ispd2015/${design_name}/gp_${design_name}.def

set tech_penalty 0.00
set coeff_to_init_rho 2.0

legalize \
  -tech_penalty ${tech_penalty} \
  -coeff_to_init_rho ${coeff_to_init_rho} \
  -max_iter 20000 \
  -min_iter_to_call_adaptive_search 100 \
  -coeff_to_neglect_displace 12.0 \
  -qp_refine \
  -use_adaptive_search \
  -sorting_algorithm "CUB"

write_def ispd2015/${design_name}/${design_name}_test.def
