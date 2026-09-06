#ifndef CMD_LEGALIZE_H
#define CMD_LEGALIZE_H

#include "skyline/SkyLine.h"
#include "legalizer/Legalizer.h"
#include "CmdCore.h"

namespace skyline
{
  class CmdLegalize : public TclCmd
  {
    public:

      CmdLegalize(const char* cmd_name) : TclCmd(cmd_name) 
      {
        addOption("-constraint", 1);
        addOption("-tech_penalty", 1);
        addOption("-max_iter", 1);
        addOption("-std_admm", 0);
        addOption("-coeff_to_init_rho", 1);
        addOption("-coeff_to_init_dual", 1);

        addOption("-qp_refine", 0);
        addOption("-size_file", 1);

        addOption("-use_adaptive_search", 0);
        addOption("-min_iter_to_call_adaptive_search", 1);
        addOption("-ovf_ratio_to_call_adaptive_search", 1);
        addOption("-coeff_to_neglect_displace", 1);

        addOption("-random_sorting_seed", 1);
        addOption("-sorting_algorithm", 1);
        addOption("-sorting_policy", 1);
      }

      void execute() override
      {
        auto sky = skyline::SkyLine::getStaticPtr();
        auto skyplace_lg = sky->getSkyPlaceLG();

        auto opt_constraint = getOptionByName("-constraint");
        std::string file_path_constraint = std::string();
        if(opt_constraint->isValid() == true)
        {
          const auto& arguments = opt_constraint->getArguments();
          file_path_constraint = arguments[0];
          skyplace_lg->setICCAD17Constraint(file_path_constraint);
          skyplace_lg->setFlagICCAD2017(true);
        }

        auto opt_tech_penalty = getOptionByName("-tech_penalty");
        if(opt_tech_penalty->isValid() == true)
        {
          const auto& arguments = opt_tech_penalty->getArguments();
          double tech_penalty = std::stod(arguments[0]);
          skyplace_lg->setTechPenalty(tech_penalty);
        }

        auto opt_max_iter = getOptionByName("-max_iter");
        if(opt_max_iter->isValid() == true)
        {
          const auto& arguments = opt_max_iter->getArguments();
          int max_admm_iter = std::stoi(arguments[0]);
          skyplace_lg->setMaxAdmmIter(max_admm_iter);
        }

        auto opt_std_admm = getOptionByName("-std_admm");
        if(opt_std_admm->isValid() == true)
          skyplace_lg->setStandardAdmm(true);

        auto opt_coeff_to_init_rho = getOptionByName("-coeff_to_init_rho");
        if(opt_coeff_to_init_rho->isValid() == true)
        {
          const auto& arguments = opt_coeff_to_init_rho->getArguments();
          double coeff_to_init_rho = std::stod(arguments[0]);
          skyplace_lg->setCoeffToInitRho(coeff_to_init_rho);
        }

        auto opt_coeff_to_init_dual = getOptionByName("-coeff_to_init_dual");
        if(opt_coeff_to_init_dual->isValid() == true)
        {
          const auto& arguments = opt_coeff_to_init_dual->getArguments();
          double coeff_to_init_dual = std::stod(arguments[0]);
          skyplace_lg->setCoeffToInitDual(coeff_to_init_dual);
        }

        auto opt_use_adaptive_search = getOptionByName("-use_adaptive_search");
        if(opt_use_adaptive_search->isValid() == true)
          skyplace_lg->setAdaptiveSearch(true);

        auto opt_min_iter_to_call_adaptive_search = getOptionByName("-min_iter_to_call_adaptive_search");
        if(opt_min_iter_to_call_adaptive_search->isValid() == true)
        {
          const auto& arguments = opt_min_iter_to_call_adaptive_search->getArguments();
          int min_iter_to_call_adaptive_search = std::stoi(arguments[0]);
          skyplace_lg->setMinIterToCallAdaptiveSearch(min_iter_to_call_adaptive_search);
        }

        auto opt_ovf_ratio_to_call_adaptive_search = getOptionByName("-ovf_ratio_to_call_adaptive_search");
        if(opt_ovf_ratio_to_call_adaptive_search->isValid() == true)
        {
          const auto& arguments = opt_ovf_ratio_to_call_adaptive_search->getArguments();
          double ovf_ratio_to_call_adaptive_search = std::stod(arguments[0]);
          skyplace_lg->setOvfRatioToCallAdaptiveSearch(ovf_ratio_to_call_adaptive_search);
        }

        auto opt_coeff_to_neglect_displace = getOptionByName("-coeff_to_neglect_displace");
        if(opt_coeff_to_neglect_displace->isValid() == true)
        {
          const auto& arguments = opt_coeff_to_neglect_displace->getArguments();
          double coeff_to_neglect_displace = std::stod(arguments[0]);
          skyplace_lg->setCoeffToNeglectDisplace(coeff_to_neglect_displace);
        }

        auto opt_qp_refine = getOptionByName("-qp_refine");
        if(opt_qp_refine->isValid() == true)
          skyplace_lg->setQpRefine(true);

        auto opt_sorting_algorithm = getOptionByName("-sorting_algorithm");
        if(opt_sorting_algorithm->isValid() == true)
        {
          const auto& arguments = opt_sorting_algorithm->getArguments();
          std::string sorting_algo = std::string(arguments[0]);
          skyplace_lg->setSortingAlgorithm(sorting_algo);
        }
        
        auto opt_sorting_policy = getOptionByName("-sorting_policy");
        if(opt_sorting_policy->isValid() == true)
        {
          const auto& arguments = opt_sorting_policy->getArguments();
          std::string sorting_policy = std::string(arguments[0]);
          skyplace_lg->setSortingPolicy(sorting_policy);
        }

        auto opt_random_sorting_seed = getOptionByName("-random_sorting_seed");
        if(opt_random_sorting_seed->isValid() == true)
        {
          const auto& arguments = opt_random_sorting_seed->getArguments();
          int seed = std::stoi(arguments[0]);
          skyplace_lg->setRandomSortingSeed(seed);
        }

        auto opt_size_file = getOptionByName("-size_file");
        if(opt_size_file->isValid() == true)
        {
          const auto& arguments = opt_size_file->getArguments();
          std::string file_path_size = arguments[0];
          skyplace_lg->setSizeByFile(file_path_size);
        }

        skyplace_lg->run();
      }
  };
}

#endif
