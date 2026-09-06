#ifndef CMD_GLOBAL_PLACE_H
#define CMD_GLOBAL_PLACE_H

#include "skyline/SkyLine.h"
#include "skyplace/SkyPlace.h"

#include "CmdCore.h"

namespace skyline
{
  class CmdGlobalPlace : public TclCmd
  {
    public:

      CmdGlobalPlace(const char* cmd_name) : TclCmd(cmd_name)
      {
        addOption("-target_density",     1);
        addOption("-target_overflow",    1);
        addOption("-initial_place",      1);
        addOption("-initial_place_only", 0);
        addOption("-initial_lambda",     1);
        addOption("-initial_gamma_inv",  1);
        addOption("-solver",             1);
        addOption("-plot",               1);
      }

      void execute() override
      {
        auto sky = skyline::SkyLine::getStaticPtr();
        auto skyplace_gp = sky->getSkyPlaceGP();

        auto opt_ovf = getOptionByName("-target_overflow");
        if(opt_ovf->isValid() == true)
        {
          const auto& arguments = opt_ovf->getArguments();
          double target_ovf = std::stod(arguments[0]);
          skyplace_gp->setTargetOverflow(target_ovf);
        }

        auto opt_den = getOptionByName("-target_density");
        if(opt_den->isValid() == true)
        {
          const auto& arguments = opt_den->getArguments();
          double target_den = std::stod(arguments[0]);
          skyplace_gp->setTargetDensity(target_den);
        }

        auto opt_ip = getOptionByName("-initial_place");
        if(opt_ip->isValid() == true)
        {
          const auto& arguments = opt_ip->getArguments();
          std::string ip_method = arguments[0];
          skyplace_gp->setInitialPlaceMethod(ip_method);
        }

        auto opt_ip_only = getOptionByName("-initial_place_only");
        if(opt_ip_only->isValid() == true)
          skyplace_gp->setInitialPlaceOnly();

        auto opt_init_lambda = getOptionByName("-initial_lambda");
        if(opt_init_lambda->isValid() == true)
        {
          const auto& arguments = opt_init_lambda->getArguments();
          float init_lambda = std::stof(arguments[0]);
          skyplace_gp->setInitLambda(init_lambda);
        }

        auto opt_init_gamma_inv = getOptionByName("-initial_gamma_inv");
        if(opt_init_gamma_inv->isValid() == true)
        {
          const auto& arguments = opt_init_gamma_inv->getArguments();
          float init_gamma_inv = std::stof(arguments[0]);
          skyplace_gp->setInitGammaInv(init_gamma_inv);
        }

        auto opt_solver = getOptionByName("-solver");
        if(opt_solver->isValid() == true)
        {
          const auto& arguments = opt_solver->getArguments();
          std::string solver_type = arguments[0];
          skyplace_gp->setSolver(solver_type);
        }

        auto opt_plot = getOptionByName("-plot");
        if(opt_plot->isValid() == true)
        {
          const auto& arguments = opt_plot->getArguments();
          std::string plot_path = arguments[0];
          skyplace_gp->setPlotMode(plot_path);
        }

        skyplace_gp->run();
      }
  };
}

#endif
