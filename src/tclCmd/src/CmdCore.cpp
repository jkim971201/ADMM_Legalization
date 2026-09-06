#include <cstdlib>
#include <cassert>
#include <tcl.h>

#include "tclCmd/CmdCore.h"

#include "tclCmd/CmdReadLef.hpp"
#include "tclCmd/CmdReadDef.hpp"
#include "tclCmd/CmdReadVerilog.hpp"
#include "tclCmd/CmdReadBookshelf.hpp"
#include "tclCmd/CmdSetTopModule.hpp"
#include "tclCmd/CmdWriteBookshelf.hpp"
#include "tclCmd/CmdWriteDef.hpp"
#include "tclCmd/CmdGlobalPlace.hpp"
#include "tclCmd/CmdLegalize.hpp"
#include "tclCmd/CmdDisplay.hpp"

namespace skyline
{

// This is a callback function when the tcl interpreter meets a command.
int cmdCbk(ClientData clientData, Tcl_Interp* interp, int objc, struct Tcl_Obj* const* objv)
{
  const char* cmd_name = Tcl_GetString(objv[0]);
  skyline::TclCmd* cmd_ptr = TclCmdList::getCmdByName(cmd_name);

  if(cmd_ptr == nullptr)
  {
    printf("Unknown Command : %s\n", cmd_name);
    return TCL_ERROR;
  }

  bool next_is_option_argument = false;
  TclCmdOption* cur_option = nullptr;
  for(int cnt = 1; cnt < objc; cnt++) 
  {
    struct Tcl_Obj* obj  = objv[cnt];
    const char* obj_char = Tcl_GetString(obj);
    std::string obj_str  = std::string(obj_char);

    if(obj_char[0] == '-')
    {
      cur_option = cmd_ptr->getOptionByName(obj_str);
      if(cur_option == nullptr)
      {
        printf("Unknown Option : %s\n", obj_char);
        return TCL_ERROR;
      }
      else if(next_is_option_argument == true)
      {
        printf("Syntax error in option argument\n");
        return TCL_ERROR;
      }
      else
      {
        cur_option->setValid();
        int num_option_arg = cur_option->getNumArg();
        assert(num_option_arg < 2); 
        if(num_option_arg > 0)
          next_is_option_argument = true;
        // does not support multiple argument. will be fixed.
      }
    }
    else
    {
      if(next_is_option_argument == true)
      {
        if(cur_option == nullptr)
        {
          printf("Give option to %s\n", cur_option->getName().c_str());
          return TCL_ERROR;
        }

        cur_option->addArgument(obj_str);
        cur_option = nullptr;
        next_is_option_argument = false;
      }
      else
      {
        // In this case, object is the argument of this command.
        cmd_ptr->giveArgument(obj_str);
      }
    }
  }

  bool cmd_validity = cmd_ptr->checkValid();
  if(cmd_validity == false)
  {
    printf("Syntax error in command %s\n", cmd_ptr->name().c_str());
    return TCL_ERROR;
  }

  cmd_ptr->execute();
  return TCL_OK;
}

// Static Class Implementation
std::unordered_map<std::string, std::unique_ptr<TclCmd>> TclCmdList::name2Cmd_;

void
TclCmdList::addTclCmd(Tcl_Interp* interp, std::unique_ptr<TclCmd> newCmd)
{
  Tcl_CreateObjCommand(interp, newCmd->name().c_str(), cmdCbk, nullptr, nullptr);
  // We don't use Client Data and Delete Callback function
  name2Cmd_[newCmd->name()] = std::move(newCmd);
}

TclCmd*
TclCmdList::getCmdByName(const char* name)
{
  std::string nameStr = std::string(name);
  auto itr = name2Cmd_.find(nameStr);

  if(itr != name2Cmd_.end())
    return itr->second.get(); // return raw pointer
  else
    return nullptr;
}

bool
TclCmd::checkValid() const
{
  bool is_valid = true;
  for(const auto& [name, option] : name_to_option_)
  {
    if(option->isValid() == true)
    {
      int num_opt_arg = option->getNumArg();
      if(num_opt_arg != option->getArguments().size())
      {
        is_valid = false;
        break;
      }
    }
  }

  return is_valid;
}

void
CmdCore::initTclCmds(Tcl_Interp* interp)
{
  // DB-related API
  TclCmdList::addTclCmd(interp, std::make_unique<CmdReadLef>("read_lef"));
  TclCmdList::addTclCmd(interp, std::make_unique<CmdReadDef>("read_def"));
  TclCmdList::addTclCmd(interp, std::make_unique<CmdReadVerilog>("read_verilog"));
  TclCmdList::addTclCmd(interp, std::make_unique<CmdReadBookshelf>("read_bookshelf"));
  TclCmdList::addTclCmd(interp, std::make_unique<CmdSetTopModule>("set_top_module"));

  TclCmdList::addTclCmd(interp, std::make_unique<CmdWriteDef>("write_def"));
  TclCmdList::addTclCmd(interp, std::make_unique<CmdWriteBookshelf>("write_bookshelf"));

  // GUI-related API
  TclCmdList::addTclCmd(interp, std::make_unique<CmdDisplay>("display"));

  // Engine-related API
  TclCmdList::addTclCmd(interp, std::make_unique<CmdGlobalPlace>("global_place"));
  TclCmdList::addTclCmd(interp, std::make_unique<CmdLegalize>("legalize"));
}

} // namespace skyline
