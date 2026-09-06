#ifndef CMD_CORE_H
#define CMD_CORE_H

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <unordered_map>

extern "C" {
struct Tcl_Interp;
}

namespace skyline
{

class TclCmdOption
{
  public:

    TclCmdOption(std::string name, int num_arg) 
      : name_(name), num_arg_(num_arg), is_valid_(false) {}

    ~TclCmdOption() {}

    void setValid() { is_valid_ = true; }
    void addArgument(std::string arg) { arguments_.push_back(arg); }

    const std::string& getName() { return name_; }

    int getNumArg() const { return num_arg_; }

    bool isValid() const { return is_valid_; }

    std::vector<std::string>& getArguments() { return arguments_; }

  private:

    std::string name_;
    int num_arg_;
    bool is_valid_;
    std::vector<std::string> arguments_;
};

class TclCmd
{
 public:

  TclCmd(const char* cmd_name) { name_ = std::string(cmd_name); }

  ~TclCmd() {}

  const std::string& name() { return name_; }

  void giveArgument(std::string& arg) { arg_ = arg; }

  TclCmdOption* getOptionByName(const std::string& opt_name) 
  {
    auto find_opt_itr = name_to_option_.find(opt_name);
    return find_opt_itr == name_to_option_.end() 
      ? nullptr : find_opt_itr->second.get();
  }

  bool checkValid() const;

  virtual void execute() = 0; // Pure Virtual Function

 protected:

  void addOption(std::string option_name, int num_arg)
  {
    std::vector<std::string> arguments(num_arg);
    std::unique_ptr<TclCmdOption> new_option 
      = std::make_unique<TclCmdOption>(option_name, num_arg);
    name_to_option_[option_name] = std::move(new_option);
  }

  std::string name_; // Command name
  std::string arg_;  // Argument of this command (do not confuse with argument of option)
  std::unordered_map<std::string, std::unique_ptr<TclCmdOption>> name_to_option_;
};

class TclCmdList
{
  public:
    static void addTclCmd(Tcl_Interp* interp, std::unique_ptr<TclCmd> newCmd);
    static TclCmd* getCmdByName(const char* name);

  private:
    static std::unordered_map<std::string, std::unique_ptr<TclCmd>> name2Cmd_;
};

class CmdCore
{
  public:
    // Tcl Interface
    static void initTclCmds(Tcl_Interp* interp);
};

}

#endif
