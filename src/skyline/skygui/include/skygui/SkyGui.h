#ifndef SKY_GUI_H
#define SKY_GUI_H

#include <memory>

#include "db/dbDatabase.h"

using namespace db;

namespace gui
{

class MainWindow;
class LayoutScene;
class GuiConfig;

class SkyGui
{
  public:

    SkyGui(std::shared_ptr<dbDatabase> db);
    ~SkyGui();

    void display();

  private:

    std::unique_ptr<MainWindow>   main_window_;
    std::shared_ptr<GuiConfig>    gui_config_;
    std::shared_ptr<dbDatabase>   db_;
};

}

#endif
