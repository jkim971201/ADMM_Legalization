#include "skygui/SkyGui.h"
#include "MainWindow.h"
#include "GuiConfig.h"
#include "LayoutScene.h"

#include <QApplication>

extern int    cmd_argc;
extern char** cmd_argv; 

namespace gui
{

SkyGui::SkyGui(std::shared_ptr<dbDatabase> db) 
  : db_(db) {}

SkyGui::~SkyGui() {}

void
SkyGui::display()
{
  // Must constrcut QApplication before a QWidget
  QApplication app(cmd_argc, cmd_argv);
  gui_config_ = std::make_shared<GuiConfig>(db_);
  MainWindow window;
  window.setDatabase(db_);
  window.setGuiConfig(gui_config_);
  window.init();
  window.show();
  int exit_code = app.exec();
}

}
