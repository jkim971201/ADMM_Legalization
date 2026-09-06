#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QKeyEvent>
#include <QColor>

#include "db/dbDatabase.h"

#include "LayoutView.h"
#include "LayoutScene.h"

using namespace db;

namespace gui
{

class GuiConfig;

class MainWindow : public QMainWindow
{
  Q_OBJECT

  public:

    MainWindow(QWidget *parent = nullptr);

    ~MainWindow();
  
    void init();
    void setDatabase(std::shared_ptr<dbDatabase> db) { db_ = db; }
    void setGuiConfig(std::shared_ptr<GuiConfig> config) { config_ = config; }

    std::shared_ptr<LayoutScene> getScene() { return layout_scene_; }
    std::shared_ptr<LayoutView>  getView()  { return layout_view_;  }

    void keyPressEvent(QKeyEvent* event) override;

  private:

    std::shared_ptr<dbDatabase>  db_;
    std::shared_ptr<GuiConfig>   config_;
    std::shared_ptr<LayoutScene> layout_scene_;
    std::shared_ptr<LayoutView>  layout_view_;

    void createMenu();
    void createDock();
    void createToolBar();
    void createItem();
};

}

#endif // MAINWINDOW_H
