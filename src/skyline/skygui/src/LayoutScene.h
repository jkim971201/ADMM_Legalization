#ifndef LAYOUT_SCENE_H
#define LAYOUT_SCENE_H

#include <memory>
#include <QColor>
#include <QPainter>
#include <QGraphicsScene>

#include "db/dbDatabase.h"
#include "skygui/SkyGui.h"

using namespace db;

namespace gui
{

class GuiConfig;

class LayoutScene : public QGraphicsScene
{
  Q_OBJECT

  public:
    
    LayoutScene(
      QObject* parent,
      std::shared_ptr<dbDatabase> db,
      std::shared_ptr<GuiConfig> config);

    void setConfig(std::shared_ptr<GuiConfig> config);
    void createGuiDie();
    void createGuiRow();
    void createGuiInst();
    void createGuiIO();
    void createGuiNet();
    void createGuiBlockage();
    void createGuiTrackGrid();
    void createGuiRegion();

    void expandScene();

  private:

    std::shared_ptr<dbDatabase> db_;
    std::shared_ptr<GuiConfig> config_;
};

}

#endif
