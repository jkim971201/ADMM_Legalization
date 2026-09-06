#include <cassert>
#include <iostream>
#include <regex>

#include "db/dbTech.h"
#include "db/dbDesign.h"
#include "db/dbDatabase.h"
#include "db/dbDie.h"
#include "db/dbInst.h"
#include "db/dbNet.h"
#include "db/dbITerm.h"
#include "db/dbBTerm.h"
#include "db/dbMTerm.h"
#include "db/dbRow.h"

#include "LayoutScene.h"
#include "GuiConfig.h"

#include "gui_item/GuiDie.h"
#include "gui_item/GuiRow.h"
#include "gui_item/GuiInst.h"
#include "gui_item/GuiIO.h"
#include "gui_item/GuiPin.h"
#include "gui_item/GuiNet.h"
#include "gui_item/GuiBlockage.h"
#include "gui_item/GuiTrackGrid.h"
#include "gui_item/GuiRegion.h"

namespace gui
{

LayoutScene::LayoutScene(
  QObject* parent,
  std::shared_ptr<dbDatabase> db,
  std::shared_ptr<GuiConfig> config) 
  : db_(db), config_(config)
{
  this->setBackgroundBrush(Qt::black);
}

void
LayoutScene::createGuiDie()
{
  const auto die = db_->getDesign()->getDie();
  GuiDie* die_gui = new GuiDie(die);

  double die_lx_micron = config_->getDieLxMicron();
  double die_ly_micron = config_->getDieLyMicron();

  double die_dx_micron = config_->getDieDxMicron();
  double die_dy_micron = config_->getDieDyMicron();

  die_gui->setRect( QRectF(die_lx_micron, die_ly_micron, 
                           die_dx_micron, die_dy_micron) );

  this->addItem(die_gui);
}

void
LayoutScene::createGuiRow()
{
  const double dbu = static_cast<double>(db_->getTech()->getDbu());

  for(auto row : db_->getDesign()->getRows())
  {
    GuiRow* row_gui = new GuiRow(row);

    double rowLx = static_cast<double>(row->lx()) / dbu;
    double rowLy = static_cast<double>(row->ly()) / dbu;
    double rowDx = static_cast<double>(row->dx()) / dbu;
    double rowDy = static_cast<double>(row->dy()) / dbu;
  
    row_gui->setRect( QRectF(rowLx, rowLy, rowDx, rowDy) );
    this->addItem(row_gui);
  }
}

void
LayoutScene::createGuiInst()
{
  const double dbu = static_cast<double>(db_->getTech()->getDbu());

  for(const auto inst : db_->getDesign()->getInsts())
  {
    GuiInst* inst_gui = new GuiInst(config_, inst);
  
    double cellLx = static_cast<double>(inst->lx()) / dbu;
    double cellLy = static_cast<double>(inst->ly()) / dbu;
    double cellDx = static_cast<double>(inst->dx()) / dbu;
    double cellDy = static_cast<double>(inst->dy()) / dbu;
  
    inst_gui->setRect( QRectF(cellLx, cellLy, cellDx, cellDy) );
    this->addItem(inst_gui);

    // NOTE : Both PIN and OBS will be converted to GuiPin.
    for(const auto iterm : inst->getITerms())
    {
      const auto mterm = iterm->getMTerm();
      const auto macro = inst->macro();
  
      for(const auto port : mterm->ports())
      {
        // MASTERSLICE is a nonrouting layer.
        // LEF DEF REF (May 2017) Page 45
        const dbLayer* layer = port->layer();
        if(layer->type() == RoutingType::MASTERSLICE)
          continue;
  
        GuiPin* cell_pin = new GuiPin(config_, false, inst, layer, port->getShape());
        inst_gui->addGuiPin(cell_pin);
      }
  
      for(const auto obs : macro->getObs())
      {
        const dbLayer* layer = obs->layer();
        if(obs->layer()->type() == RoutingType::MASTERSLICE)
          continue;
  
        GuiPin* cell_obs = new GuiPin(config_, true , inst, layer, obs->getShape());
        cell_obs->setZValue(layer->index() + 1);
        inst_gui->addGuiPin(cell_obs);
      }
    }
  }
}

void
LayoutScene::createGuiIO()
{
  for(const auto bterm : db_->getDesign()->getBTerms())
  {
    GuiIO* io_gui  = new GuiIO(config_, bterm);
    this->addItem(io_gui);
  }
}

void
LayoutScene::createGuiNet()
{
  for(const auto net : db_->getDesign()->getNets())
  {
    if(net->getWire()->getSegments().empty())
      continue;

    GuiNet* net_gui  = new GuiNet(config_, net);
    this->addItem(net_gui);
  }
}

void
LayoutScene::createGuiBlockage()
{
  for(const auto blk : db_->getDesign()->getBlockages())
  {
    if(blk->isPlacementBlockage() == true)
      continue;

    GuiBlockage* blk_gui  = new GuiBlockage(config_, blk);
    this->addItem(blk_gui);
  }
}

void
LayoutScene::createGuiTrackGrid()
{
  for(const auto grid : db_->getDesign()->getTrackGrids())
  {
    GuiTrackGrid* grid_gui  = new GuiTrackGrid(config_, grid);
    this->addItem(grid_gui);
  }
}

void
LayoutScene::createGuiRegion()
{
  for(const auto region : db_->getDesign()->getRegions())
  {
		for(auto& box : region->getRegionBoxes())
		{
      GuiRegion* region_gui = new GuiRegion(config_, region, &box);
      this->addItem(region_gui);
		}
  }
}

void
LayoutScene::expandScene()
{
  QRectF rect = sceneRect();
  double sceneW = rect.width();
  double sceneH = rect.height();

  // Make 20%  blank margin along boundary
  rect.adjust(-0.2 * sceneW, -0.2 * sceneH, 
              +0.2 * sceneW, +0.2 * sceneH);

  this->setSceneRect(rect);
}

}
