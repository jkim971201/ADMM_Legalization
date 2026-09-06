#include <random>

#include "GuiConfig.h"

#include "db/dbDie.h"

namespace gui
{

GuiConfig::GuiConfig(std::shared_ptr<dbDatabase> db)
{
  const auto tech   = db->getTech();
  const auto design = db->getDesign();

  // We keep die coordinates 
  // to share the boundary information of the layout
  // with all the gui items.
  dbu_ = tech->getDbu();
  const double dbu_double = static_cast<double>(dbu_);
  const auto die = design->getDie();
  dieLxMicron_ = die->lx() / dbu_double;
  dieLyMicron_ = die->ly() / dbu_double;
  dieUxMicron_ = die->ux() / dbu_double;
  dieUyMicron_ = die->uy() / dbu_double;

  initLayerColor(tech);
  initGroupColor(design);
}

void
GuiConfig::initLayerColor(std::shared_ptr<dbTech> tech)
{
  QColor color;
  int routingLayerIdx = 0; // Index for ROUTING LAYER
  int cutLayerIdx     = 1; // Index for     CUT LAYER
  // Cut Layer matches to its upper layer.
  for(auto layer : tech->getLayers())
  {
    auto type = layer->type();
    if(type == RoutingType::ROUTING)
    {
      if(routingLayerIdx < QCOLOR_ARRAY.size())
      {
        color = QCOLOR_ARRAY.at(routingLayerIdx);
        routingLayerIdx++;
      }
      else 
        color = QColor(240, 255, 255); // color name : azure 
    }
    else if(type == RoutingType::CUT)
    {
      if(cutLayerIdx < QCOLOR_ARRAY.size())
      {
        color = QCOLOR_ARRAY.at(cutLayerIdx);
        cutLayerIdx++;
      }
      else 
        color = QColor(240, 255, 255); // color name : azure 
    }
    else
      color = QColor(240, 255, 255); // color name : azure 

    layer2Color_[layer] = color;
  }
}

void
GuiConfig::initGroupColor(std::shared_ptr<dbDesign> design)
{
  const auto& db_groups = design->getGroups();
  const int num_groups = static_cast<int>(db_groups.size());
  
	// Fixed seed to see fixed color for every run.

  const int group_opacity = 150; // default opacity = 255
  for(int i = 0; i < num_groups; i++)
  {
    std::mt19937 rand_gen_1(3 * i + 0);
    std::mt19937 rand_gen_2(3 * i + 1);
    std::mt19937 rand_gen_3(3 * i + 2);
    std::uniform_int_distribution<int> dist(0, 255);
    int r_val = dist(rand_gen_1);
    int g_val = dist(rand_gen_2);
    int b_val = dist(rand_gen_3);
    QColor color(r_val, g_val, b_val, group_opacity);
    group2Color_[db_groups.at(i)] = color;
  }
}

const QColor
GuiConfig::getLayerColor(const dbLayer* layer)
{
  auto itr = layer2Color_.find(layer);
  if(itr == layer2Color_.end())
    assert(0); // TODO : Exception Handling
  return itr->second;
}

const QColor
GuiConfig::getGroupColor(const dbGroup* group)
{
  auto itr = group2Color_.find(group);
  if(itr == group2Color_.end())
    assert(0); // TODO : Exception Handling
  return itr->second;
}

}
