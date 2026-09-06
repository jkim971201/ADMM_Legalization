#ifndef GUI_CONFIG_H
#define GUI_CONFIG_H

#include <QColor>

#include <memory>
#include <array>
#include <unordered_map>

#include "db/dbDatabase.h"
#include "db/dbTech.h"
#include "db/dbDesign.h"
#include "db/dbLayer.h"
#include "db/dbRegion.h"
#include "db/dbGroup.h"

namespace gui
{

using namespace db;

// Follow the color scheme of Innovus 23
// Color name list : w3.org/TR/SVG11/types.html#ColorKeywords
const std::array<QColor, 14> QCOLOR_ARRAY =
{
  QColor(  0,   0, 255), // 0  M1  : blue
  QColor(255,   0,   0), // 1  M2  : red
  QColor(  0, 255, 127), // 2  M3  : spring green
  QColor(139,   0,   0), // 3  M4  : dark red
  QColor(128,   0,   0), // 5  M5  : maroon 
  QColor(255, 165,   0), // 5  M6  : orange
  QColor(255,   0, 255), // 6  M7  : magenta
  QColor(  0, 255, 255), // 7  M8  : cyan
  QColor(240, 255, 255), // 8  ??  : Azure (TCAP LAYER D8 TCAP)
  QColor(160,  82,  45), // 9  M9  : sienna
  QColor(255, 255,   0), // 10 M10 : yellow 
  QColor(  0, 255,   0), // 11 M11 : green 
  QColor(255,  20, 255), // 12 M12 : magenta
  QColor(255,   0, 255)  // 13 M13 : fuchsia
}; // <- This is for Layer

class GuiConfig
{
  public:
  
    GuiConfig(std::shared_ptr<dbDatabase> db);

    const QColor getLayerColor(const dbLayer* layer);
    const QColor getGroupColor(const dbGroup* group);

    int dbu() const { return dbu_; }

    double getDieLxMicron() const { return dieLxMicron_; }
    double getDieLyMicron() const { return dieLyMicron_; }
    double getDieUxMicron() const { return dieUxMicron_; }
    double getDieUyMicron() const { return dieUyMicron_; }

    double getDieCxMicron() const { return (dieLxMicron_ + dieUxMicron_) / 2.0; }
    double getDieCyMicron() const { return (dieLyMicron_ + dieUyMicron_) / 2.0; }

    double getDieDxMicron() const { return dieUxMicron_ - dieLxMicron_; }
    double getDieDyMicron() const { return dieUyMicron_ - dieLyMicron_; }

  private: // should be private?

    void initLayerColor(std::shared_ptr<dbTech> tech);
    void initGroupColor(std::shared_ptr<dbDesign> design);

    int dbu_;
    std::unordered_map<const dbLayer*, QColor> layer2Color_;
    std::unordered_map<const dbGroup*, QColor> group2Color_;

    // NOTE : in micron
    double dieLxMicron_;
    double dieLyMicron_;
    double dieUxMicron_;
    double dieUyMicron_;
};

}

#endif
