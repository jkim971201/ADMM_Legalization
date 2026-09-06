#ifndef GUI_REGION_H
#define GUI_REGION_H

#include "db/dbRegion.h"

#include "GuiItem.h"
#include "GuiRect.h"

using namespace db;

namespace gui
{

class GuiRegion : public GuiRect
{
  public:

    GuiRegion(std::shared_ptr<GuiConfig> cfg,
              const dbRegion* const region,
              const dbBox* const box);

    QRectF boundingRect() const override;

    void paint(QPainter* painter, 
               const QStyleOptionGraphicsItem* option, 
               QWidget* widget) override;

  private:

    const dbBox* box_;
    const dbRegion* region_;
};

}

#endif
