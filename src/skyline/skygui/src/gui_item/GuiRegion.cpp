#include <QStyleOptionGraphicsItem>

#include "GuiRegion.h"

namespace gui
{


GuiRegion::GuiRegion(std::shared_ptr<GuiConfig> cfg,
                     const dbRegion* const region,
                     const dbBox* const box)
  : region_(region), box_(box)
{
  setConfig(cfg);
  this->setZValue(0);

  const double dbu = static_cast<double>(getConfig()->dbu());

  int region_lx_dbu = box->lx();
  int region_ly_dbu = box->ly();
  int region_ux_dbu = box->ux();
  int region_uy_dbu = box->uy();

  double region_lx_micron = region_lx_dbu / dbu;
  double region_ly_micron = region_ly_dbu / dbu;
  double region_ux_micron = region_ux_dbu / dbu;
  double region_uy_micron = region_uy_dbu / dbu;

  rect_ = QRectF(region_lx_micron, 
                 region_ly_micron, 
                 region_ux_micron - region_lx_micron, 
                 region_uy_micron - region_ly_micron);
}

QRectF
GuiRegion::boundingRect() const
{
  return rect_;
}

void
GuiRegion::paint(QPainter* painter, 
               const QStyleOptionGraphicsItem* option,
               QWidget* widget)
{
  const qreal lod 
    = option->levelOfDetailFromTransform(painter->worldTransform());
  
  auto db_group = region_->getGroup();
  auto color = getConfig()->getGroupColor(db_group);

  getPen().setColor(color);
  getBrush().setColor(color);
  getBrush().setStyle(Qt::BrushStyle::Dense2Pattern);
  //getBrush().setStyle(Qt::BrushStyle::DiagCrossPattern);

  GuiRect::paint(painter, option, widget);
}

}
