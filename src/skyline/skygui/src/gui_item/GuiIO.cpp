#include <cassert>
#include <cmath>
#include <iostream>

#include <QStyleOptionGraphicsItem>
#include "GuiIO.h"

namespace gui
{

GuiIOPort::GuiIOPort(std::shared_ptr<GuiConfig> cfg,
                     const dbBTermPort* port)
  : port_(port)
{
  setConfig(cfg);
  this->setZValue(port->layer()->index() + 1);

  const double dbu = static_cast<double>(getConfig()->dbu());

  int ioLx_dbu = port_->lx();
  int ioLy_dbu = port_->ly();
  int ioUx_dbu = port_->ux();
  int ioUy_dbu = port_->uy();
  int ioDx_dbu = std::max(port_->dx(), 1);
  int ioDy_dbu = std::max(port_->dy(), 1); 
  // Zero size wil make bug in the QPoint

  double ioLx_micron = ioLx_dbu / dbu;
  double ioLy_micron = ioLy_dbu / dbu;
  double ioUx_micron = ioUx_dbu / dbu;
  double ioUy_micron = ioUy_dbu / dbu;
  double ioDx_micron = ioDx_dbu / dbu;
  double ioDy_micron = ioDy_dbu / dbu;

  lx_ = ioLx_micron;
  ly_ = ioLy_micron;
  ux_ = ioLx_micron + ioDx_micron;
  uy_ = ioLy_micron + ioDy_micron;

  rect_ = QRectF(lx_, ly_, ux_ - lx_, uy_ - ly_);
}

void
GuiIOPort::paint(QPainter* painter, 
                 const QStyleOptionGraphicsItem* option,
                 QWidget* widget)
{
  const qreal lod 
    = option->levelOfDetailFromTransform(painter->worldTransform());
  
  const auto layer = port_->layer();
  auto color = getConfig()->getLayerColor(layer);

  getPen().setColor(color);
  getBrush().setColor(color);
  getBrush().setStyle(Qt::BrushStyle::DiagCrossPattern);

  GuiRect::paint(painter, option, widget);

  qreal io_lx = rect_.left();
  qreal io_ly = rect_.top();
  qreal io_ux = rect_.right();
  qreal io_uy = rect_.bottom();

  qreal io_cx = rect_.center().x();
  qreal io_cy = rect_.center().y();

  qreal len = std::min(std::abs(rect_.width()), 
                       std::abs(rect_.height())) * 1.0;

  if(lod < 10.0)
    len = len * 10.0;
  else 
    len = len;

  qreal p1X = 0.0;
  qreal p1Y = 0.0;
  qreal p2X = 0.0;
  qreal p2Y = 0.0;
  qreal p3X = 0.0;
  qreal p3Y = 0.0;

  // Couter clock-wise
  // 0.8660 ~= sqrt(3)/2
  constexpr qreal sqrt3_2 = std::sqrt(3) / 2.0;
  double die_cx_micron = getConfig()->getDieCxMicron();
  double die_cy_micron = getConfig()->getDieCyMicron();

  double die_lx_micron = getConfig()->getDieLxMicron();
  double die_ly_micron = getConfig()->getDieLyMicron();
  double die_ux_micron = getConfig()->getDieUxMicron();
  double die_uy_micron = getConfig()->getDieUyMicron();

  // Type 1 : ^
  // Type 2 : V
  // Type 3 : >
  // Type 4 : <

  int triangle_type = 0;
  if(io_cy < die_cy_micron)
  {
    double dist_to_bottom = std::abs(io_cy - die_ly_micron);
    if(io_cx < die_cx_micron)
    {
      // Lower Left
      double dist_to_left = std::abs(io_cx - die_lx_micron);
      triangle_type = dist_to_left > dist_to_bottom ? 1 : 3;
    }
    else
    {
      // Lower Right
      double dist_to_right = std::abs(io_cx - die_ux_micron);
      triangle_type = dist_to_right > dist_to_bottom ? 1 : 4;
    }
  }
  else
  {
    double dist_to_top = std::abs(io_cy - die_uy_micron);
    if(io_cx < die_cx_micron)
    {
      // Upper Left
      double dist_to_left = std::abs(io_cx - die_lx_micron);
      triangle_type = dist_to_left > dist_to_top ? 2 : 3;
    }
    else
    {
      // Upper Right
      double dist_to_right = std::abs(io_cx - die_ux_micron);
      triangle_type = dist_to_right > dist_to_top ? 2 : 4;
    }
  }

  switch(triangle_type)
  {
    case 1 :
    {
      p1X = io_cx;
      p1Y = io_ly;
      p2X = p1X - len * 0.5;
      p2Y = p1Y - len * sqrt3_2;
      p3X = p1X + len * 0.5;
      p3Y = p1Y - len * sqrt3_2;

      ux_ += len * 0.5;
      lx_ -= len * 0.5;
      ly_ -= len * sqrt3_2;
      break;
    }
    case 2 :
    {
      p1X = io_cx;
      p1Y = io_uy;
      p2X = p1X + len * 0.5;
      p2Y = p1Y + len * sqrt3_2;
      p3X = p1X - len * 0.5;
      p3Y = p1Y + len * sqrt3_2;

      ux_ += len * 0.5;
      lx_ -= len * 0.5;
      uy_ += len * sqrt3_2;
      break;
    }
    case 3 : 
    {
      p1X = io_lx;
      p1Y = io_cy;
      p2X = p1X - len * sqrt3_2;
      p2Y = p1Y + len * 0.5;
      p3X = p1X - len * sqrt3_2;
      p3Y = p1Y - len * 0.5;

      lx_ -= len * sqrt3_2;
      uy_ += len * 0.5;
      ly_ -= len * 0.5;
      break;
    }
    case 4 :
    {
      p1X = io_ux;
      p1Y = io_cy;
      p2X = p1X + len * sqrt3_2;
      p2Y = p1Y + len * 0.5;
      p3X = p1X + len * sqrt3_2;
      p3Y = p1Y - len * 0.5;

      ux_ += len * sqrt3_2;
      uy_ += len * 0.5;
      uy_ -= len * 0.5;
      break;
    }
    default:
      break;
  }

  //getPen().setColor(QColor(255, 215, 0)); // color name : gold
  //getBrush().setColor(QColor(255, 215, 0)); 

  getPen().setColor(color);
  getBrush().setColor(color);
  getBrush().setStyle(Qt::BrushStyle::SolidPattern);

  painter->setPen(getPen());
  painter->setBrush(getBrush());

  QPainterPath path;
  path.moveTo(p1X, p1Y);
  path.lineTo(p2X, p2Y);
  path.lineTo(p3X, p3Y);
  path.lineTo(p1X, p1Y);

  color.setAlphaF(0.8); // Set Transparency
  getBrush().setColor(color);
  painter->drawPath(path);

  drawIOName(painter, color, lod, triangle_type);
}

void
GuiIOPort::drawIOName(QPainter* painter, const QColor& color, qreal lod, int triangle_type)
{
  const qreal rect_width = rect_.width();
  const qreal rect_height = rect_.height();

  const qreal scale_adjust = 1.0 / lod;
  qreal text_font_size = rect_height;

  QString name(port_->getBTerm()->name().c_str());

  QPen pen = painter->pen();
  pen.setColor(color);
  painter->setPen(pen);

  painter->save();

  qreal orient_x = rect_.center().x();
  qreal orient_y = rect_.center().y();

  // delta below will be restored
  bool rotate90 = false;
  switch(triangle_type)
  {
    case 1:
    {
      text_font_size = rect_width * lod;
      orient_x = rect_.left();
      orient_y = rect_.bottom() - 2 * rect_height;
      rotate90 = true;
      break;
    }
    case 2:
    {
      text_font_size = rect_width * lod;
      orient_x = rect_.left();
      orient_y = rect_.bottom() + 2 * rect_height;
      rotate90 = true;
      break;
    }
    case 3:
    {
      text_font_size = rect_height * lod;
      orient_x = rect_.left() - 2 * rect_width;
      orient_y = rect_.top();
      break;
    }
    case 4:
    {
      text_font_size = rect_height * lod;
      orient_x = rect_.right() + 2 * rect_width;
      orient_y = rect_.top();
      break;
    }
    default:
      break;
  }

  QFont font = painter->font();
  font.setPointSizeF(text_font_size);
  painter->setFont(font);

  QRect text_bbox = painter->fontMetrics().boundingRect(name);
  if(triangle_type == 2)
    orient_y += text_bbox.width() * scale_adjust;
  if(triangle_type == 3)
    orient_x -= text_bbox.width() * scale_adjust;

  painter->translate(orient_x, orient_y);
  painter->scale(scale_adjust, -scale_adjust);

  if(rotate90 == true)
    painter->rotate(90);

  painter->drawText(0, 0, name);

  painter->restore();
}

QRectF
GuiIOPort::boundingRect() const
{
  return rect_;
}

/* GuiIO */
GuiIO::GuiIO(std::shared_ptr<GuiConfig> cfg,
             const dbBTerm* io)
  : io_(io)
{
  setConfig(cfg);

  lx_ = std::numeric_limits<double>::max();
  ly_ = std::numeric_limits<double>::max();
  ux_ = std::numeric_limits<double>::min();
  uy_ = std::numeric_limits<double>::min();

  const double dbu = static_cast<double>(config_->dbu());

  for(const auto port : io_->ports())
  {
    GuiIOPort* gui_port 
      = new GuiIOPort(config_, port);

    gui_ports_.push_back(gui_port);

    const auto bbox = gui_port->boundingRect();

    lx_ = std::min(lx_, bbox.left());
    ly_ = std::min(ly_, bbox.bottom());
    ux_ = std::max(ux_, bbox.right());
    uy_ = std::max(uy_, bbox.top());
  }
}

QRectF
GuiIO::boundingRect() const
{
  // QGrahicsView decide to re-paint a QGraphicsItem when
  // its boundingRect() is inside the scene.
  // If boundingRect() is not computed properly,
  // an item can disappear while repainting.
  return QRectF(lx_, ly_, ux_ - lx_, uy_ - ly_);
}

void
GuiIO::paint(QPainter* painter, 
             const QStyleOptionGraphicsItem* option,
             QWidget* widget)
{
  for(auto port : gui_ports_)
    port->paint(painter, option, widget);
}

}
