#include "Painter.h"

#include "SkyPlaceDB.h"
#include "object/GPObject.h"

#include "db/dbInst.h"

#include <random>
#include <QImage>
#include <QPainter>
#include <QApplication>

extern int    cmd_argc;
extern char** cmd_argv; 

namespace skyplace 
{

inline QRectF convertGPCellToQRectF(float scale, const GPCell* const cell)
{
  return QRectF(scale * cell->lx(), scale * cell->ly(), 
                scale * cell->dx(), scale * cell->dy());
}

inline QLineF getOrientLine(const QRectF& rect, dbInst* inst)
{
  qreal inst_width  = std::abs(rect.width());
  qreal inst_height = std::abs(rect.height());

  qreal delta_y = inst_height / 4.0;
  qreal delta_x = delta_y > inst_width ? 
                  inst_width / 4.0 : delta_y;

  QPointF p1;
  QPointF p2;

  switch(inst->orient())
  {
    case Orient::N  :
    case Orient::FW :
    {
      QPointF tl = rect.topLeft();
      p1 = tl + QPointF(0.0, delta_y);
      p2 = tl + QPointF(delta_x, 0.0);
      break;
    }
    case Orient::S  :
    case Orient::FE :
    {
      QPointF br = rect.bottomRight();
      p1 = br + QPointF(0.0, -delta_y);
      p2 = br + QPointF(-delta_x, 0.0);
      break;
    }
    case Orient::W  :
    case Orient::FN :
    {
      QPointF tr = rect.topRight();
      p1 = tr + QPointF(0.0, delta_y);
      p2 = tr + QPointF(-delta_x, 0.0);
      break;
    }
    case Orient::E  :
    case Orient::FS :
    {
      QPointF bl = rect.bottomLeft();
      p1 = bl + QPointF(0.0, -delta_y);
      p2 = bl + QPointF(delta_x, 0.0);
      break;
    }
    default:
      break;
  }

  QLineF line(p1, p2);
  return line;
}

Painter::Painter(std::shared_ptr<SkyPlaceDB> db)
  : db_(db) 
{
  qapp_ = std::make_unique<QApplication>(cmd_argc, cmd_argv);
}

Painter::~Painter() {}

void 
Painter::drawDieRect(float k_scale, QPainter* painter) const
{
  painter->save();
  float die_lx = db_->die()->lx();
  float die_ly = db_->die()->ly();
  float die_ux = db_->die()->ux();
  float die_uy = db_->die()->uy();

  float die_dx = db_->die()->dx();
  float die_dy = db_->die()->dy();

  // Draw Die
  QPen pen_for_die;
  pen_for_die.setColor(Qt::gray);
  pen_for_die.setStyle(Qt::PenStyle::DashDotLine);
  painter->setPen(pen_for_die);

  QRectF die_rect(k_scale * die_lx, k_scale * die_ly, k_scale * die_dx, k_scale * die_dy);
  painter->drawLine(die_rect.topLeft(), die_rect.topRight());
  painter->drawLine(die_rect.bottomLeft(), die_rect.bottomRight());
  painter->drawLine(die_rect.topLeft(), die_rect.bottomLeft());
  painter->drawLine(die_rect.topRight(), die_rect.bottomRight());
  painter->restore();
}

void
Painter::drawFixedRect(float k_scale, QPainter* painter) const
{
  // Draw Fixed Instances
  QBrush brush_for_macro(QColor(220, 220, 220), Qt::BrushStyle::Dense6Pattern);
  painter->setBrush(brush_for_macro);

  QPen pen_for_macro;
  pen_for_macro.setJoinStyle(Qt::PenJoinStyle::BevelJoin);
  pen_for_macro.setStyle(Qt::PenStyle::SolidLine);
  pen_for_macro.setColor(QColor(220, 220, 220));

  double line_width = 2.0 * pen_for_macro.widthF();
  pen_for_macro.setWidthF(line_width);

  painter->setPen(pen_for_macro);

  for(const auto fixed_cell : db_->fixedCells())
  {
    auto macro_rect = convertGPCellToQRectF(k_scale, fixed_cell);
    painter->drawLine(macro_rect.topLeft(), macro_rect.topRight());
    painter->drawLine(macro_rect.bottomLeft(), macro_rect.bottomRight());
    painter->drawLine(macro_rect.topLeft(), macro_rect.bottomLeft());
    painter->drawLine(macro_rect.topRight(), macro_rect.bottomRight());
    painter->drawRect(macro_rect);

    auto orient_line = getOrientLine(macro_rect, fixed_cell->dbInstPtr());
    painter->drawLine(orient_line);
  }
}

void
Painter::drawInstRect(float k_scale, QPainter* painter) const
{
  // Draw Standard Cells
  QBrush brush_for_std(Qt::gray, Qt::BrushStyle::Dense6Pattern);
  painter->setBrush(brush_for_std);

  QPen pen_for_std;
  pen_for_std.setJoinStyle(Qt::PenJoinStyle::BevelJoin);
  pen_for_std.setStyle(Qt::PenStyle::SolidLine);

  const int num_cluster = db_->numCluster();
  for(const auto movable_inst : db_->movableCells())
  {
    if(movable_inst->isFiller() == true)
      continue;

    int cluster_id = movable_inst->clusterID();
    auto find_color = clusterId2Color_.find(cluster_id);
    QColor inst_color = num_cluster > 0 ? find_color->second : Qt::gray;
    inst_color.setAlphaF(0.75);
    pen_for_std.setColor(inst_color);
    painter->setPen(pen_for_std);

    auto std_rect = convertGPCellToQRectF(k_scale, movable_inst);
    painter->drawLine(std_rect.topLeft(), std_rect.topRight());
    painter->drawLine(std_rect.bottomLeft(), std_rect.bottomRight());
    painter->drawLine(std_rect.topLeft(), std_rect.bottomLeft());
    painter->drawLine(std_rect.topRight(), std_rect.bottomRight());
    painter->drawRect(std_rect);
  }
}

void
Painter::saveImage(int iter, float hpwl, float overflow)
{
  float die_dx = db_->die()->ux();
  float die_dy = db_->die()->uy();

  const float k_scale_x = 500 / die_dx;
  const float k_scale_y = 500 / die_dy;
  const float k_scale = std::min(k_scale_x, k_scale_y);
  constexpr float k_margin = 0.2;

  int image_size_x = static_cast<int>(k_scale * (1 + k_margin) * die_dx);
  int image_size_y = static_cast<int>(k_scale * (1 + k_margin) * die_dy);

  QImage image(image_size_x, image_size_y, QImage::Format_RGB32);
  image.fill(Qt::black);
  QPainter painter(&image); 
  // When image is null (isNull() returns true),
  // QPainter will be not active.

  float offset_x = k_scale * k_margin * die_dx / 2.0;
  float offset_y = k_scale * k_margin * die_dy / 2.0;
  painter.translate(offset_x, offset_y);

  drawDieRect(k_scale, &painter);
  
  drawFixedRect(k_scale, &painter);

  drawInstRect(k_scale, &painter);

  auto get_padded_string = [] (size_t size, const std::string& str)
  {
    if(str.size() < size)
      return std::string(size - str.size(), ' ') + str;
    else
      return str;
  };

  std::string info;
  info += "Iter: ";
  info += get_padded_string(5, std::to_string(iter));
  info += "   HPWL: ";
  info += get_padded_string(13, std::to_string(static_cast<int64_t>(hpwl)) );
  info += "   Overflow: ";
  info += std::to_string(overflow);

  QFont font = painter.font();
  font.setBold(true);
  const qreal text_font_size = k_scale * 0.025 * die_dy;
  font.setPointSizeF(text_font_size);
  painter.setFont(font);
  painter.scale(1.0, -1.0);
  auto pen_for_text = painter.pen();
  pen_for_text.setColor(Qt::white);
  painter.setPen(pen_for_text);

  QString metric_qstring(info.c_str());
  painter.drawText(0, offset_y / 2.0, metric_qstring);

  std::ostringstream iter_string;
  iter_string << std::setw(4) << std::setfill('0') << iter;

  std::string jpeg_name = iter_string.str() + ".jpeg";
  std::filesystem::path file_name = path_to_save_image_ / jpeg_name;
  QImage image_to_save = image.mirrored(false, true);
  image_to_save.save(file_name.string().c_str(), "jpeg", 50);
}

void 
Painter::preparePlot(const std::string& plot_path)
{
  path_to_save_image_ = std::filesystem::path(plot_path);

  if(std::filesystem::exists(path_to_save_image_) == false)
  {
    //if(std::filesystem::is_directory(path_to_save_image_) == false)
    std::string command = "mkdir -p " + path_to_save_image_.string();
    std::system(command.c_str());
  }

  const int num_cluster = db_->numCluster();
  const int opacity = 255;
  for(int i = 0; i < num_cluster; i++)
  {
    std::mt19937 rand_gen_1(3 * i + 0);
    std::mt19937 rand_gen_2(3 * i + 1);
    std::mt19937 rand_gen_3(3 * i + 2);
    std::uniform_int_distribution<int> dist(0, 255);
    int r_val = dist(rand_gen_1);
    int g_val = dist(rand_gen_2);
    int b_val = dist(rand_gen_3);
    QColor color(r_val, g_val, b_val, opacity);
    clusterId2Color_[i] = color;
  }
}

} 
