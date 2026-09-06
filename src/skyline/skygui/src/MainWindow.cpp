#include <QMenu>
#include <QAction>
#include <QMenuBar>
#include <QToolBar>
#include <QListWidget>
#include <QTreeWidget>
#include <QHeaderView>
#include <QDockWidget>
#include <QStandardItem>
#include <QStatusBar>
#include <QDebug>
#include <QDesktopWidget>
#include <QGuiApplication>
#include <QScreen>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QPainter>
#include <QColor>

#include <cassert>
#include <iostream>

#include "MainWindow.h"
#include "GuiConfig.h"

#include "db/dbTech.h"
#include "db/dbTypes.h"
#include "db/dbDesign.h"

static void loadResources()
{
  Q_INIT_RESOURCE(resource);
}

namespace gui
{

MainWindow::MainWindow(QWidget *parent) {}

MainWindow::~MainWindow() {}

void
MainWindow::init()
{
  assert(db_ != nullptr);

  loadResources();

  setWindowTitle("SkyLine");

  // Scene
  layout_scene_ = std::make_shared<LayoutScene>(this, db_, config_);

  // View
  layout_view_ = std::make_shared<LayoutView>();
  layout_view_->setScene(layout_scene_.get());
  setCentralWidget(layout_view_.get());

  // Menu Bar
  createMenu();

  // Dock
  createDock();

  // Tool Bar
  createToolBar();

  // Status Bar
  statusBar()->showMessage(tr("Ready"));

  // This line is only for Qt5
  // I don't know why, but this returns merged size when using multiple monitors.
  QSize screenSize = QGuiApplication::primaryScreen()->size();

  // if-else to handle the problem above.
  QSize size = (screenSize.width() > 5000) ? screenSize * 0.4 
                                           : screenSize * 0.8;
  resize(size);

  // Draw Objects
  createItem();
}

void
MainWindow::createMenu()
{
  QMenu* menu;

  QFont font = menuBar()->font();
  font.setPointSize(12);
  menuBar()->setFont( font );

  menu = menuBar()->addMenu(tr("&File"));
  menu = menuBar()->addMenu(tr("&View"));
  menu = menuBar()->addMenu(tr("&Help"));
}

void
MainWindow::createDock()
{
  struct DockRow
  {
    QStandardItem* name;      // First  Item
    QStandardItem* color_box; // Second Item
    QStandardItem* visible;   // Third  Item
  };

  // DockWidget always...
  // 1. appear ahead of main window.
  // 2. separable from main window.
  QDockWidget* dock_widget = new QDockWidget(tr("GUI Control"), this);
  // dock_widget->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

  QFont dock_font = dock_widget->font();
  dock_font.setPointSizeF(13.0);
  dock_widget->setFont(dock_font);

  QTreeView* tree_view = new QTreeView(dock_widget);

  const int num_dock_column = 3;
  QStandardItemModel* tree_model = new QStandardItemModel(0, num_dock_column, dock_widget);
  tree_view->setModel(tree_model);
  tree_view->setContextMenuPolicy(Qt::CustomContextMenu);

  // Set Headers
  tree_model->setHorizontalHeaderLabels({"Name", "C", "V"});

  QHeaderView* header = tree_view->header();
  for(int column_index = 0; column_index < num_dock_column; column_index++)
  {
    auto resize_option = column_index == 0 ? QHeaderView::Stretch
                                           : QHeaderView::ResizeToContents;
    // auto resize_option = QHeaderView::ResizeToContents;
    header->setSectionResizeMode(column_index, resize_option);
  }
  
  // Make Top Rows
  QStandardItem* root = tree_model->invisibleRootItem();
  
  // 1. Layers
  DockRow layer_row_parent =
    { new QStandardItem(QString("Layers")),
      new QStandardItem(),
      new QStandardItem() };
  layer_row_parent.visible->setCheckable(true);

  root->appendRow({layer_row_parent.name, 
                   layer_row_parent.color_box,
                   layer_row_parent.visible});

  for(const auto layer : db_->getTech()->getLayers())
  {
    if(layer->type() != RoutingType::ROUTING)
      continue;

    QList<QStandardItem*> layer_row_child;
    QStandardItem* name = new QStandardItem(QString(layer->name().c_str()));
    layer_row_child << name;

    QColor layer_color = config_->getLayerColor(layer);
    QPixmap pix_map(20, 20);
    pix_map.fill(layer_color);
    QIcon color_icon = QIcon(pix_map);

    QStandardItem* color_box = new QStandardItem(color_icon, QString());
    color_box->setCheckable(false);
    layer_row_child << color_box;

    QStandardItem* visible = new QStandardItem();
    visible->setCheckable(true);
    layer_row_child << visible;

    layer_row_parent.name->appendRow(layer_row_child);
  }

  // 2. Instances
  DockRow insts_row_parent =
    { new QStandardItem(QString("Instances")),
      new QStandardItem(),
      new QStandardItem() };
  insts_row_parent.visible->setCheckable(true);

  root->appendRow({insts_row_parent.name,
                   insts_row_parent.color_box,
                   insts_row_parent.visible});

  // 3. Tracks
  DockRow tracks_row_parent =
    { new QStandardItem(QString("Tracks")),
      new QStandardItem(),
      new QStandardItem() };
  tracks_row_parent.visible->setCheckable(true);

  root->appendRow({tracks_row_parent.name,
                   tracks_row_parent.color_box,
                   tracks_row_parent.visible});

  // 4. Rows 
  DockRow rows_row_parent =
    { new QStandardItem(QString("Rows")),
      new QStandardItem(),
      new QStandardItem() };
  rows_row_parent.visible->setCheckable(true);

  root->appendRow({rows_row_parent.name,
                   rows_row_parent.color_box,
                   rows_row_parent.visible});

  // 5. Groups 
  DockRow group_row_parent =
    { new QStandardItem(QString("Groups")),
      new QStandardItem(),
      new QStandardItem() };
  group_row_parent.visible->setCheckable(true);

  root->appendRow({group_row_parent.name, 
                   group_row_parent.color_box,
                   group_row_parent.visible});

  for(const auto group : db_->getDesign()->getGroups())
  {
    QList<QStandardItem*> group_row_child;
    QStandardItem* name = new QStandardItem(QString(group->getName().data()));
    group_row_child << name;

    QColor group_color = config_->getGroupColor(group);
    QPixmap pix_map(20, 20);
    pix_map.fill(group_color);
    QIcon color_icon = QIcon(pix_map);

    QStandardItem* color_box = new QStandardItem(color_icon, QString());
    color_box->setCheckable(false);
    group_row_child << color_box;

    QStandardItem* visible = new QStandardItem();
    visible->setCheckable(true);
    group_row_child << visible;

    group_row_parent.name->appendRow(group_row_child);
  }

  tree_view->expandAll();
  dock_widget->setWidget(tree_view);
  addDockWidget(Qt::RightDockWidgetArea, dock_widget);
}

void
MainWindow::createToolBar()
{
  QToolBar* toolBar;

  QAction* zoomIn  = new QAction(QIcon(":/zoom_in.png") , tr("Zoom In") , this);
  QAction* zoomOut = new QAction(QIcon(":/zoom_out.png"), tr("Zoom Out"), this);
  QAction* zoomFit = new QAction(QIcon(":/zoom_fit.png"), tr("Zoom Fit"), this);

  toolBar = addToolBar(tr("Tool Bar"));
  toolBar->addAction( zoomIn  );
  toolBar->setStatusTip(tr("Zoom In Layout View"));
  connect(zoomIn, SIGNAL(triggered()), layout_view_.get(), SLOT(zoomIn_slot()));

  toolBar->addAction( zoomOut );
  toolBar->setStatusTip(tr("Zoom Out Layout View"));
  connect(zoomOut, SIGNAL(triggered()), layout_view_.get(), SLOT(zoomOut_slot()));

  toolBar->addAction( zoomFit );
  toolBar->setStatusTip(tr("Zoom Fit Layout View"));
  connect(zoomFit, SIGNAL(triggered()), layout_view_.get(), SLOT(zoomFit_slot()));
}

void
MainWindow::createItem()
{
  layout_scene_->createGuiDie();
  layout_scene_->createGuiRow();
  layout_scene_->createGuiInst();
  layout_scene_->createGuiIO();
  layout_scene_->createGuiNet();
  layout_scene_->createGuiBlockage();
  // layout_scene_->createGuiTrackGrid();
  layout_scene_->createGuiRegion();
  layout_scene_->expandScene();
}

void
MainWindow::keyPressEvent(QKeyEvent* event)
{
  if(event->key() == Qt::Key_Q)
    close();
}

}
