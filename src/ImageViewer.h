#pragma once

#include <QtWidgets/QMainWindow>
#include <QtWidgets>
#include "ui_ImageViewer.h"
#include "ViewerWidget.h"
#include "IPmodul.h"

class ImageViewer : public QMainWindow
{
	Q_OBJECT

public:
	ImageViewer(QWidget* parent = Q_NULLPTR);

private:
	Ui::ImageViewerClass* ui;
	ViewerWidget* vW;

	QSettings settings;
	QMessageBox msgBox;

	//ImageViewer Events
	void closeEvent(QCloseEvent* event);

	//Image functions
	bool openImage(QString filename);
	bool saveImage(QString filename);
	bool invertColors();

private slots:
	void on_actionOpen_triggered();
	void on_actionSave_as_triggered();
	void on_actionExit_triggered();
	void on_actionInvert_triggered();

	void on_actionPrint_histogram_triggered();
	void on_actionFSHS_triggered();
	void on_actionEKV_HIST_triggered();

	void on_actionConvolution_triggered();
	void on_actionExplicit_Heat_Eq_triggered();
	void on_actionImplicit_Heat_Eq_triggered();

	void on_actionPerona_Malik_model_triggered();

	void on_actionMCF_triggered();

	void on_actionDistance_function_triggered();

	void on_checkBox_UseGMCF_clicked(bool isChecked);

	void on_pushButton_mirrorTest_clicked();
};
