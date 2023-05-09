#include "ImageViewer.h"

ImageViewer::ImageViewer(QWidget* parent)
	: QMainWindow(parent), ui(new Ui::ImageViewerClass)
{
	ui->setupUi(this);
	vW = new ViewerWidget(QSize(500, 500));
	ui->scrollArea->setWidget(vW);

	ui->scrollArea->setBackgroundRole(QPalette::Dark);
	ui->scrollArea->setWidgetResizable(true);
	ui->scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
	ui->scrollArea->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);

	vW->setObjectName("ViewerWidget");
}

//ImageViewer Events
void ImageViewer::closeEvent(QCloseEvent* event)
{
	if (QMessageBox::Yes == QMessageBox::question(this, "Close Confirmation", "Are you sure you want to exit?", QMessageBox::Yes | QMessageBox::No))
	{
		event->accept();
	}
	else {
		event->ignore();
	}
}

//Image functions
bool ImageViewer::openImage(QString filename)
{
	QImage loadedImg(filename);
	if (!loadedImg.isNull()) {
		return vW->setImage(loadedImg);
	}
	return false;
}
bool ImageViewer::saveImage(QString filename)
{
	QFileInfo fi(filename);
	QString extension = fi.completeSuffix();

	QImage* img = vW->getImage();
	return img->save(filename, extension.toStdString().c_str());
}

bool ImageViewer::invertColors()
{
	if (vW->isEmpty()) {
		return false;
	}

	uchar* data = vW->getData();

	int row = vW->getImage()->bytesPerLine();
	int depth = vW->getImage()->depth();

	for (int i = 0; i < vW->getImgHeight(); i++)
	{
		for (int j = 0; j < vW->getImgWidth(); j++)
		{
			//Grayscale
			if (depth == 8) {
				vW->setPixel(j, i, static_cast<uchar>(255 - data[i * row + j]));
			}
			//RGBA
			else {
				uchar r = static_cast<uchar>(255 - data[i * row + j * 4]);
				uchar g = static_cast<uchar>(255 - data[i * row + j * 4 + 1]);
				uchar b = static_cast<uchar>(255 - data[i * row + j * 4 + 2]);
				vW->setPixel(j, i, r, g, b);
			}
		}
	}
	vW->update();
	return true;
}

//Slots
void ImageViewer::on_actionOpen_triggered()
{
	QString folder = settings.value("folder_img_load_path", "").toString();

	QString fileFilter = "Image data (*.bmp *.gif *.jpg *.jpeg *.png *.pbm *.pgm *.ppm .*xbm .* xpm);;All files (*)";
	QString fileName = QFileDialog::getOpenFileName(this, "Load image", folder, fileFilter);
	if (fileName.isEmpty()) { return; }

	QFileInfo fi(fileName);
	settings.setValue("folder_img_load_path", fi.absoluteDir().absolutePath());

	if (!openImage(fileName)) {
		msgBox.setText("Unable to open image.");
		msgBox.setIcon(QMessageBox::Warning);
		msgBox.exec();
	}
}
void ImageViewer::on_actionSave_as_triggered()
{
	QString folder = settings.value("folder_img_save_path", "").toString();

	QString fileFilter = "Image data (*.bmp *.gif *.jpg *.jpeg *.png *.pbm *.pgm *.ppm .*xbm .* xpm);;All files (*)";
	QString fileName = QFileDialog::getSaveFileName(this, "Save image", folder, fileFilter);
	if (!fileName.isEmpty()) {
		QFileInfo fi(fileName);
		settings.setValue("folder_img_save_path", fi.absoluteDir().absolutePath());

		if (!saveImage(fileName)) {
			msgBox.setText("Unable to save image.");
			msgBox.setIcon(QMessageBox::Warning);
		}
		else {
			msgBox.setText(QString("File %1 saved.").arg(fileName));
			msgBox.setIcon(QMessageBox::Information);
		}
		msgBox.exec();
	}
}
void ImageViewer::on_actionExit_triggered()
{
	this->close();
}

void ImageViewer::on_actionInvert_triggered()
{
	invertColors();
}

void ImageViewer::on_actionPrint_histogram_triggered()
{
	IPmodul ipmodul;
	uint* hist = nullptr;
	
	ipmodul.computeHistogramData(vW->getData(), vW->getImage()->bytesPerLine(), vW->getImage()->width(), vW->getImage()->height());
	hist = ipmodul.getHistogram();
	
	for (int i = 0; i < 256; i++)
	{
		printf("hist[%d]: %d\n", i, hist[i]);
	}
	printf("\n");
}

void ImageViewer::on_actionFSHS_triggered()
{
	if (vW->isEmpty()) {
		return;
	}
	
	IPmodul ipmodul;

	ipmodul.FSHS(vW->getData(), vW->getImage()->bytesPerLine(), vW->getImgWidth(), vW->getImgHeight());

	vW->update();
}

void ImageViewer::on_actionEKV_HIST_triggered()
{
	if (vW->isEmpty()) {
		return;
	}

	IPmodul ipmodul;
	
	ipmodul.EKV_HIST(vW->getData(), vW->getImage()->bytesPerLine(), vW->getImgWidth(), vW->getImgHeight());
	
	vW->update();
}

void ImageViewer::on_actionConvolution_triggered()
{
	if (vW->isEmpty()) {
		return;
	}

	IPmodul ipmodul;
	ConvolutionKernel kernel;

	double d[25] = { 0.000664574, 0.00600398, 0.01241,  0.00600398, 0.000664574,
					 0.00600398,  0.0542418,  0.112116, 0.0542418,  0.00600398,
					 0.01241,     0.112116,   0.234237, 0.112116,   0.01241,
					 0.00600398,  0.0542418,  0.112116, 0.0542418,  0.00600398,
					 0.000664574, 0.00600398, 0.01241,  0.00600398, 0.000664574
					};
	kernel.setKernel(5, d);
	//kernel.printKernel();

	// compute convolution
	uchar* newImg = ipmodul.convolution(vW->getData(), vW->getImage()->bytesPerLine(), vW->getImgWidth(), vW->getImgHeight(), &kernel);

	// copy new image values
	for (int i = 0; i < vW->getImgHeight(); i++)
	{
		for (int j = 0; j < vW->getImgWidth(); j++)
		{
			vW->getData()[i * vW->getBytesPerLine() + j] = newImg[i * vW->getImgWidth() + j];
		}
	}

	// free memory
	delete[] newImg;

	vW->update();
	//QString fileName = QInputDialog::getText(this, tr("Image export"), tr("File name:"));
	//fileName.prepend("../temp/");
	//IPmodul::exportToPGM(fileName.toStdString(), vW->getImgWidth(), vW->getImgHeight(), 255, newImg);

	//QMessageBox::information(this, "File export", "Image export done.", QMessageBox::Ok);
}

void ImageViewer::on_actionExplicit_Heat_Eq_triggered()
{
	if (vW->isEmpty()) {
		return;
	}

	IPmodul ip;
	int timeSteps = ui->spinBox_HeatEqTimeSteps->value();
	double tau = ui->doubleSpinBox_HeatEqTau->value();
	double h = 1.0;
	
	uchar* newImg = ip.filtrationExplicitHeatEq(vW->getData(), vW->getBytesPerLine(), vW->getImgWidth(), vW->getImgHeight(), tau, h, timeSteps);

	// copy new image values
	for (int i = 0; i < vW->getImgHeight(); i++)
	{
		for (int j = 0; j < vW->getImgWidth(); j++)
		{
			vW->getData()[i * vW->getBytesPerLine() + j] = newImg[i * vW->getImgWidth() + j];
		}
	}

	// free memory
	delete[] newImg;

	vW->update();
}

void ImageViewer::on_actionImplicit_Heat_Eq_triggered()
{
	if (vW->isEmpty()) {
		return;
	}

	IPmodul ip;
	//ip.printMsg = false;
	int timeSteps = ui->spinBox_HeatEqTimeSteps->value();
	double tau = ui->doubleSpinBox_HeatEqTau->value();

	uchar* newImg = ip.filtrationImplicitHeatEq(vW->getData(), vW->getBytesPerLine(), vW->getImgWidth(), vW->getImgHeight(), tau, timeSteps);

	//IPmodul::exportToPGM("../temp/lena_impl", vW->getImgWidth(), vW->getImgHeight(), 255, newImg);

	// copy new image values
	for (int i = 0; i < vW->getImgHeight(); i++)
	{
		for (int j = 0; j < vW->getImgWidth(); j++)
		{
			vW->getData()[i * vW->getBytesPerLine() + j] = newImg[i * vW->getImgWidth() + j];
		}
	}

	// free memory
	delete[] newImg;

	vW->update();
}

void ImageViewer::on_actionPerona_Malik_model_triggered()
{
	if (vW->isEmpty()) {
		return;
	}

	IPmodul ip;
	//ip.printMsg = false;
	int timeSteps = ui->spinBox_nonLinearDiffTimeSteps->value();
	double tau = ui->doubleSpinBox_NonLinearDiffTau->value();
	double sigma = ui->doubleSpinBox_NonLinearDiffSigma->value();
	double K = ui->doubleSpinBox_NonLinearDiffParameterK->value();
	
	uchar* newImg = ip.filtrationSemiImplicitPeronaMalik(vW->getData(), vW->getBytesPerLine(), vW->getImgWidth(), vW->getImgHeight(), sigma, tau, K, timeSteps);

	//IPmodul::exportToPGM("../temp/lena_PM", vW->getImgWidth(), vW->getImgHeight(), 255, newImg);

	// copy new image values
	for (int i = 0; i < vW->getImgHeight(); i++)
	{
		for (int j = 0; j < vW->getImgWidth(); j++)
		{
			vW->getData()[i * vW->getBytesPerLine() + j] = newImg[i * vW->getImgWidth() + j];
		}
	}

	// free memory
	delete[] newImg;

	vW->update();
}

void ImageViewer::on_actionMCF_triggered()
{
	if (vW->isEmpty()) {
		return;
	}

	IPmodul ip;
	//ip.printMsg = false;
	int timeSteps = ui->spinBox_GMCFTimeSteps->value();
	double tau = ui->doubleSpinBox_GMCFTau->value();
	double sigma = ui->doubleSpinBox_GMCFSigma->value();
	double K = 0.0;

	if (ui->checkBox_UseGMCF->isChecked())
		K = ui->doubleSpinBox_GMCFParameterK->value();
	
	uchar* newImg = nullptr;

	if (ui->checkBox_UseBiCGStab->isChecked())
	{
		auto start = std::chrono::high_resolution_clock::now();

		newImg = ip.filtrationSemiImplicitGMCF_BiCGStab(vW->getData(), vW->getBytesPerLine(), vW->getImgWidth(), vW->getImgHeight(), sigma, tau, K, timeSteps);

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

		printf("BiCGStab: %.4lf seconds\n\n", (double)duration.count() / 1000000.0);
	}
	else
	{
		auto start = std::chrono::high_resolution_clock::now();

		newImg = ip.filtrationSemiImplicitGMCF_SOR(vW->getData(), vW->getBytesPerLine(), vW->getImgWidth(), vW->getImgHeight(), sigma, tau, K, timeSteps);

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

		printf("SOR: %.4lf seconds\n\n", (double)duration.count() / 1000000.0);

	}

	//IPmodul::exportToPGM("../temp/lena_PM", vW->getImgWidth(), vW->getImgHeight(), 255, newImg);

	if (newImg == nullptr)
	{
		return;
	}

	// copy new image values
	for (int i = 0; i < vW->getImgHeight(); i++)
	{
		for (int j = 0; j < vW->getImgWidth(); j++)
		{
			vW->getData()[i * vW->getBytesPerLine() + j] = newImg[i * vW->getImgWidth() + j];
		}
	}

	// free memory
	delete[] newImg;

	vW->update();
}

void ImageViewer::on_actionDistance_function_triggered()
{
	printf("signed distance function\n");

	if (vW->isEmpty()) {
		return;
	}

	IPmodul ip;

	double* data = ip.signedDistanceFunction(vW->getData(), vW->getBytesPerLine(), vW->getImgWidth(), vW->getImgHeight());

	std::string fileName = "../temp/distFuncTest_edges";
	
	IPmodul::exportSgnFunction(fileName, vW->getImgWidth(), vW->getImgHeight(), data);
}

void ImageViewer::on_checkBox_UseGMCF_clicked(bool isChecked)
{
	if (!isChecked)
	{
		ui->doubleSpinBox_GMCFParameterK->setEnabled(false);
	}
	else
	{
		ui->doubleSpinBox_GMCFParameterK->setEnabled(true);
	}
}

void ImageViewer::on_pushButton_mirrorTest_clicked()
{
	if (vW->getImage() == nullptr)
		return;
	
	IPmodul ipmodul;
	uint N = QInputDialog::getInt(this, tr("Mirror pixels"), tr("Number of pixels:"), 2, 1);
	
	//ipmodul.pixelsMirror(d, 4, 4, N);
	bool result = ipmodul.pixelsMirror(vW->getData(), vW->getImage()->bytesPerLine(), vW->getImage()->width(), vW->getImage()->height(), N);

	printf("mirror pixels result: %d\n", result);
	

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < ipmodul.getImgWidth(); j++)
		{
			printf("%.1lf ", ipmodul.getImgData()[i * ipmodul.getImgWidth() + j]);
		}
		printf("\n\n\n");
	}

	uchar* data = ipmodul.pixelsUnmirror(N);
	if (data == nullptr)
		printf("unmirror unsuccessful\n");
	else
		printf("unmirror successful\n");
	
	//IPmodul::exportToPGM("../temp/unmirrorTest", vW->getImgWidth(), vW->getImgHeight(), 255, data);
	IPmodul::exportToPGM("../temp/mirrorTest", ipmodul.getImgWidth(), ipmodul.getImgHeight(), 255, ipmodul.getImgData(), false);
}
