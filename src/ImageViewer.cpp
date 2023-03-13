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
			vW->getData()[i * vW->getImage()->bytesPerLine() + j] = newImg[i * vW->getImgWidth() + j];
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

	uchar* newImg = ip.filtrationExplicitHeatEq(vW->getData(), vW->getBytesPerLine(), vW->getImgWidth(), vW->getImgHeight(), 0.15, 1.0, 10);

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

void ImageViewer::on_pushButton_mirrorTest_clicked()
{
	if (vW->getImage() == nullptr)
		return;
	
	IPmodul ipmodul;
	uint N = QInputDialog::getInt(this, tr("Mirror pixels"), tr("Number of pixels:"), 2, 1);
	
	//ipmodul.pixelsMirror(d, 4, 4, N);
	bool result = ipmodul.pixelsMirror(vW->getData(), vW->getImage()->bytesPerLine(), vW->getImage()->width(), vW->getImage()->height(), N);

	printf("mirror pixels result: %d\n", result);
	
	uchar* data = ipmodul.pixelsUnmirror(N);
	if (data == nullptr)
		printf("unmirror unsuccessful\n");
	else
		printf("unmirror successful\n");
	
	//IPmodul::exportToPGM("../temp/test", vW->getImgWidth(), vW->getImgHeight(), 255, data);
	IPmodul::exportToPGM("../temp/test2", ipmodul.getImgWidth(), ipmodul.getImgHeight(), 255, ipmodul.getImgData(), false);
}
