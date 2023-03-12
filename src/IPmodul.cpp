#include "IPmodul.h"

bool ConvolutionKernel::setKernel(int kernelSize, double* kernel)
{
	m_kernelSize = kernelSize;

	if (m_pKernel != nullptr) // delete old kernel data
	{
		delete[] m_pKernel; // delete old kernel
	}

	// allocate space for new kernel
	m_pKernel = new double[kernelSize * kernelSize] {0.0};
	//m_pKernel = static_cast<double*>(calloc(kernelSize * kernelSize, sizeof(double)));
	if (m_pKernel == nullptr) // check if allocation was successful
		return false;

	// copy new kernel data
	for (int i = 0; i < kernelSize * kernelSize; i++)
	{
		m_pKernel[i] = kernel[i];
	}

	return true;
}

void ConvolutionKernel::printKernel()
{
	for (int i = 0; i < m_kernelSize; i++)
	{
		for (int j = 0; j < m_kernelSize; j++)
		{
			printf("%.4f\t", m_pKernel[i * m_kernelSize + j]);
		}
		printf("\n");
	}
}


IPmodul::IPmodul()
{

}


IPmodul::~IPmodul()
{
	delete[] m_pImgLocalData;
}

bool IPmodul::pixelsMirror(uchar* originalImgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const int padding)
{
	int indexNew = 0, indexOld = 0;
	int temp = 0;

	// check if there is already some image data stored
	if (m_pImgLocalData != nullptr)
	{
		delete[] m_pImgLocalData; // if there is, delete old data
		m_pImgLocalData = nullptr;
	}
	
	// compute new size
	m_imgWidth = imgWidth + 2 * padding;
	m_imgHeight = imgHeight + 2 * padding;
	size_t size = (size_t)m_imgWidth * m_imgHeight;
	
	m_pImgLocalData = new double[size] {0.0}; // allocate memory
	//m_pImgLocalData = (double*)calloc(size, sizeof(double)); // allocate memory

	if (m_pImgLocalData == nullptr) // check, if allocation was successful
		return false;

	// copy old image data
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			indexNew = (i + padding) * m_imgWidth + (j + padding);
			indexOld = i * bytesPerLine + j;

			m_pImgLocalData[indexNew] = static_cast<double>(originalImgData[indexOld]);
		}
	}

	// mirror over Upper and Lower edges
	temp = 1;
	for (int i = 0; i < padding; i++)
	{
		for (int j = padding; j < m_imgWidth - padding; j++)
		{
			// upper egde
			indexOld = (i + padding) * m_imgWidth + j;
			indexNew = (i + padding - temp) * m_imgWidth + j;
			m_pImgLocalData[indexNew] = m_pImgLocalData[indexOld];

			// lower edge
			indexOld = (m_imgHeight - i - padding - 1) * m_imgWidth + j;
			indexNew = (m_imgHeight - i - padding + temp - 1) * m_imgWidth + j;
			m_pImgLocalData[indexNew] = m_pImgLocalData[indexOld];
		}
		temp += 2;
	}

	// mirror over Left and Right edges
	for (int i = 0; i < m_imgHeight; i++)
	{
		temp = 1;
		for (int j = 0; j < padding; j++)
		{
			// left edge
			indexOld = i * m_imgWidth + (j + padding);
			indexNew = i * m_imgWidth + (j + padding - temp);
			m_pImgLocalData[indexNew] = m_pImgLocalData[indexOld];

			// right edge
			indexOld = i * m_imgWidth + (m_imgWidth - padding - 1 - j);
			indexNew = i * m_imgWidth + (m_imgWidth - padding - 1 - j + temp);
			m_pImgLocalData[indexNew] = m_pImgLocalData[indexOld];

			temp += 2;
		}
	}

	// print data
	//for (int i = 0; i < m_imgHeight; i++)
	//{
	//	for (int j = 0; j < m_imgWidth; j++)
	//	{
	//		indexNew = i * m_imgWidth + j;
	//
	//		printf("%.0f\t", m_pImgLocalData[indexNew]);
	//	}
	//	printf("\n");
	//}

	return true;
}

uchar* IPmodul::pixelsUnmirror(int padding)
{
	if (m_pImgLocalData == nullptr) // check, if there is some image to be cropped
		return nullptr;

	// calculate new size
	uint newWidth  = m_imgWidth - 2 * padding;
	uint newHeight = m_imgHeight - 2 * padding;
	size_t size = (size_t)newWidth * newHeight;

	uchar* pImgData = new uchar[size]{ 0 };
	//uchar* pImgData = (uchar*)calloc(size, sizeof(uchar)); // allocate new memory
	
	if (pImgData == nullptr) // check, if allocation successful
		return nullptr;

	// copy only the requested data based on the "padding"
	int indexOld = 0, indexNew = 0;
	int iNew = 0, jNew = 0;
	for (int i = padding; i < m_imgHeight - padding; i++)
	{
		jNew = 0;
		for (int j = padding; j < m_imgWidth - padding; j++)
		{
			indexOld = i * m_imgWidth + j;
			indexNew = iNew * newWidth + jNew;
			pImgData[indexNew] = static_cast<uchar>(m_pImgLocalData[indexOld] + 0.5);
			
			jNew++;
		}
		iNew++;
	}

	return pImgData;
}

void IPmodul::computeHistogramData(uchar* originalImgData, const int bytesPerLine, const int imgWidth, const int imgHeight)
{
	// remove old histogram values
	for (int i = 0; i < 256; i++)
	{
		m_histogram[i] = 0;
	}

	int index = 0;
	uchar value = 0;
	m_minValue = INT_MAX;
	m_maxValue = INT_MIN;

	// compute image histogram
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			// read pixel value
			index = i * bytesPerLine + j;
			value = originalImgData[index];

			// increase corresponding place in histogram
			m_histogram[value]++;

			// find min/max values of the given image
			if (value < m_minValue) m_minValue = value;
			if (value > m_maxValue) m_maxValue = value;

		}
	}

	int pixels = imgWidth * imgHeight;
	// compute normalized histogram
	for (int i = 0; i < 256; i++)
	{
		// compute relative frequency
		m_histogramNormalized[i] = static_cast<double>(m_histogram[i]) / pixels;
	}

	// compute cumulative normalized histogram
	m_histogramCumulative[0] = m_histogramNormalized[0];
	for (int i = 1; i < 256; i++)
	{
		m_histogramCumulative[i] = m_histogramCumulative[i - 1] + m_histogramNormalized[i];
	}
}

bool IPmodul::FSHS(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight)
{
	if (imgData == nullptr)
		return false;

	// compute histogram for the given image, function also finds min and max values
	computeHistogramData(imgData, bytesPerLine, imgWidth, imgHeight);

	int index = 0;
	double temp = 0;
	uchar scaledValue = 0;
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			index = i * bytesPerLine + j;
			// scale values from original range [m_minValue, m_maxValue] to [0, 255]
			temp = static_cast<double>(imgData[index]);
			scaledValue = static_cast<uchar>(((temp - m_minValue) / (m_maxValue - m_minValue)) * 255 + 0.5);
			imgData[index] = scaledValue;
		}
	}

	return true;
}

bool IPmodul::EKV_HIST(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight)
{
	if (imgData == nullptr)
		return false;

	// compute histogram data for the given image
	computeHistogramData(imgData, bytesPerLine, imgWidth, imgHeight);

	int index = 0;
	int temp = 0;
	uchar scaledValue = 0;
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			index = i * bytesPerLine + j;
			
			temp = static_cast<int>(imgData[index]); // get pixel value
			scaledValue = static_cast<uchar>(255.0 * m_histogramCumulative[temp] + 0.5);
			
			imgData[index] = scaledValue;
		}
	}

	return true;
}

uchar* IPmodul::convolution(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, ConvolutionKernel* convolutionKernel)
{
	uchar* newImgData = nullptr;
	
	double* kernel = nullptr;
	int kernelSize = -1;

	double newValue = 0.0;
	double imgValue = 0.0;
	int imgIndex = -1;
	int imgNewIndex = -1;
	int kernelIndex = -1;
	int iNew = 0, jNew = 0;
	uchar scaledValue = 0;

	// find the necessary padding for pixels mirroring
	int padding = convolutionKernel->kernelSize() / 2;

	// extend image data based on kernel size
	if (!pixelsMirror(imgData, bytesPerLine, imgWidth, imgHeight, padding))
		return nullptr;

	// allocate space for new image data
	newImgData = new uchar[(size_t)imgWidth * imgHeight]{ 0 };
	//newImgData = static_cast<uchar*>(calloc(imgWidth * imgHeight, sizeof(uchar*)));
	if (newImgData == nullptr) // check, if allocation successful
		return nullptr;

	// save kernel info
	kernel     = convolutionKernel->kernel();
	kernelSize = convolutionKernel->kernelSize();

	// perform convolution with given kernel
	// original image saved in m_pImgLocalData at [padding, imgHeight - padding]x[padding, imgWidth - padding]
	// iterate through all pixels of the original data
	for (int i = padding; i < m_imgHeight - padding; i++)
	{
		jNew = 0;
		for (int j = padding; j < m_imgWidth - padding; j++)
		{
			newValue = 0.0;

			// at (i,j) iterate through all kernel elements
			for (int k = -padding; k <= padding; k++)
			{
				for (int l = -padding; l <= padding; l++)
				{
					// compute corresponding indices for image and kernel data
					imgIndex    = (i + k) * m_imgWidth + (j + l);
					kernelIndex = (k + padding) * kernelSize + (l + padding);

					// compute new value by multipling pixel data with corresponding kernel weight
					imgValue  = m_pImgLocalData[imgIndex];
					newValue += imgValue * kernel[kernelIndex];
				}
			} // end of kernel iterations

			// compute index in new image
			imgNewIndex = iNew * imgWidth + jNew;

			// check, wheather newValue is not out of bounds [0, 255]
			if (newValue > 255.0) newValue = 255.0;
			if (newValue < 0.0)   newValue = 0.0;

			scaledValue = static_cast<uchar>(newValue + 0.5);
			newImgData[imgNewIndex] = scaledValue;

			jNew++;
		}
		
		iNew++;
	}

	// return new image data
	return newImgData;
}

bool IPmodul::exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, double* imgData, bool scaleData)
{
	printf("Exporting image to pgm...");
	FILE* fp = nullptr;
	fp = fopen((fileName + ".pgm").c_str(), "w+");
	if (fp == nullptr)
		return false;

	unsigned char scaledValue = 0;
	int dataSize = imgWidth * imgHeight;
	fprintf(fp, "P2\n%d %d\n%d\n", imgWidth, imgHeight, maxValue);
	
	if (scaleData)
	{
		for (size_t i = 0; i < dataSize; i++)
		{
			scaledValue = static_cast<unsigned char>(imgData[i] * maxValue + 0.5);
			fprintf(fp, "%d ", scaledValue);

			if ((i + 1) % 70 == 0)
				fprintf(fp, "\n");

			if ((i + 1) % (dataSize / 10) == 0)
				printf("\rExporting image to pgm... %d%% done", 10 * ((int)i + 1) / (dataSize / 10));
		}
	}
	else
	{
		for (size_t i = 0; i < dataSize; i++)
		{
			scaledValue = static_cast<unsigned char>(imgData[i] + 0.5);
			fprintf(fp, "%d ", scaledValue);

			if ((i + 1) % 70 == 0)
				fprintf(fp, "\n");

			if ((i + 1) % (dataSize / 10) == 0)
				printf("\rExporting image to pgm... %d%% done", 10 * ((int)i + 1) / (dataSize / 10));
		}
	}
	printf("\n");
	fclose(fp);

	return true;
}

bool IPmodul::exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, uchar* imgData)
{
	printf("Exporting image to pgm...");
	FILE* fp = nullptr;
	fp = fopen((fileName + ".pgm").c_str(), "w+");
	if (fp == nullptr)
		return false;

	int dataSize = imgWidth * imgHeight;
	fprintf(fp, "P2\n%d %d\n%d\n", imgWidth, imgHeight, maxValue);
	for (size_t i = 0; i < dataSize; i++)
	{
		fprintf(fp, "%d ", imgData[i]);

		if ((i + 1) % 70 == 0)
			fprintf(fp, "\n");

		if ((i + 1) % (dataSize / 10) == 0)
		{
			printf("\rExporting image to pgm... %d%% done", 10 * ((int)i + 1) / (dataSize / 10));
			//_sleep(500);
		}
	}
	printf("\n");
	fclose(fp);

	return true;
}

bool IPmodul::ExportToPPM(std::string fileName, int width, int height, int maxValue, float* r, float* g, float* b)
{
	printf("Exporting image to ppm...\n");
	FILE* fp = nullptr;
	fp = fopen((fileName + ".ppm").c_str(), "w+");
	if (fp == nullptr)
		return false;

	unsigned char scaledR = 0, scaledG = 0, scaledB = 0;
	size_t dataSize = width * height;
	fprintf(fp, "P3\n%d %d\n%d\n", width, height, maxValue);
	for (size_t i = 0; i < dataSize; i++)
	{
		scaledR = static_cast<unsigned char>(r[i] * maxValue + 0.5);
		scaledG = static_cast<unsigned char>(g[i] * maxValue + 0.5);
		scaledB = static_cast<unsigned char>(b[i] * maxValue + 0.5);

		fprintf(fp, "%d %d %d\t", scaledR, scaledG, scaledB);

		if ((i + 1) % 70 == 0)
			fprintf(fp, "\n");

		if ((i + 1) % (dataSize / 10) == 0)
			printf("\rExporting image to ppm... %d%% done", 10 * (i + 1) / (dataSize / 10));
	}
	fclose(fp);

	return true;
}