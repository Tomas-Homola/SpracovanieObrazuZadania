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


void IPmodul::updateEdges(const int padding)
{
	int temp = 0;
	int indexOld = 0;
	int indexNew = 0;

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
}

void IPmodul::filtrationPeronaMalik_UpdateEdges(const int padding)
{
	int temp = 0;
	int indexOld = 0;
	int indexNew = 0;

	// mirror over Upper and Lower edges
	temp = 1;
	for (int i = 0; i < padding; i++)
	{
		for (int j = padding; j < m_imgWidth - padding; j++)
		{
			// upper egde
			indexOld = (i + padding) * m_imgWidth + j;
			indexNew = (i + padding - temp) * m_imgWidth + j;
			m_uSigma[indexNew] = m_uSigma[indexOld];

			// lower edge
			indexOld = (m_imgHeight - i - padding - 1) * m_imgWidth + j;
			indexNew = (m_imgHeight - i - padding + temp - 1) * m_imgWidth + j;
			m_uSigma[indexNew] = m_uSigma[indexOld];
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
			m_uSigma[indexNew] = m_uSigma[indexOld];

			// right edge
			indexOld = i * m_imgWidth + (m_imgWidth - padding - 1 - j);
			indexNew = i * m_imgWidth + (m_imgWidth - padding - 1 - j + temp);
			m_uSigma[indexNew] = m_uSigma[indexOld];

			temp += 2;
		}
	}
}

void IPmodul::filtrationPeronaMalik_LinearDiffusion(double sigma)
{
	int padding = 1;
	int imgWidth = m_imgWidth - 2 * padding;
	int imgHeight = m_imgHeight - 2 * padding;

	double* b = new double[imgWidth * imgHeight] {0.0}; // right side of the system
	double* phi = m_uSigma; // = m_uSigma

	double newValue = 0.0;
	int indexC = -1, indexN = -1, indexS = -1, indexW = -1, indexE = -1;

	double temp = 0.0;

	// copy values from m_ImgLocalData to m_uSigma and vector b
	int i = 0, j = 0;
	for (int I = padding; I < m_imgHeight - padding; I++)
	{
		j = 0;
		for (int J = padding; J < m_imgWidth - padding; J++)
		{
			// copy internal data from m_pImgLocalData to m_uSigma
			indexC = I * m_imgWidth + J;
			m_uSigma[indexC] = m_pImgLocalData[indexC];

			// set right side to u0 = m_uSigma (original img data)
			b[i * imgWidth + j] = m_uSigma[indexC];

			j++;
		}
		i++;
	}

	// mirror edges to m_uSigma
	filtrationPeronaMalik_UpdateEdges(padding);

	// SOR variables
	double omega = 1.25;
	const int MAX_ITER = 1000;
	const double TOL = 1.0E-6;
	int iter = 0;
	double rez = 0.0;
	double sigmaSOR = 0.0; // SOR algorithm variable
	double Aii = 1.0 + 4 * sigma; // = 1 + 4*sigma
	double Aij = -sigma; // = -sigma

	// iterate through time steps
	for (int t = 0; t < 1; t++)
	{
		i = 0; j = 0;
		
		// tu sa asi bude diat SOR
		iter = 0;
		rez = 0.0;
		do
		{
			iter++;

			// iterate over all image pixels
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;
					indexN = (I - 1) * m_imgWidth + J;
					indexS = (I + 1) * m_imgWidth + J;
					indexW = I * m_imgWidth + J - 1;
					indexE = I * m_imgWidth + J + 1;

					sigmaSOR = 0.0;
					sigmaSOR += Aij * (phi[indexN] + phi[indexS] + phi[indexE] + phi[indexW]);

					phi[indexC] = (1.0 - omega) * phi[indexC] + (omega / Aii) * (b[i * imgWidth + j] - sigmaSOR);

					j++;
				}
				i++;
			}

			// update edges to m_uSigma
			filtrationPeronaMalik_UpdateEdges(padding);

			// compute residuals

			rez = 0.0;
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;
					indexN = (I - 1) * m_imgWidth + J;
					indexS = (I + 1) * m_imgWidth + J;
					indexW = I * m_imgWidth + J - 1;
					indexE = I * m_imgWidth + J + 1;

					// compute Ax - b
					temp = (Aii * phi[indexC] + Aij * (phi[indexN] + phi[indexS] + phi[indexE] + phi[indexW])) - b[i * imgWidth + j];
					rez += temp * temp;

					j++;
				}
				i++;
			}
			rez = sqrt(rez);
			//printf("time %d -> SOR iter %d: rez = %.8lf\n", t, iter, rez);
			if (rez < TOL) // if residuals small enough -> break
				break;

		} while (iter <= MAX_ITER); // END SOR

		//printf("PM linear filtration -> SOR stop iter %d: rez = %.8lf\n", iter, rez);
	}

	// free memory from vector b
	delete[] b;
}

void IPmodul::filtrationPeronaMalik_ComputeGradientsNormSquared(int padding)
{
	double gradX = 0.0;
	double gradY = 0.0;

	double h = 1.0;

	int NW = 0, N = 0, NE = 0;
	int W  = 0, p = 0, E  = 0;
	int SW = 0, S = 0, SE = 0;

	// iterate over all internal img pixels of m_uSigma
	for (int I = padding; I < m_imgHeight - padding; I++)
	{
		for (int J = padding; J < m_imgWidth - padding; J++)
		{
			// compute correct indices
			NW = (I - 1) * m_imgWidth + (J - 1);
			N  = (I - 1) * m_imgWidth + J;
			NE = (I - 1) * m_imgWidth + (J + 1);

			W = I * m_imgWidth + (J - 1);
			p = I * m_imgWidth + J;
			E = I * m_imgWidth + (J + 1);

			SW = (I + 1) * m_imgWidth + (J - 1);
			S  = (I + 1) * m_imgWidth + J;
			SE = (I + 1) * m_imgWidth + (J + 1);

			// EAST edge
			gradX = (m_uSigma[E] - m_uSigma[p]) / h;
			gradY = (m_uSigma[N] + m_uSigma[NE] - m_uSigma[S] - m_uSigma[SE]) / (4.0 * h);

			m_gradientsNormSquared[p].gradSqEast = gradX * gradX + gradY * gradY;

			// NORTH edge
			gradX = (m_uSigma[W] + m_uSigma[NW] - m_uSigma[E] - m_uSigma[NE]) / (4.0 * h);
			gradY = (m_uSigma[N] - m_uSigma[p]) / h;

			m_gradientsNormSquared[p].gradSqNorth = gradX * gradX + gradY * gradY;

			// WEST edge
			gradX = (m_uSigma[W] - m_uSigma[p]) / h;
			gradY = (m_uSigma[S] + m_uSigma[SW] - m_uSigma[N] - m_uSigma[NW]) / (4.0 * h);

			m_gradientsNormSquared[p].gradSqWest = gradX * gradX + gradY * gradY;

			// SOUTH edge
			gradX = (m_uSigma[E] + m_uSigma[SE] - m_uSigma[W] - m_uSigma[SW]) / (4.0 * h);
			gradY = (m_uSigma[S] - m_uSigma[p]) / h;

			m_gradientsNormSquared[p].gradSqSouth = gradX * gradX + gradY * gradY;

		}
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

bool IPmodul::pixelsMirrorDouble(double* originalImgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const int padding)
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

			m_pImgLocalData[indexNew] = originalImgData[indexOld];
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

uchar* IPmodul::filtrationExplicitHeatEq(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double tau, const double h, const int timeSteps)
{
	// mirror pixels
	int padding = 1;
	pixelsMirror(imgData, bytesPerLine, imgWidth, imgHeight, padding);

	size_t size = (size_t)imgWidth * imgHeight;
	uchar* resultData = new uchar[size]{ 0 };
	double* newData = new double[size] {0.0};
	
	double newValue = 0.0;
	uchar scaledValue = 0;
	int indexC = -1, indexN = -1, indexS = -1, indexW = -1, indexE = -1, indexNew = -1;
	int iNew = 0, jNew = 0;

	double sum = 0.0;
	double scaledValueD = 0.0;
	uchar origValue = 0;

	// compute mean value of the original image
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			origValue = imgData[i * bytesPerLine + j];
			scaledValueD = static_cast<double>(origValue);
			sum += scaledValueD;
		}
	}

	if (printMsg)
		printf("Original image mean value: %.10lf\n", sum / size);

	// iterate through time steps
	for (int t = 0; t < timeSteps; t++)
	{
		iNew = 0; jNew = 0;
		sum = 0;

		// iterate through extended image
		for (int i = padding; i < m_imgHeight - padding; i++)
		{
			for (int j = padding; j < m_imgWidth - padding; j++)
			{
				indexC = i * m_imgWidth + j;
				indexN = (i - 1) * m_imgWidth + j;
				indexS = (i + 1) * m_imgWidth + j;
				indexW = i * m_imgWidth + j - 1;
				indexE = i * m_imgWidth + j + 1;

				newValue = (1.0 - 4.0 * ((tau) / (h * h))) * m_pImgLocalData[indexC] + (tau / (h * h)) * (m_pImgLocalData[indexN] + m_pImgLocalData[indexS] + m_pImgLocalData[indexW] + m_pImgLocalData[indexE]);

				indexNew = iNew * imgWidth + jNew;
				newData[indexNew] = newValue;

				jNew++;
			}
			iNew++;
			jNew = 0;
		}

		// compute mean value of the image
		for (size_t i = 0; i < size; i++)
		{
			sum += newData[i];
		}

		if (printMsg)
			printf("Time step %d filtered image EXPLICIT mean value: %.10lf\n", t, sum / size);

		if (t != (timeSteps - 1))
			pixelsMirrorDouble(newData, imgWidth, imgWidth, imgHeight, padding);
	}

	if (printMsg)
		printf("\n\n");
	
	// cast double data to uchar
	for (size_t i = 0; i < size; i++)
	{
		resultData[i] = static_cast<uchar>(newData[i] + 0.5);
	}

	// release memory
	delete[] newData;

	return resultData;
}


uchar* IPmodul::filtrationImplicitHeatEq(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double tau, const int timeSteps)
{
	// mirror pixels
	int padding = 1;
	pixelsMirror(imgData, bytesPerLine, imgWidth, imgHeight, padding);

	size_t size = (size_t)imgWidth * imgHeight; // size of the original image
	uchar* resultData = nullptr; // return data in uchar after filtration
	double* b = new double[size] {0.0}; // right side of the system
	double* phi = m_pImgLocalData; // = m_pImgLocalData

	double newValue = 0.0;
	int indexC = -1, indexN = -1, indexS = -1, indexW = -1, indexE = -1, indexNew = -1;

	double sum = 0.0;
	double temp = 0.0;
	uchar origValue = 0;

	// SOR
	double omega = 1.25;
	const int MAX_ITER = 1000;
	const double TOL = 1.0E-6;
	int iter = 0;
	double rez = 0.0;
	double sigmaSOR = 0.0; // SOR algorithm variable
	double Aii = 1.0 + 4 * tau; // = 1 + 4*tau
	double Aij = -tau; // = -tau
	
	// compute mean value of the original image
	// TODO: prepisat na cyklus cez m_pLocalImgData, nech sa nemusi robit znovu static_cast
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			origValue = imgData[i * bytesPerLine + j];
			temp = static_cast<double>(origValue);
			sum += temp;

			// set right side to u0 = original image data
			b[i * imgWidth + j] = temp;
		}
	}

	if (printMsg)
		printf("Original image mean value: %.10lf\n", sum / size);

	int i = 0, j = 0;
	// iterate through time steps
	for (int t = 0; t < timeSteps; t++)
	{
		i = 0; j = 0;
		sum = 0;

		// tu sa asi bude diat SOR
		iter = 0;
		rez = 0.0;
		do
		{
			iter++;

			// iterate over all image pixels
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;
					indexN = (I - 1) * m_imgWidth + J;
					indexS = (I + 1) * m_imgWidth + J;
					indexW = I * m_imgWidth + J - 1;
					indexE = I * m_imgWidth + J + 1;

					sigmaSOR = 0.0;
					sigmaSOR += Aij * (phi[indexN] + phi[indexS] + phi[indexE] + phi[indexW]);

					phi[indexC] = (1.0 - omega) * phi[indexC] + (omega / Aii) * (b[i * imgWidth + j] - sigmaSOR);

					j++;
				}
				i++;
			}

			// compute residuals
			updateEdges(padding);

			rez = 0.0;
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;
					indexN = (I - 1) * m_imgWidth + J;
					indexS = (I + 1) * m_imgWidth + J;
					indexW = I * m_imgWidth + J - 1;
					indexE = I * m_imgWidth + J + 1;

					// compute Ax - b
					temp = (Aii * phi[indexC] + Aij * (phi[indexN] + phi[indexS] + phi[indexE] + phi[indexW])) - b[i * imgWidth + j];
					rez += temp * temp;

					j++;
				}
				i++;
			}
			rez = sqrt(rez);
			//printf("time %d -> SOR iter %d: rez = %.8lf\n", t, iter, rez);
			if (rez < TOL)
				break;

		} while (iter <= MAX_ITER); // END SOR

		printf("time %d -> SOR stop iter %d: rez = %.8lf\n", t, iter, rez);
		
		// compute mean value of the image
		sum = 0.0;
		i = 0;
		for (int I = padding; I < m_imgHeight - padding; I++)
		{
			j = 0;
			for (int J = padding; J < m_imgWidth - padding; J++)
			{
				indexC = I * m_imgWidth + J;
				sum += phi[indexC];

				// update b vector as the solution from the last time iteration
				b[i * imgWidth + j] = phi[indexC];
				j++;
			}
			i++;
		}
		
		if (printMsg)
			printf("Time step %d filtered image IMPLICIT mean value: %.10lf\n", t, sum / size);

	}
	if (printMsg)
		printf("\n\n");
	
	// free memory from vector b
	delete[] b;

	resultData = pixelsUnmirror(padding);
	return resultData;
}

uchar* IPmodul::filtrationSemiImplicitPeronaMalik(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double sigma, const double tau, const double K, const int timeSteps)
{
	// mirror pixels
	int padding = 1;
	pixelsMirror(imgData, bytesPerLine, imgWidth, imgHeight, padding);

	size_t size = (size_t)imgWidth * imgHeight; // size of the original image
	uchar* resultData = nullptr; // return data in uchar after filtration

	double* b = new double[size] {0.0}; // right side of the system
	double* phi = m_pImgLocalData; // = m_pImgLocalData
	m_uSigma = new double[(size_t)m_imgWidth * m_imgHeight] {0.0}; // allocate space for uSigma, same as m_pImgLocalData

	// resize vector for gradient norms squared
	m_gradientsNormSquared.resize((size_t)m_imgWidth * m_imgHeight, GradientNormSquared());

	// indices
	int indexC = -1, indexN = -1, indexS = -1, indexW = -1, indexE = -1, indexNew = -1;

	// auxilliary variables
	double newValue = 0.0;
	double sum = 0.0;
	double temp = 0.0;
	uchar origValue = 0;

	// system matrix coefficients
	double g_N = 0.0; // diffusion ceofficient over North edge
	double g_S = 0.0; // diffusion ceofficient over South edge
	double g_E = 0.0; // diffusion ceofficient over East  edge
	double g_W = 0.0; // diffusion ceofficient over West  edge

	double Aii   = 0.0; // = 1 + tau * (g_N + g_S + g_E + g_W)
	double Aij_N = 0.0; // = - tau * g_N
	double Aij_S = 0.0; // = - tau * g_S
	double Aij_E = 0.0; // = - tau * g_E
	double Aij_W = 0.0; // = - tau * g_W

	// compute mean value of the original image
	// TODO: prepisat na cyklus cez m_pLocalImgData, nech sa nemusi robit znovu static_cast
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			origValue = imgData[i * bytesPerLine + j];
			temp = static_cast<double>(origValue);
			sum += temp;

			// set right side to u0 = original image data
			b[i * imgWidth + j] = temp;
		}
	}

	if (printMsg)
		printf("Original image mean value: %.12lf\n", sum / size);


	// SOR variables
	double omega = 1.25;
	const int MAX_ITER = 1000;
	const double TOL = 1.0E-6;
	int iter = 0;
	double rez = 0.0;
	double sigmaSOR = 0.0; // SOR algorithm sigma variable

	int i = 0, j = 0;
	// iterate through time steps
	for (int t = 1; t <= timeSteps; t++)
	{
		i = 0; j = 0;
		sum = 0;

		// perform one time step of linear diffusion with step sigma
		filtrationPeronaMalik_LinearDiffusion(sigma);

		// compute gradient norms squared from m_uSigma
		filtrationPeronaMalik_ComputeGradientsNormSquared(padding);

		// tu sa asi bude diat SOR
		iter = 0;
		rez = 0.0;
		do
		{
			iter++;

			// iterate over all image pixels
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;
					indexN = (I - 1) * m_imgWidth + J;
					indexS = (I + 1) * m_imgWidth + J;
					indexW = I * m_imgWidth + J - 1;
					indexE = I * m_imgWidth + J + 1;

					// compute system matrix coefficients
					g_N = diffCoefFunction1(K, m_gradientsNormSquared[indexC].gradSqNorth);
					g_S = diffCoefFunction1(K, m_gradientsNormSquared[indexC].gradSqSouth);
					g_E = diffCoefFunction1(K, m_gradientsNormSquared[indexC].gradSqEast);
					g_W = diffCoefFunction1(K, m_gradientsNormSquared[indexC].gradSqWest);

					Aii   = 1 + tau * (g_N + g_S + g_E + g_W);
					Aij_N = -tau * g_N;
					Aij_S = -tau * g_S;
					Aij_E = -tau * g_E;
					Aij_W = -tau * g_W;

					sigmaSOR = Aij_N * phi[indexN] + Aij_S * phi[indexS] + Aij_E * phi[indexE] + Aij_W * phi[indexW];

					phi[indexC] = (1.0 - omega) * phi[indexC] + (omega / Aii) * (b[i * imgWidth + j] - sigmaSOR);

					j++;
				}
				i++;
			}

			// compute residuals
			updateEdges(padding);

			rez = 0.0;
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;
					indexN = (I - 1) * m_imgWidth + J;
					indexS = (I + 1) * m_imgWidth + J;
					indexW = I * m_imgWidth + J - 1;
					indexE = I * m_imgWidth + J + 1;

					// compute system matrix coefficients
					g_N = diffCoefFunction1(K, m_gradientsNormSquared[indexC].gradSqNorth);
					g_S = diffCoefFunction1(K, m_gradientsNormSquared[indexC].gradSqSouth);
					g_E = diffCoefFunction1(K, m_gradientsNormSquared[indexC].gradSqEast);
					g_W = diffCoefFunction1(K, m_gradientsNormSquared[indexC].gradSqWest);

					Aii = 1 + tau * (g_N + g_S + g_E + g_W);
					Aij_N = -tau * g_N;
					Aij_S = -tau * g_S;
					Aij_E = -tau * g_E;
					Aij_W = -tau * g_W;

					// compute Ax - b
					temp = (Aii * phi[indexC] + Aij_N * phi[indexN] + Aij_S * phi[indexS] + Aij_E * phi[indexE] + Aij_W * phi[indexW]) - b[i * imgWidth + j];
					rez += temp * temp;

					j++;
				}
				i++;
			}
			rez = sqrt(rez);
			//printf("time %d -> SOR iter %d: rez = %.8lf\n", t, iter, rez);
			if (rez < TOL)
				break;

		} while (iter <= MAX_ITER); // END SOR

		printf("time %d -> SOR stop iter %d: rez = %.8lf\n", t, iter, rez);

		// compute mean value of the image
		sum = 0.0;
		i = 0;
		for (int I = padding; I < m_imgHeight - padding; I++)
		{
			j = 0;
			for (int J = padding; J < m_imgWidth - padding; J++)
			{
				indexC = I * m_imgWidth + J;
				sum += phi[indexC];

				// update b vector as the solution from the last time iteration
				b[i * imgWidth + j] = phi[indexC];
				j++;
			}
			i++;
		}

		if (printMsg)
			printf("Time step %d filtered image PERONA-MALIK mean value: %.12lf\n", t, sum / size);

	}
	if (printMsg)
		printf("\n\n");

	delete[] b; // clear space for vector b
	delete[] m_uSigma; // clear space for m_uSigma
	m_gradientsNormSquared.clear(); // clear vector for gradient norms squared

	// return filtered data
	resultData = pixelsUnmirror(padding);
	return resultData;
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

			if ((i + 1) % imgWidth == 0)
			{
				fprintf(fp, "\n");
			}

			if ((i + 1) % (dataSize / 10) == 0)
				printf("\rExporting image to pgm... %d%% done", 10 * ((int)i + 1) / (dataSize / 10));
		}
	}
	else
	{
		for (size_t i = 0; i < imgHeight; i++)
		{
			for (size_t j = 0; j < imgWidth; j++)
			{
				scaledValue = static_cast<unsigned char>(imgData[i * imgWidth + j] + 0.5);
				fprintf(fp, "%d ", scaledValue);
			}
			fprintf(fp, "\n");
		}

		/*for (size_t i = 0; i < dataSize; i++)
		{
			scaledValue = static_cast<unsigned char>(imgData[i] + 0.5);
			fprintf(fp, "%d ", scaledValue);

			if ((i + 1) % imgWidth == 0)
			{
				fprintf(fp, "\n");
				printf("i+1: %d\n",i);
			}

			if ((i + 1) % (dataSize / 10) == 0)
				printf("\rExporting image to pgm... %d%% done", 10 * ((int)i + 1) / (dataSize / 10));
		}*/
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
			printf("\rExporting image to ppm... %llu%% done", 10 * (i + 1) / (dataSize / 10));
	}
	fclose(fp);

	return true;
}