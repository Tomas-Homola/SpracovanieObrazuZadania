#pragma once

#include <string>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned char uchar;
typedef unsigned int uint;

class ConvolutionKernel
{
private:
	int m_kernelSize = -1;
	double* m_pKernel = nullptr;

public:
	
	ConvolutionKernel() {}

	/// <summary>
	/// Creates a square convolution kernel that computes new values as the arithmetic mean from surrounding pixels.
	/// </summary>
	/// <param name="kernelSize">->...</param>
	ConvolutionKernel(int kernelSize) : m_kernelSize(kernelSize)
	{
		m_pKernel = static_cast<double*>(calloc(kernelSize * kernelSize, sizeof(double)));
		if (m_pKernel == nullptr)
			return;

		for (int i = 0; i < kernelSize*kernelSize; i++)
		{
			m_pKernel[i] = 1.0 / (kernelSize * kernelSize);
		}
	}

	
	~ConvolutionKernel() { free(m_pKernel); }

	int kernelSize() { return m_kernelSize; }
	double* kernel() { return m_pKernel; }

	bool setKernel(int kernelSize, double* kernel);

	void printKernel();
};

class IPmodul
{
private:
	//################# Variables #################//
	
	double* m_pImgLocalData = nullptr;

	uint m_imgWidth  = 0;
	uint m_imgHeight = 0;

	uint m_histogram[256] = { 0 };
	int m_minValue = -1;
	int m_maxValue = -1;

	double m_histogramNormalized[256] = { 0.0 };
	double m_histogramCumulative[256] = { 0.0 };
	//################# Methods #################//

	
public:
	/// <summary>
	/// Empty constructor
	/// </summary>
	IPmodul();

	/// <summary>
	/// Destructor
	/// </summary>
	~IPmodul();

	//################# Get functions #################//
	
	// Returns pointer to memory where the image data are stored.
	double* getImgData() { return m_pImgLocalData; }
	
	// Get width of the locally stored image.
	uint getImgWidth() { return m_imgWidth; }

	// Get height of the locally stored image.
	uint getImgHeight() { return m_imgHeight; }

	// Returns pointer to image histogram
	uint* getHistogram() { return m_histogram; }

	//################# Image Processing functions #################//

	/// <summary>
	/// Extend image by copying N pixels over its edges.
	/// </summary>
	/// <param name="originalImgData">-> input image data</param>
	/// <param name="bytesPerLine">-> ...</param>
	/// <param name="imgWidth">-> input image width</param>
	/// <param name="imgHeight">-> input image height</param>
	/// <param name="padding">-> number of pixels to mirror</param>
	/// <returns>True if successful, false otherwise.</returns>
	bool pixelsMirror(uchar* originalImgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const int padding);

	/// <summary>
	/// Crop extended image based on the given padding.
	/// </summary>
	/// <param name="padding">-> ...</param>
	/// <returns>Pointer to memory where cropped image data are stored if successful, nullptr otherwise.</returns>
	uchar* pixelsUnmirror(int padding);

	/// <summary>
	/// Computes histogram data for the given image: histogram, max and min values, normalized histogram and cumulative normalized histogram.
	/// </summary>
	/// <param name="originalImgData"></param>
	/// <param name="bytesPerLine"></param>
	/// <param name="imgWidth"></param>
	/// <param name="imgHeight"></param>
	void computeHistogramData(uchar* originalImgData, const int bytesPerLine, const int imgWidth, const int imgHeight);

	/// <summary>
	/// Performs Full Scale Histogram Stretch on the given image.
	/// </summary>
	/// <param name="imgData">-> ...</param>
	/// <param name="bytesPerLine">-> ...</param>
	/// <param name="imgWidth">-> ...</param>
	/// <param name="imgHeight">-> ...</param>
	/// <returns>True if successful, false otherwise.</returns>
	bool FSHS(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight);

	/// <summary>
	/// Performs Histogram Equalization on the given image.
	/// </summary>
	/// <param name="imgData"></param>
	/// <param name="bytesPerLine"></param>
	/// <param name="imgWidth"></param>
	/// <param name="imgHeight"></param>
	/// <returns>True if successful, false otherwise.</returns>
	bool EKV_HIST(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight);


	uchar* convolution(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, ConvolutionKernel* convolutionKernel);

	//################# Image Export functions #################//

	static bool exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, double* imgData, bool scaleData = true);
	static bool exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, uchar* imgData);

};