#pragma once

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <omp.h>

#define EPSILON 0.000000001

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
		//m_pKernel = static_cast<double*>(calloc(kernelSize * kernelSize, sizeof(double)));
		m_pKernel = new double[(size_t)kernelSize * kernelSize] {0.0};
		if (m_pKernel == nullptr)
			return;

		for (int i = 0; i < kernelSize*kernelSize; i++)
		{
			m_pKernel[i] = 1.0 / (kernelSize * kernelSize);
		}
	}

	
	~ConvolutionKernel() { delete[] m_pKernel; }

	int kernelSize() { return m_kernelSize; }
	double* kernel() { return m_pKernel; }

	bool setKernel(int kernelSize, double* kernel);

	void printKernel();
};

inline double diffCoefFunction1(double K, double normGradSquared) { return 1.0 / (1.0 + K * normGradSquared); }

class IPmodul
{
private:

	//################# Structs #################//

	struct GradientNormSquared
	{
		double gradSqNorth = 1.0;
		double gradSqSouth = 1.0;
		double gradSqEast = 1.0;
		double gradSqWest = 1.0;
	};

	struct AllGradientsNormSquared
	{
		double uSigmaGradSqNorth = 1.0;
		double uSigmaGradSqSouth = 1.0;
		double uSigmaGradSqEast  = 1.0;
		double uSigmaGradSqWest  = 1.0;

		double uOrigGradNorth = 1.0;
		double uOrigGradSouth = 1.0;
		double uOrigGradEast = 1.0;
		double uOrigGradWest = 1.0;

		double uMeanGrad = 1.0;
	};

	struct MatrixCoefs
	{
		double Aii = 1.0;
		double Aij_N = 0.0;
		double Aij_S = 0.0;
		double Aij_E = 0.0;
		double Aij_W = 0.0;
	};

	//################# Variables #################//
	
	int m_threads = 2;

	double* m_pImgLocalData = nullptr;
	double* m_uSigma = nullptr;

	uint m_imgWidth  = 0;
	uint m_imgHeight = 0;

	uint m_histogram[256] = { 0 };
	int m_minValue = -1;
	int m_maxValue = -1;

	double m_histogramNormalized[256] = { 0.0 };
	double m_histogramCumulative[256] = { 0.0 };

	// Vector of system matrix coefficients
	std::vector<MatrixCoefs> m_matrixCoefs = {};

	double m_MCF_epsilon = 0.01; // 0.01

	//################# Methods #################//
	
	// update local image data edges
	void updateEdges(const int padding);

	void updateEdges(const int padding, double* data);

	//------------------Perona-Malik filtration---------------------

	// update edges for m_uSigma
	void uSigma_UpdateEdges(const int padding);

	// Performs 1 linear diffusion step for Perona-Malik with parameter sigma, stores values to m_uSigma.
	void uSigma_LinearDiffusion(double sigma);

	// compute system matrix coefficients for each pixel of the image for Perona-Malik model
	void PM_computeMatrixCoefs(int padding, double tau, double K);

	//------------------MCF/GMCF filtration---------------------

	// compute system matrix coefficients for each pixel of the image
	void GMCF_computeMatrixCoefs(int padding, double epsilon, double tau, double K);

	// compute v = A.p in BiCGStab
	void BiCGStab_compute_v_Ap(int imgWidth, int imgHeight, int padding, double* v, double* p);

	// compute t = A.s in BiCGStab
	void BiCGStab_compute_t_As(int imgWidth, int imgHeight, int padding, double* t, double* s);

	//------------------Distance function---------------------

	//
	double M(double* d, int i, int j, int p, int q);

	//
	void SignDistFun_setData(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, int* Fn, int* inOutPixels, int* toUpdate);
public:

	// Print messages from filtration functions.
	bool printMsg = true;

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

	// Returns pointer to normalized image histogram
	double* getHistogramNormalized() { return m_histogramNormalized; }

	// Returns pointer to normalized image histogram
	double* getHistogramCummulative() { return m_histogramCumulative; }

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
	/// Extend image by copying N pixels over its edges.
	/// </summary>
	/// <param name="originalImgData">-> input image data</param>
	/// <param name="bytesPerLine">-> ...</param>
	/// <param name="imgWidth">-> input image width</param>
	/// <param name="imgHeight">-> input image height</param>
	/// <param name="padding">-> number of pixels to mirror</param>
	/// <returns>True if successful, false otherwise.</returns>
	bool pixelsMirrorDouble(double* originalImgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const int padding);

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

	/// <summary>
	/// Performs convolution with the specified convolution kernel on the given image.
	/// </summary>
	/// <param name="imgData"></param>
	/// <param name="bytesPerLine"></param>
	/// <param name="imgWidth"></param>
	/// <param name="imgHeight"></param>
	/// <param name="convolutionKernel"></param>
	/// <returns>Pointer to new image data if successful, nullptr otherwise.</returns>
	uchar* convolution(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, ConvolutionKernel* convolutionKernel);

	/// <summary>
	/// Performs filtration via explicit scheme for non-stationary heat equation.
	/// </summary>
	/// <param name="imgData"></param>
	/// <param name="bytesPerLine"></param>
	/// <param name="imgWidth"></param>
	/// <param name="imgHeight"></param>
	/// <param name="tau"></param>
	/// <param name="h"></param>
	/// <param name="timeSteps"></param>
	/// <returns></returns>
	uchar* filtrationExplicitHeatEq(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double tau, const double h, const int timeSteps);

	/// <summary>
	/// Performs filtration via implicit scheme for non-stationary heat equation.
	/// </summary>
	/// <param name="imgData"></param>
	/// <param name="bytesPerLine"></param>
	/// <param name="imgWidth"></param>
	/// <param name="imgHeight"></param>
	/// <param name="tau"></param>
	/// <param name="timeSteps"></param>
	/// <returns></returns>
	uchar* filtrationImplicitHeatEq(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double tau, const int timeSteps);

	/// <summary>
	/// Performs filtration via semi-implicit scheme for regularized Perona-Malik equation.
	/// </summary>
	/// <param name="imgData"></param>
	/// <param name="bytesPerLine"></param>
	/// <param name="imgWidth"></param>
	/// <param name="imgHeight"></param>
	/// <param name="sigma">-> step for linear diffusion.</param>
	/// <param name="tau">-> step for non-linear diffusion</param>
	/// <param name="K">-> </param>
	/// <param name="timeSteps"></param>
	/// <returns></returns>
	uchar* filtrationSemiImplicitPeronaMalik(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double sigma, const double tau, const double K, const int timeSteps);
	
	/// <summary>
	/// Performs filtration via semi-implicit Geodesic Mean Curvature Flow scheme.
	/// </summary>
	/// <param name="imgData">-> </param>
	/// <param name="bytesPerLine">-> </param>
	/// <param name="imgWidth">-> </param>
	/// <param name="imgHeight">-> </param>
	/// <param name="sigma">-> </param>
	/// <param name="tau">-> </param>
	/// <param name="K">-> if == 0 ==> Mean Curvature Flow scheme.</param>
	/// <param name="timeSteps"></param>
	/// <returns>...</returns>
	uchar* filtrationSemiImplicitGMCF_SOR(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double sigma, const double tau, const double K, const int timeSteps);

	/// <summary>
	/// Performs filtration via semi-implicit Geodesic Mean Curvature Flow scheme.
	/// </summary>
	/// <param name="imgData">-> </param>
	/// <param name="bytesPerLine">-> </param>
	/// <param name="imgWidth">-> </param>
	/// <param name="imgHeight">-> </param>
	/// <param name="sigma">-> </param>
	/// <param name="tau">-> </param>
	/// <param name="K">-> if == 0 ==> Mean Curvature Flow scheme.</param>
	/// <param name="timeSteps"></param>
	/// <returns>...</returns>
	uchar* filtrationSemiImplicitGMCF_BiCGStab(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double sigma, const double tau, const double K, const int timeSteps);


	/// <summary>
	/// 
	/// </summary>
	/// <param name="imgData"></param>
	/// <param name="bytesPerLine"></param>
	/// <param name="imgWidth"></param>
	/// <param name="imgHeight"></param>
	/// <returns></returns>
	double* signedDistanceFunction(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight);


	//################# Image Export functions #################//

	static bool exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, double* imgData, bool scaleData = true);
	static bool exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, uchar* imgData);
	static bool ExportToPPM(std::string fileName, int width, int height, int maxValue, float* r, float* g, float* b);

	static bool exportSgnFunction(std::string fileName, int width, int height, double* data);
};