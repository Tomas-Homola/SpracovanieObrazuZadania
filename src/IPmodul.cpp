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

void IPmodul::updateEdges(const int padding, double* data)
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
			data[indexNew] = data[indexOld];

			// lower edge
			indexOld = (m_imgHeight - i - padding - 1) * m_imgWidth + j;
			indexNew = (m_imgHeight - i - padding + temp - 1) * m_imgWidth + j;
			data[indexNew] = data[indexOld];
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
			data[indexNew] = data[indexOld];

			// right edge
			indexOld = i * m_imgWidth + (m_imgWidth - padding - 1 - j);
			indexNew = i * m_imgWidth + (m_imgWidth - padding - 1 - j + temp);
			data[indexNew] = data[indexOld];

			temp += 2;
		}
	}
}

void IPmodul::uSigma_UpdateEdges(const int padding)
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

void IPmodul::uSigma_LinearDiffusion(double sigma)
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
	uSigma_UpdateEdges(padding);

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
			uSigma_UpdateEdges(padding);

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

void IPmodul::PM_computeMatrixCoefs(int padding, double tau, double K)
{
	double gradX = 0.0;
	double gradY = 0.0;
	double normSq = 0.0;

	double h = 1.0;

	int NW = 0, N = 0, NE = 0;
	int W  = 0, p = 0, E  = 0;
	int SW = 0, S = 0, SE = 0;

	// diffusion coefficients
	double g_N = 0.0; // diffusion ceofficient over North edge
	double g_S = 0.0; // diffusion ceofficient over South edge
	double g_E = 0.0; // diffusion ceofficient over East  edge
	double g_W = 0.0; // diffusion ceofficient over West  edge

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

			//----- EAST edge -----//
			gradX = (m_uSigma[E] - m_uSigma[p]) / h;
			gradY = (m_uSigma[N] + m_uSigma[NE] - m_uSigma[S] - m_uSigma[SE]) / (4.0 * h);

			normSq = gradX * gradX + gradY * gradY;
			g_E = diffCoefFunction1(K, normSq);

			//----- NORTH edge -----//
			gradX = (m_uSigma[W] + m_uSigma[NW] - m_uSigma[E] - m_uSigma[NE]) / (4.0 * h);
			gradY = (m_uSigma[N] - m_uSigma[p]) / h;

			normSq = gradX * gradX + gradY * gradY;
			g_N = diffCoefFunction1(K, normSq);

			//----- WEST edge -----//
			gradX = (m_uSigma[W] - m_uSigma[p]) / h;
			gradY = (m_uSigma[S] + m_uSigma[SW] - m_uSigma[N] - m_uSigma[NW]) / (4.0 * h);

			normSq = gradX * gradX + gradY * gradY;
			g_W = diffCoefFunction1(K, normSq);

			//----- SOUTH edge -----//
			gradX = (m_uSigma[E] + m_uSigma[SE] - m_uSigma[W] - m_uSigma[SW]) / (4.0 * h);
			gradY = (m_uSigma[S] - m_uSigma[p]) / h;

			normSq = gradX * gradX + gradY * gradY;
			g_S = diffCoefFunction1(K, normSq);


			m_matrixCoefs[p].Aii   = 1.0 + tau * (g_N + g_S + g_E + g_W);
			m_matrixCoefs[p].Aij_N = - tau * g_N;
			m_matrixCoefs[p].Aij_S = - tau * g_S;
			m_matrixCoefs[p].Aij_E = - tau * g_E;
			m_matrixCoefs[p].Aij_W = - tau * g_W;
		}
	}
}

void IPmodul::GMCF_computeMatrixCoefs(int padding, double epsilon, double tau, double K)
{
	double gradX = 0.0;
	double gradY = 0.0;
	double normSq = 0.0;
	double gradMean = 0.0;

	double h = 1.0;

	int NW = 0, N = 0, NE = 0;
	int W  = 0, p = 0, E  = 0;
	int SW = 0, S = 0, SE = 0;

	// system matrix coefficients
	double g_N = 0.0; // diffusion ceofficient over North edge
	double g_S = 0.0; // diffusion ceofficient over South edge
	double g_E = 0.0; // diffusion ceofficient over East  edge
	double g_W = 0.0; // diffusion ceofficient over West  edge

	double gradNorm_N = 0.0; // ||grad u||_epsilon over North edge
	double gradNorm_S = 0.0; // ||grad u||_epsilon over South edge
	double gradNorm_E = 0.0; // ||grad u||_epsilon over East  edge
	double gradNorm_W = 0.0; // ||grad u||_epsilon over West  edge

	double* dataOrig = m_pImgLocalData;

	double eps2 = epsilon * epsilon;

	// iterate over all internal img pixels of m_uSigma and m_pImgLocalData
	for (int I = padding; I < m_imgHeight - padding; I++)
	{
		for (int J = padding; J < m_imgWidth - padding; J++)
		{
			// compute correct indices
			NW = (I - 1) * m_imgWidth + (J - 1);
			N = (I - 1) * m_imgWidth + J;
			NE = (I - 1) * m_imgWidth + (J + 1);

			W = I * m_imgWidth + (J - 1);
			p = I * m_imgWidth + J;
			E = I * m_imgWidth + (J + 1);

			SW = (I + 1) * m_imgWidth + (J - 1);
			S = (I + 1) * m_imgWidth + J;
			SE = (I + 1) * m_imgWidth + (J + 1);


			//----- EAST edge -----//
			// uSigma
			gradX = (m_uSigma[E] - m_uSigma[p]) / h;
			gradY = (m_uSigma[N] + m_uSigma[NE] - m_uSigma[S] - m_uSigma[SE]) / (4.0 * h);

			normSq = gradX * gradX + gradY * gradY;
			g_E = diffCoefFunction1(K, normSq);

			// original data
			gradX = (dataOrig[E] - dataOrig[p]) / h;
			gradY = (dataOrig[N] + dataOrig[NE] - dataOrig[S] - dataOrig[SE]) / (4.0 * h);

			gradNorm_E = sqrt(eps2 + gradX * gradX + gradY * gradY);
			gradMean += gradNorm_E;


			//----- NORTH edge -----//
			// uSigma
			gradX = (m_uSigma[W] + m_uSigma[NW] - m_uSigma[E] - m_uSigma[NE]) / (4.0 * h);
			gradY = (m_uSigma[N] - m_uSigma[p]) / h;

			normSq = gradX * gradX + gradY * gradY;
			g_N = diffCoefFunction1(K, normSq);

			// original data
			gradX = (dataOrig[W] + dataOrig[NW] - dataOrig[E] - dataOrig[NE]) / (4.0 * h);
			gradY = (dataOrig[N] - dataOrig[p]) / h;

			gradNorm_N = sqrt(eps2 + gradX * gradX + gradY * gradY);
			gradMean += gradNorm_N;


			//----- WEST edge -----//
			// uSigma
			gradX = (m_uSigma[W] - m_uSigma[p]) / h;
			gradY = (m_uSigma[S] + m_uSigma[SW] - m_uSigma[N] - m_uSigma[NW]) / (4.0 * h);

			normSq = gradX * gradX + gradY * gradY;
			g_W = diffCoefFunction1(K, normSq);

			// original data
			gradX = (dataOrig[W] - dataOrig[p]) / h;
			gradY = (dataOrig[S] + dataOrig[SW] - dataOrig[N] - dataOrig[NW]) / (4.0 * h);

			gradNorm_W = sqrt(eps2 + gradX * gradX + gradY * gradY);
			gradMean += gradNorm_W;


			//----- SOUTH edge -----//
			// uSigma
			gradX = (m_uSigma[E] + m_uSigma[SE] - m_uSigma[W] - m_uSigma[SW]) / (4.0 * h);
			gradY = (m_uSigma[S] - m_uSigma[p]) / h;

			normSq = gradX * gradX + gradY * gradY;
			g_S = diffCoefFunction1(K, normSq);

			// original data
			gradX = (dataOrig[E] + dataOrig[SE] - dataOrig[W] - dataOrig[SW]) / (4.0 * h);
			gradY = (dataOrig[S] - dataOrig[p]) / h;

			gradNorm_S = sqrt(eps2 + gradX * gradX + gradY * gradY);
			gradMean += gradNorm_S;


			//----- Mean gradient -----//
			gradMean = gradMean / 4.0;
			
			m_matrixCoefs[p].Aii = 1.0 + tau * gradMean * (g_N / gradNorm_N + g_S / gradNorm_S + g_E / gradNorm_E + g_W / gradNorm_W);
			m_matrixCoefs[p].Aij_N = -tau * gradMean * (g_N / gradNorm_N);
			m_matrixCoefs[p].Aij_S = -tau * gradMean * (g_S / gradNorm_S);
			m_matrixCoefs[p].Aij_E = -tau * gradMean * (g_E / gradNorm_E);
			m_matrixCoefs[p].Aij_W = -tau * gradMean * (g_W / gradNorm_W);

			gradMean = 0.0;
		}
	}
}

void IPmodul::BiCGStab_compute_v_Ap(int imgWidth, int imgHeight, int padding, double* v, double* p)
{
	int I = 0;
	omp_set_dynamic(0);
#pragma omp parallel num_threads(m_threads)
	{
		int J = 0;

		int indexC = 0;
		int indexN = 0;
		int indexS = 0;
		int indexW = 0;
		int indexE = 0;

		int ic = 0;
		int in = 0;
		int is = 0;
		int iw = 0;
		int ie = 0;

		double Aii = 0.0;
		double Aij_N = 0.0;
		double Aij_S = 0.0;
		double Aij_E = 0.0;
		double Aij_W = 0.0;

		//i = 0;
		// iterate through inner part 
#pragma omp parallel for
		for (I = padding; I < m_imgHeight - padding; I++)
		{
			//j = 0;
			for (J = padding; J < m_imgWidth - padding; J++)
			{
				/*if (J < 4)
				{
					printf_s("%d z vlakna %d\n", I, omp_get_thread_num());
				}*/

				indexC = I * m_imgWidth + J;
				indexN = (I - 1) * m_imgWidth + J;
				indexS = (I + 1) * m_imgWidth + J;
				indexW = I * m_imgWidth + J - 1;
				indexE = I * m_imgWidth + J + 1;

				ic = (I - padding) * imgWidth + (J - padding);

				// compute system matrix coefficients
				Aii = m_matrixCoefs[indexC].Aii;
				Aij_N = m_matrixCoefs[indexC].Aij_N;
				Aij_S = m_matrixCoefs[indexC].Aij_S;
				Aij_E = m_matrixCoefs[indexC].Aij_E;
				Aij_W = m_matrixCoefs[indexC].Aij_W;

				v[ic] = Aii * p[indexC] + Aij_N * p[indexN] + Aij_S * p[indexS] + Aij_E * p[indexE] + Aij_W * p[indexW];

				//j++;
			}
			//i++;
		}
	}
}

void IPmodul::BiCGStab_compute_t_As(int imgWidth, int imgHeight, int padding, double* t, double* s)
{
	//int i = 0, j = 0;
	int I = 0;
	omp_set_dynamic(0);
#pragma omp parallel num_threads(m_threads)
	{
		int J = 0;

		int indexC = 0;
		int indexN = 0;
		int indexS = 0;
		int indexW = 0;
		int indexE = 0;

		int ic = 0;
		int in = 0;
		int is = 0;
		int iw = 0;
		int ie = 0;

		double Aii = 0.0;
		double Aij_N = 0.0;
		double Aij_S = 0.0;
		double Aij_E = 0.0;
		double Aij_W = 0.0;

		//i = 0;
		// iterate through inner part 
#pragma omp parallel for
		for (I = padding; I < m_imgHeight - padding; I++)
		{
			//j = 0;
			for (J = padding; J < m_imgWidth - padding; J++)
			{
				/*if (J < 4)
				{
					printf_s("%d z vlakna %d\n", I, omp_get_thread_num());
				}*/

				indexC = I * m_imgWidth + J;
				indexN = (I - 1) * m_imgWidth + J;
				indexS = (I + 1) * m_imgWidth + J;
				indexW = I * m_imgWidth + J - 1;
				indexE = I * m_imgWidth + J + 1;

				ic = (I - padding) * imgWidth + (J - padding);

				// compute system matrix coefficients
				Aii   = m_matrixCoefs[indexC].Aii;
				Aij_N = m_matrixCoefs[indexC].Aij_N;
				Aij_S = m_matrixCoefs[indexC].Aij_S;
				Aij_E = m_matrixCoefs[indexC].Aij_E;
				Aij_W = m_matrixCoefs[indexC].Aij_W;

				t[ic] = Aii * s[indexC] + Aij_N * s[indexN] + Aij_S * s[indexS] + Aij_E * s[indexE] + Aij_W * s[indexW];

				//j++;
			}
			//i++;
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

	double* b = new double[size] {}; // right side of the system
	double* phi = m_pImgLocalData; // = m_pImgLocalData
	m_uSigma = new double[(size_t)m_imgWidth * m_imgHeight] {}; // allocate space for uSigma, same as m_pImgLocalData

	// resize vector for gradient norms squared
	m_matrixCoefs.resize((size_t)m_imgWidth * m_imgHeight, MatrixCoefs());

	// indices
	int indexC = -1, indexN = -1, indexS = -1, indexW = -1, indexE = -1, indexNew = -1;

	// auxilliary variables
	double newValue = 0.0;
	double sum = 0.0;
	double temp = 0.0;
	uchar origValue = 0;

	// system matrix coefficients
	double Aii   = 0.0; // = 1 + tau * (g_N + g_S + g_E + g_W)
	double Aij_N = 0.0; // = - tau * g_N
	double Aij_S = 0.0; // = - tau * g_S
	double Aij_E = 0.0; // = - tau * g_E
	double Aij_W = 0.0; // = - tau * g_W

	// compute mean value of the original image
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
	const double TOL = 1.0E-8;
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
		uSigma_LinearDiffusion(sigma);

		// compute system matrix coefficients for all pixels
		//uSigma_ComputeGradientsNormSquared(padding);
		PM_computeMatrixCoefs(padding, tau, K);

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

					// get system matrix coefficients
					Aii   = m_matrixCoefs[indexC].Aii;
					Aij_N = m_matrixCoefs[indexC].Aij_N;
					Aij_S = m_matrixCoefs[indexC].Aij_S;
					Aij_E = m_matrixCoefs[indexC].Aij_E;
					Aij_W = m_matrixCoefs[indexC].Aij_W;

					sigmaSOR = Aij_N * phi[indexN] + Aij_S * phi[indexS] + Aij_E * phi[indexE] + Aij_W * phi[indexW];

					phi[indexC] = (1.0 - omega) * phi[indexC] + (omega / Aii) * (b[i * imgWidth + j] - sigmaSOR);

					j++;
				}
				i++;
			}

			updateEdges(padding);

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

					// get system matrix coefficients
					Aii = m_matrixCoefs[indexC].Aii;
					Aij_N = m_matrixCoefs[indexC].Aij_N;
					Aij_S = m_matrixCoefs[indexC].Aij_S;
					Aij_E = m_matrixCoefs[indexC].Aij_E;
					Aij_W = m_matrixCoefs[indexC].Aij_W;

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

		// print final SOR iter and its residuum
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
	m_matrixCoefs.clear(); // clear vector for gradient norms squared

	// return filtered data
	resultData = pixelsUnmirror(padding);
	return resultData;
}

uchar* IPmodul::filtrationSemiImplicitGMCF_SOR(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double sigma, const double tau, const double K, const int timeSteps)
{
	// mirror pixels
	int padding = 1;
	pixelsMirror(imgData, bytesPerLine, imgWidth, imgHeight, padding);

	size_t size = (size_t)imgWidth * imgHeight; // size of the original image
	size_t size2 = (size_t)m_imgWidth * m_imgHeight; // size of the original image
	uchar* resultData = nullptr; // return data in uchar after filtration

	double* b = new double[size] {}; // right side of the system
	double* phi = m_pImgLocalData; // = m_pImgLocalData
	m_uSigma = new double[size2] {}; // allocate space for uSigma, same as m_pImgLocalData

	// resize vector for system matrix coefficients
	//m_GMCF_GradientNorms.resize((size_t)m_imgWidth * m_imgHeight, AllGradientsNormSquared());
	m_matrixCoefs.resize(size2, MatrixCoefs());

	// indices
	int indexC = -1, indexN = -1, indexS = -1, indexW = -1, indexE = -1, indexNew = -1;

	// auxilliary variables
	//double newValue = 0.0;
	double temp = 0.0;
	uchar origValue = 0;

	// system matrix coefficients
	double Aii   = 0.0; // = 1 + tau * (g_N + g_S + g_E + g_W)
	double Aij_N = 0.0; // = - tau * g_N
	double Aij_S = 0.0; // = - tau * g_S
	double Aij_E = 0.0; // = - tau * g_E
	double Aij_W = 0.0; // = - tau * g_W

	// set right hand for the first time iteration
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			origValue = imgData[i * bytesPerLine + j];
			temp = static_cast<double>(origValue);
			
			// set right side to u0 = original image data
			b[i * imgWidth + j] = temp;
		}
	}

	/*if (printMsg)
		printf("Original image mean value: %.12lf\n", sum / size);*/


	// SOR variables
	double epsilon = m_MCF_epsilon; // regularization parameter for this model
	double omega = 1.25;
	const int MAX_ITER = 9999;
	const double TOL = 1.0E-3;
	int iter = 0;
	double rez = 0.0;
	double sigmaSOR = 0.0; // SOR algorithm sigma variable

	int i = 0, j = 0;

	// iterate through time steps
	for (int t = 1; t <= timeSteps; t++)
	{
		i = 0; j = 0;

		// perform one time step of linear diffusion with step sigma
		uSigma_LinearDiffusion(sigma);

		// compute gradient norms
		//GMCF_computeAllGradients(padding, epsilon);
		GMCF_computeMatrixCoefs(padding, epsilon, tau, K);

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

					// get system matrix coefficients
					Aii   = m_matrixCoefs[indexC].Aii;
					Aij_N = m_matrixCoefs[indexC].Aij_N;
					Aij_S = m_matrixCoefs[indexC].Aij_S;
					Aij_E = m_matrixCoefs[indexC].Aij_E;
					Aij_W = m_matrixCoefs[indexC].Aij_W;

					sigmaSOR = Aij_N * phi[indexN] + Aij_S * phi[indexS] + Aij_E * phi[indexE] + Aij_W * phi[indexW];

					phi[indexC] = (1.0 - omega) * phi[indexC] + (omega / Aii) * (b[i * imgWidth + j] - sigmaSOR);

					j++;
				}
				i++;
			}

			updateEdges(padding);

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

					// get system matrix coefficients
					Aii = m_matrixCoefs[indexC].Aii;
					Aij_N = m_matrixCoefs[indexC].Aij_N;
					Aij_S = m_matrixCoefs[indexC].Aij_S;
					Aij_E = m_matrixCoefs[indexC].Aij_E;
					Aij_W = m_matrixCoefs[indexC].Aij_W;

					// compute Ax - b
					temp = (Aii * phi[indexC] + Aij_N * phi[indexN] + Aij_S * phi[indexS] + Aij_E * phi[indexE] + Aij_W * phi[indexW]) - b[i * imgWidth + j];
					rez += temp * temp;

					j++;
				}
				i++;
			}
			rez = sqrt(rez);

			if (iter % 1000 == 0)
				printf("time %4d -> SOR iter %d: rez = %.8lf\n", t, iter, rez);
			
			if (rez < TOL)
				break;

		} while (iter <= MAX_ITER); // END SOR

		// print final SOR iter and its residuum
		printf("time %d -> SOR stop iter %d: rez = %.8lf\n", t, iter, rez);

		// update b vector
		i = 0;
		for (int I = padding; I < m_imgHeight - padding; I++)
		{
			j = 0;
			for (int J = padding; J < m_imgWidth - padding; J++)
			{
				indexC = I * m_imgWidth + J;

				// update b vector as the solution from the last time iteration
				b[i * imgWidth + j] = phi[indexC];
				j++;
			}
			i++;
		}

		if (printMsg)
			printf("Time step %d filtered image GCMF done\n", t);

	}
	if (printMsg)
		printf("\n\n");

	delete[] b; // clear space for vector b
	delete[] m_uSigma; // clear space for m_uSigma
	m_matrixCoefs.clear(); // clear vector for gradient norms squared

	// return filtered data
	resultData = pixelsUnmirror(padding);
	return resultData;
}

uchar* IPmodul::filtrationSemiImplicitGMCF_BiCGStab(uchar* imgData, const int bytesPerLine, const int imgWidth, const int imgHeight, const double sigma, const double tau, const double K, const int timeSteps)
{
	// mirror pixels
	int padding = 1;
	pixelsMirror(imgData, bytesPerLine, imgWidth, imgHeight, padding);

	size_t size = (size_t)imgWidth * imgHeight; // size of the original image
	size_t size2 = (size_t)m_imgWidth * m_imgHeight; // size of the mirrored image
	uchar* resultData = nullptr; // return data in uchar after filtration

	double* b = new double[size] {}; // right side of the system
	double* xSol = m_pImgLocalData; // = m_pImgLocalData

	m_uSigma = new double[size2] {}; // allocate space for uSigma, same as m_pImgLocalData

	// resize vector for gradient norms squared
	m_matrixCoefs.resize(size2, MatrixCoefs());

	// indices
	int indexC = -1, indexN = -1, indexS = -1, indexW = -1, indexE = -1, indexNew = -1;
	int ic = -1, in = -1, is = -1, iw = -1, ie = -1;

	// auxilliary variables
	//double newValue = 0.0;
	double temp = 0.0, temp2 = 0.0;
	uchar origValue = 0;

	// system matrix coefficients
	double Aii   = 0.0; // = 1 + tau * (g_N + g_S + g_E + g_W)
	double Aij_N = 0.0; // = - tau * g_N
	double Aij_S = 0.0; // = - tau * g_S
	double Aij_E = 0.0; // = - tau * g_E
	double Aij_W = 0.0; // = - tau * g_W

	int i = 0, j = 0;

	// set right hand for the first time iteration
	for (i = 0; i < imgHeight; i++)
	{
		for (j = 0; j < imgWidth; j++)
		{
			origValue = imgData[i * bytesPerLine + j];
			temp = static_cast<double>(origValue);

			// set right side to u0 = original image data
			b[i * imgWidth + j] = temp;
		}
	}

	// BiCGStab variables
	double* r_hat = new double[size] {};  // vektor tilde{r} = b - A.x^0; [size]
	double* r     = new double[size] {};  // vektor pre rezidua [size]
	double* p     = new double[size2] {}; // mirrored [size2]
	double* v     = new double[size] {};  //  [size]
	double* s     = new double[size2] {}; // mirrored [size2]
	double* t     = new double[size] {};  //  [size]

	double beta = 0.0;
	double rhoNew = 1.0;
	double rhoOld = 0.0;
	double alpha = 1.0;
	double omega = 1.0;

	double tempDot = 0.0;
	double tempDot2 = 0.0;
	double sNorm = 0.0;

	double epsilon = m_MCF_epsilon; // regularization parameter for this model
	int MAX_ITER = 9999;
	double TOL = 1.0E-4;
	int iter = 1;

	double rezNorm = 0.0;

	// iterate through time steps
	for (int time = 1; time <= timeSteps; time++)
	{
		i = 0; j = 0;

		// perform one time step of linear diffusion with step sigma
		uSigma_LinearDiffusion(sigma);

		// compute gradient norms
		GMCF_computeMatrixCoefs(padding, epsilon, tau, K);

		// set "r" and "r_hat" vectors to right hand side vector "b"
		rezNorm = 0.0;
		for (int k = 0; k < size; k++)
		{
			r[k] = b[k];
			r_hat[k] = b[k];

			// compute initial ||r||
			rezNorm += r[k] * r[k];
		}

		printf("||r0||: %.10lf\n", sqrt(rezNorm));
		rezNorm = 0.0;


		// tu sa asi bude diat BiCGStab
		iter = 1;
		rezNorm = 0.0;
		
		beta = 0.0;
		rhoNew = 1.0;
		rhoOld = 0.0;
		alpha = 1.0;
		omega = 1.0;
		
		// always start from zeros
		for (int k = 0; k < size2; k++)
		{
			xSol[k] = 0.0;
		}

		do
		{
			rhoOld = rhoNew; // save previous rho_{i-2}
			rhoNew = 0.0; // compute new rho_{i-1}
			for (int k = 0; k < size; k++)
			{
				temp = r_hat[k];
				temp2 = r[k];
				rhoNew += r_hat[k] * r[k];
			}
				

			// check if rhoNew > 0
			if (abs(rhoNew) < DBL_EPSILON)
			{
				delete[] b; // clear space for vector b
				delete[] m_uSigma; // clear space for m_uSigma
				// clear BiCGStab variables
				delete[] r_hat;
				delete[] r;
				delete[] p;
				delete[] v;
				delete[] s;
				delete[] t;

				m_matrixCoefs.clear();

				printf("Time step %d: BiCGStab failed -> rhoNew < 0\n", time);
				return nullptr;
			}

			beta = (rhoNew / rhoOld) * (alpha / omega);
			
			// update vector p^(i)
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;

					ic = i * imgWidth + j;

					p[indexC] = r[ic] + beta * (p[indexC] - omega * v[ic]);

					j++;
				}
				i++;
			}
			
			updateEdges(padding, p); // update edges of p

			// compute vector v = A.p
			BiCGStab_compute_v_Ap(imgWidth, imgHeight, padding, v, p);

			// compute alpha
			tempDot = 0.0;
			for (int k = 0; k < size; k++)
				tempDot += r_hat[k] * v[k];

			alpha = rhoNew / tempDot;

			// compute vektor s
			sNorm = 0.0;
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;

					ic = i * imgWidth + j;

					s[indexC] = r[ic] - alpha * v[ic];
					sNorm += s[indexC] * s[indexC];
					
					j++;
				}
				i++;
			}

			updateEdges(padding, s);

			sNorm = sqrt(sNorm);
			if (sNorm < TOL) // check if ||s|| is small enough
			{
				// update solution x
				i = 0;
				for (int I = padding; I < m_imgHeight - padding; I++)
				{
					j = 0;
					for (int J = padding; J < m_imgWidth - padding; J++)
					{
						indexC = I * m_imgWidth + J;
						ic = i * imgWidth + j;

						xSol[indexC] = xSol[indexC] + alpha * p[indexC];

						j++;
					}
					i++;
				}

				//printf("BCGS stop: ||s|| is small enough, iter: %d\n", iter);
				break;
			}

			// compute vector t = A.s
			BiCGStab_compute_t_As(imgWidth, imgHeight, padding, t, s);
			
			// compute omega
			tempDot = 0.0; tempDot2 = 0.0;
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;
					ic = i * imgWidth + j;

					tempDot += t[ic] * s[indexC];
					tempDot2 += t[ic] * t[ic];

					j++;
				}
				i++;
			}
			
			omega = tempDot / tempDot2;
			
			rezNorm = 0.0;
			// update solution x and compute rez norm
			i = 0;
			for (int I = padding; I < m_imgHeight - padding; I++)
			{
				j = 0;
				for (int J = padding; J < m_imgWidth - padding; J++)
				{
					indexC = I * m_imgWidth + J;
					ic = i * imgWidth + j;

					xSol[indexC] = xSol[indexC] + alpha * p[indexC] + omega * s[indexC]; // update solution x

					r[ic] = s[indexC] - omega * t[ic]; // compute new residuum vector
					rezNorm += r[ic] * r[ic]; // compute residuum norm

					j++;
				}
				i++;
			}

			rezNorm = sqrt(rezNorm);
			if (iter % 500 == 0)
			{
				printf("iter: %4d ||r||: %.10lf\n", iter, rezNorm);
			}

			if (rezNorm < TOL)
			{
				//printf("BCGS stop iter: ||r|| is small enough\n");
				break;
			}
			
			iter++;

		} while ((iter < MAX_ITER) && (rezNorm > TOL)); // END BiCGStab

		printf("time %d -> BiCGStab stop iter %d: rez = %.8lf\n", time, iter, rezNorm);

		for (int k = 0; k < size; k++)
		{
			r_hat[k] = 0.0;
			r[k] = 0.0;
			v[k] = 0.0;
			t[k] = 0.0;
		}

		for (int k = 0; k < size2; k++)
		{
			p[k] = 0.0;
			s[k] = 0.0;
		}

		// print final SOR iter and its residuum
		//printf("time %d -> SOR stop iter %d: rez = %.8lf\n", t, iter, rez);

		// update b vector; TODO: netreba pocitat, ak je uz posledna casova iteracia
		//double max = DBL_MIN;
		//double min = DBL_MAX;
		i = 0;
		for (int I = padding; I < m_imgHeight - padding; I++)
		{
			j = 0;
			for (int J = padding; J < m_imgWidth - padding; J++)
			{
				indexC = I * m_imgWidth + J;
				ic = i * imgWidth + j;

				// update b vector as the solution from the last time iteration
				b[ic] = xSol[indexC];
				//if (xSol[indexC] > max)
				//	max = xSol[indexC];
				//if (xSol[indexC] < min)
				//	min = xSol[indexC];

				j++;
			}
			i++;
		}

		//printf("max: %.3lf\t\tmin: %.3lf\n", max, min);
		// update edges
		updateEdges(padding);

		if (printMsg)
			printf("Time step %d filtered image GCMF done\n", time);

	}
	if (printMsg)
		printf("\n");

	delete[] b; // clear space for vector b
	delete[] m_uSigma; // clear space for m_uSigma
	// clear BiCGStab variables
	delete[] r_hat;
	delete[] r;
	delete[] p;
	delete[] v;
	delete[] s;
	delete[] t;

	m_matrixCoefs.clear(); // clear vector for gradient norms squared

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