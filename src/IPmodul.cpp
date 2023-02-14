#include "IPmodul.h"


IPmodul::IPmodul()
{

}


IPmodul::~IPmodul()
{
	delete m_imgLocalData;
}

bool IPmodul::mirrorPixels(uchar* originalImgData, const uint imgWidth, const uint imgHeight, const uint N)
{
	int indexNew = 0, indexOld = 0;
	int temp = 0;

	// allocate memory for new data
	int newWidth = imgWidth + 2 * N;
	int newHeight = imgHeight + 2 * N;
	int size = newWidth * newHeight;
	m_imgLocalData = static_cast<double*>(calloc(size, sizeof(double)));
	if (m_imgLocalData == nullptr)
		return false;

	// copy old image data
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			indexNew = (i + N) * newWidth + (j + N);
			indexOld = i * imgWidth + j;

			m_imgLocalData[indexNew] = static_cast<double>(originalImgData[indexOld]);
		}
	}

	// mirror over Upper and Lower edges
	temp = 1;
	for (int i = 0; i < N; i++)
	{
		for (int j = N; j < newWidth - N; j++)
		{
			// upper egde
			indexOld = (i + N) * newWidth + j;
			indexNew = (i + N - temp) * newWidth + j;
			m_imgLocalData[indexNew] = m_imgLocalData[indexOld];

			// lower edge
			indexOld = (newHeight - i - N - 1) * newWidth + j;
			indexNew = (newHeight - i - N + temp - 1) * newWidth + j;
			m_imgLocalData[indexNew] = m_imgLocalData[indexOld];
		}
		temp += 2;
	}


	// print data
	for (int i = 0; i < newHeight; i++)
	{
		for (int j = 0; j < newWidth; j++)
		{
			indexNew = i * newWidth + j;

			printf("%.0f\t", m_imgLocalData[indexNew]);
		}
		printf("\n");
	}

	return true;
}

bool IPmodul::exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, double* imgData)
{
	printf("Exporting image to pgm...\n");
	FILE* fp = nullptr;
	fp = fopen((fileName + ".pgm").c_str(), "w+");
	if (fp == nullptr)
		return false;

	unsigned char scaledValue = 0;
	int dataSize = imgWidth * imgHeight;
	fprintf(fp, "P2\n%d %d\n%d\n", imgWidth, imgHeight, maxValue);
	for (size_t i = 0; i < dataSize; i++)
	{
		scaledValue = static_cast<unsigned char>(imgData[i] * maxValue + 0.5);
		fprintf(fp, "%d ", scaledValue);

		if ((i + 1) % 70 == 0)
			fprintf(fp, "\n");

		if ((i + 1) % (dataSize / 10) == 0)
			printf("%d%% done\n", 10 * ((int)i + 1) / (dataSize / 10));
	}
	fclose(fp);
	printf("Export done\n");

	return true;
}

