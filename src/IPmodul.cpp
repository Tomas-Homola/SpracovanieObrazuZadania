#include "IPmodul.h"


IPmodul::IPmodul()
{

}


IPmodul::~IPmodul()
{
	free(m_imgLocalData);
}

double* IPmodul::getOriginalImgData()
{
	return nullptr;
}

bool IPmodul::mirrorPixels(uchar* originalImgData, const uint imgWidth, const uint imgHeight, const uint padding)
{
	int indexNew = 0, indexOld = 0;
	int temp = 0;

	// check if there is already some image data stored
	if (m_imgLocalData != nullptr)
	{
		free(m_imgLocalData); // if there is, delete old data
		m_imgLocalData = nullptr;
	}
	// allocate memory for new data
	int newWidth = imgWidth + 2 * padding;
	int newHeight = imgHeight + 2 * padding;
	int size = newWidth * newHeight; // compute new size
	m_padding = padding;
	m_imgLocalData = (double*)calloc(size, sizeof(double)); // allocate memory

	if (m_imgLocalData == nullptr) // check, if allocation was successful
		return false;

	// copy old image data
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			indexNew = (i + padding) * newWidth + (j + padding);
			indexOld = i * imgWidth + j;

			m_imgLocalData[indexNew] = static_cast<double>(originalImgData[indexOld]);
		}
	}

	// mirror over Upper and Lower edges
	temp = 1;
	for (int i = 0; i < padding; i++)
	{
		for (int j = padding; j < newWidth - padding; j++)
		{
			// upper egde
			indexOld = (i + padding) * newWidth + j;
			indexNew = (i + padding - temp) * newWidth + j;
			m_imgLocalData[indexNew] = m_imgLocalData[indexOld];

			// lower edge
			indexOld = (newHeight - i - padding - 1) * newWidth + j;
			indexNew = (newHeight - i - padding + temp - 1) * newWidth + j;
			m_imgLocalData[indexNew] = m_imgLocalData[indexOld];
		}
		temp += 2;
	}

	// mirror over Left and Right edges
	for (int i = 0; i < newHeight; i++)
	{
		temp = 1;
		for (int j = 0; j < padding; j++)
		{
			// left edge
			indexOld = i * newWidth + (j + padding);
			indexNew = i * newWidth + (j + padding - temp);
			m_imgLocalData[indexNew] = m_imgLocalData[indexOld];

			// right edge
			indexOld = i * newWidth + (newHeight - padding - 1 - j);
			indexNew = indexOld + temp;
			m_imgLocalData[indexNew] = m_imgLocalData[indexOld];

			temp += 2;
		}
	}

	// print data
	//for (int i = 0; i < newHeight; i++)
	//{
	//	for (int j = 0; j < newWidth; j++)
	//	{
	//		indexNew = i * newWidth + j;
	//
	//		printf("%.0f\t", m_imgLocalData[indexNew]);
	//	}
	//	printf("\n");
	//}

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

