#include "IPmodul.h"


IPmodul::IPmodul()
{

}


IPmodul::~IPmodul()
{
	free(m_pImgLocalData);
}

bool IPmodul::pixelsMirror(uchar* originalImgData, const uint imgWidth, const uint imgHeight, const uint padding)
{
	int indexNew = 0, indexOld = 0;
	int temp = 0;

	// check if there is already some image data stored
	if (m_pImgLocalData != nullptr)
	{
		free(m_pImgLocalData); // if there is, delete old data
		m_pImgLocalData = nullptr;
	}
	
	// compute new size
	//int newWidth = imgWidth + 2 * padding;
	m_imgWidth = imgWidth + 2 * padding;
	//int newHeight = imgHeight + 2 * padding;
	m_imgHeight = imgHeight + 2 * padding;
	int size = m_imgWidth * m_imgHeight;
	
	m_pImgLocalData = (double*)calloc(size, sizeof(double)); // allocate memory

	if (m_pImgLocalData == nullptr) // check, if allocation was successful
		return false;

	// copy old image data
	for (int i = 0; i < imgHeight; i++)
	{
		for (int j = 0; j < imgWidth; j++)
		{
			indexNew = (i + padding) * m_imgWidth + (j + padding);
			indexOld = i * imgWidth + j;

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
			indexOld = i * m_imgWidth + (m_imgHeight - padding - 1 - j);
			indexNew = indexOld + temp; // mozno funguje, mozno nie, neviem :(
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

uchar* IPmodul::pixelsUnmirror(uint padding)
{
	if (m_pImgLocalData == nullptr) // check, if there is some image to be cropped
		return nullptr;

	// calculate new size
	uint newWidth  = m_imgWidth - 2 * padding;
	uint newHeight = m_imgHeight - 2 * padding;
	uint size = newWidth * newHeight;

	uchar* pImgData = (uchar*)calloc(size, sizeof(uchar)); // allocate new memory
	
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

