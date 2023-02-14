#pragma once

#include <string>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned char uchar;
typedef unsigned int uint;
// asi trieda "IPmodul", ktora bude mat funkcie na Image Processing
// na vstup do konstruktora dostane data, s ktorymi sa bude pracovat
// double* m_localData
// funkcia na zrkadlenie a odzrkadlenie -> mirror(int N), rozsiri obrazok o N pixelov kazdym smerom tak, ze prekopiruje hodnoty pixelov pri hranici obrazka

class IPmodul
{
private:
	//################# Variables #################//
	
	double* m_imgLocalData = nullptr;

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

	//################# Image Processing functions #################//

	/// <summary>
	/// Extend image by copying N pixels over its edges
	/// </summary>
	/// <param name="data">-> input image data</param>
	/// <param name="imgWidth">-> input image width</param>
	/// <param name="imgHeight">-> input image height</param>
	/// <param name="N">-> number of pixels to mirror</param>
	bool mirrorPixels(uchar* originalImgData, const uint imgWidth, const uint imgHeight, const uint N);


	//################# Image Export functions #################//

	static bool exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, double* imgData);

};