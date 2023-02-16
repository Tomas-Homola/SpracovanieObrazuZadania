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
	
	double* m_pImgLocalData = nullptr;

	uint m_imgWidth  = 0;
	uint m_imgHeight = 0;

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


	//################# Image Processing functions #################//

	/// <summary>
	/// Extend image by copying N pixels over its edges.
	/// </summary>
	/// <param name="data">-> input image data</param>
	/// <param name="imgWidth">-> input image width</param>
	/// <param name="imgHeight">-> input image height</param>
	/// <param name="padding">-> number of pixels to mirror</param>
	/// <returns>True if successful, false otherwise.</returns>
	bool pixelsMirror(uchar* originalImgData, const uint imgWidth, const uint imgHeight, const uint padding);

	/// <summary>
	/// Crop extended image based on the given padding.
	/// </summary>
	/// <param name="padding">-> ...</param>
	/// <returns>Pointer to memory where cropped image data are stored if successful, nullptr otherwise.</returns>
	uchar* pixelsUnmirror(uint padding);

	//################# Image Export functions #################//

	static bool exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, double* imgData);
	static bool exportToPGM(std::string fileName, uint imgWidth, uint imgHeight, int maxValue, uchar* imgData);

};