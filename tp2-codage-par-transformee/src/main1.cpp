#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <vector>

#include "CImg.h"

#define PI 3.14159265
#define BLOCK_SIZE 8

#define ROUND( a ) ( ( (a) < 0 ) ? (int) ( (a) - 0.5 ) : (int) ( (a) + 0.5 ) )

using namespace std;
using namespace cimg_library;

CImg<> getQ(float quality) {
    CImg<> Q(8,8);
    
    Q(0,0)=16;   Q(0,1)=11;   Q(0,2)=10;   Q(0,3)=16;   Q(0,4)=24;   Q(0,5)=40;   Q(0,6)=51;   Q(0,7)=61;
    Q(1,0)=12;   Q(1,1)=12;   Q(1,2)=14;   Q(1,3)=19;   Q(1,4)=26;   Q(1,5)=58;   Q(1,6)=60;   Q(1,7)=55;
    Q(2,0)=14;   Q(2,1)=13;   Q(2,2)=16;   Q(2,3)=24;   Q(2,4)=40;   Q(2,5)=57;   Q(2,6)=69;   Q(2,7)=56;
    Q(3,0)=14;   Q(3,1)=17;   Q(3,2)=22;   Q(3,3)=29;   Q(3,4)=51;   Q(3,5)=87;   Q(3,6)=80;   Q(3,7)=62;
    Q(4,0)=18;   Q(4,1)=22;   Q(4,2)=37;   Q(4,3)=56;   Q(4,4)=68;   Q(4,5)=109;  Q(4,6)=103;  Q(4,7)=77;
    Q(5,0)=24;   Q(5,1)=35;   Q(5,2)=55;   Q(5,3)=64;   Q(5,4)=81;   Q(5,5)=104;  Q(5,6)=113;  Q(5,7)=92;
    Q(6,0)=49;   Q(6,1)=64;   Q(6,2)=78;   Q(6,3)=87;   Q(6,4)=103;  Q(6,5)=121;  Q(6,6)=120;  Q(6,7)=101;
    Q(7,0)=72;   Q(7,1)=92;   Q(7,2)=95;   Q(7,3)=98;   Q(7,4)=112;  Q(7,5)=100;  Q(7,6)=103;  Q(7,7)=99;

    return Q * quality;
}

// Encode 1 pixel
double PixelEncoder(CImg<double> sub_image, int i, int j) {
    double pixel = 0;
    
    for (int x1 = 0; x1 < BLOCK_SIZE; ++x1)
        for (int y1 = 0; y1 < BLOCK_SIZE; ++y1)
            pixel += sub_image(x1, y1) * cos(((2 * x1 + 1) * i * PI) / (2 * BLOCK_SIZE)) * cos(((2 * y1 + 1) * j * PI) / (2 * BLOCK_SIZE));
    
    return pixel * (2.0 / BLOCK_SIZE) * (i == 0 ? (1.0 / sqrt(2.0)) : 1.0) * (j == 0 ? (1.0 / sqrt(2.0)) : 1.0);
}

// Parcours chaque bloc de l'image et l'encode
void BlocEncoder(CImg<double> & originalImage, CImg<double> & comp, int i0, int j0) {
    CImg<double> sub_image = originalImage.get_crop(i0, j0, i0 + (BLOCK_SIZE - 1), j0 + (BLOCK_SIZE - 1));

    for (int i = 0; i < BLOCK_SIZE; ++i)
        for (int j = 0; j < BLOCK_SIZE; ++j)
            comp(i + i0, j + j0) = PixelEncoder(sub_image, i, j);
}

// Parcours l'image par bloc de pixels et les encode
CImg<double> JPEGEncoder(CImg<double> image, float quality) {
	CImg<double> comp(image.width(), image.height(), 1, 1, 0);
	comp = image;
    CImg<> Q = getQ(quality);

    // Centre l'image
    for(int i = 0; i < image.width(); ++i)
        for(int j = 0; j < image.height(); ++j)
            image(i,j) -= 128;

	// Parcours de l'image original par bloc de 8 pixels
	for (int i = 0; i <= image.width() - BLOCK_SIZE; i += BLOCK_SIZE)
		for (int j = 0; j <= image.height() - BLOCK_SIZE; j += BLOCK_SIZE)
            BlocEncoder(image, comp, i, j);

    // Etape de quantification
    for(int i = 0; i < image.width(); i++)
        for(int j = 0; j < image.height(); j++)
            comp(i, j) = ROUND(comp(i, j) / Q(i % BLOCK_SIZE, j % BLOCK_SIZE));

	return comp;
}

// Decode 1 pixel
double PixelDecoder(CImg<double> sub_image, int i, int j) {
    double pixel = 0;
    
    for (int x1 = 0; x1 < BLOCK_SIZE; ++x1)
        for (int y1 = 0; y1 < BLOCK_SIZE; ++y1)
            pixel += sub_image[y1 * BLOCK_SIZE + x1] * cos(((2 * i + 1) * x1 * PI) / (2 * BLOCK_SIZE)) * cos(((2 * j + 1) * y1 * PI) / (2 * BLOCK_SIZE)) * 
			        (x1 == 0 ? (1.0 / sqrt(2.0)) : 1.0) * (y1 == 0 ? (1.0 / sqrt(2.0)) : 1.0);
    
    return pixel * (2.0 / BLOCK_SIZE);
}

// Parcours chaque bloc de l'image et le decode
void BlocDecoder(CImg<double> & originalImage, CImg<double> & comp, int i0, int j0) {
    CImg<double> sub_image = originalImage.get_crop(i0, j0, i0 + (BLOCK_SIZE - 1), j0 + (BLOCK_SIZE - 1));
    
    for (int i = 0; i < BLOCK_SIZE; ++i)
        for (int j = 0; j < BLOCK_SIZE; ++j)
            comp(i + i0, j + j0) = PixelDecoder(sub_image, i, j);
}

// Parcours l'image par bloc de pixels et les encode
CImg<double> JPEGDecoder(CImg<double> image, float quality) {
	CImg<double> comp(image.width(), image.height(), 1, 1, 0);
	comp = image;
    CImg<> Q = getQ(quality);

    // Etape de quantification inverse
    for(int n = 0; n < image.width(); n++)
        for(int m = 0; m < image.height(); m++)
            comp(n,m) = ROUND(comp(n,m) * Q(n % BLOCK_SIZE, m % BLOCK_SIZE));

	// Parcours de l'image original par bloc de 8 pixels
	for (int i0 = 0; i0 <= image.width() - BLOCK_SIZE; i0 += BLOCK_SIZE)
		for (int j0 = 0; j0 <= image.height() - BLOCK_SIZE; j0 += BLOCK_SIZE)
            BlocDecoder(image, comp, i0, j0);

    // Repositionne de l'image
    for(int n = 0; n < image.width(); n++)
        for(int m = 0; m < image.height(); m++)
            comp(n,m) += 128;

	return comp;
}

// Affiche un graph de l'Ã©volution de la distorsion par rapport au facteur qualitÃ©
void display_graph(vector<double> distorsions, double start, double end) {
    CImg<double> values(1, distorsions.size(), 1, 1, 0);

    for (int i = 0; i < distorsions.size(); ++i)
        values(0, i) = distorsions[i];

    values.display_graph(NULL, 1, 1, "Facteur de qualitÃ©", start, end, "Taux de distorsion");
}

// Calcule la distorsion d'une image compressÃ©e
double JPEGDistorsion(CImg<double> image, CImg<double> compressedImage) {
    double distorsion = 0;

	for (int i = 0; i < image.width(); ++i)
		for (int j = 0; j < image.height(); ++j)
			distorsion += pow(image(i, j) - 128 - compressedImage(i, j), 2);

	return distorsion / (double)(image.width() * image.height()); 
}

// Affiche la distorsion de chaque images compressÃ©es
void JPEGDistorsions(CImg<double> image, double start, double end, double step) {
    vector<double> distorsions;

    for (double quality=start; quality <= end; quality += step)
        distorsions.push_back(JPEGDistorsion(image, JPEGEncoder(image, quality)));

    display_graph(distorsions, start, end);
}

// Affiche l'image passÃ©e en paramÃ¨tre
void display_image(CImg<double> image, const char * title) {
	CImgDisplay tmp(image, title);
    while (!tmp.is_closed())
        tmp.wait();
}

// Lance l'outil de compression et de dÃ©compression de l'image
int main(int argc, char** const argv) {
	// Recupere l'image "lena.bmp"
	CImg<unsigned char> my_image("lena.bmp");

	// Prend l'information de luminance 
	my_image.channel(0);

    int choice = 0;
    double quality, start, end, step;
	CImg<double> compressed, decompressed;

    cout << "1. Appliquer la DCT et afficher le resultat. Vous devez fournir un facteur de qualite." << endl;
    cout << "2. Appliquer la DCT puis la DCT inverse et afficher l'image decompressee. Vous devez fournir un facteur de qualite." << endl;
    cout << "3. Calculer le taux de distorsion en fonction du facteur de qualite. Vous devez fournir un facteur de qualite minimum, un facteur de qualite maximum, et un pas." << endl;
    cout << "Autre. Quitter" << endl;
    cout << "Choix : ";
    cin >> choice;

    if (choice == 1) {
        cout << "QualitÃ© : ";
        cin >> quality;
        display_image(JPEGEncoder(my_image, quality), "Compressed Image");
    }
    else if (choice == 2) {
        cout << "QualitÃ© : ";
        cin >> quality;
        display_image(my_image, "Initial Image");
        display_image(JPEGDecoder(JPEGEncoder(my_image, quality), quality), "Decompressed Image");
    }
    else if (choice == 3) {
        cout << "Borne infÃ©rieure : ";
        cin >> start;
        cout << "Borne supÃ©rieure : ";
        cin >> end;
        cout << "Pas : ";
        cin >> step;
        JPEGDistorsions(my_image, start, end, step);
    }
}