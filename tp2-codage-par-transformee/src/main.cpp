// --- Including librairies
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

// --- Including CImg.h
#include "CImg.h"

// --- Defining the round "function"
#define ROUND(a) (((a) < 0) ? (int)((a)-0.5) : (int)((a) + 0.5))
#define SIZE 8

// --- Namespaces
using namespace std;
using namespace cimg_library;

/* 
 * In order to understand how the JPEG compression/decompression is working :
 * 1 . Realize the encoding, which is DCT + Quantization
 * 2 . Realize the decoding, which is inverse DCT + inverse Quantization
 */

// --- Helpful methods
CImg<> getQuantizationMatrix(float quality)
{
    CImg<> Q(8, 8);
    Q(0, 0) = 16;
    Q(0, 1) = 11;
    Q(0, 2) = 10;
    Q(0, 3) = 16;
    Q(0, 4) = 24;
    Q(0, 5) = 40;
    Q(0, 6) = 51;
    Q(0, 7) = 61;
    Q(1, 0) = 12;
    Q(1, 1) = 12;
    Q(1, 2) = 14;
    Q(1, 3) = 19;
    Q(1, 4) = 26;
    Q(1, 5) = 58;
    Q(1, 6) = 60;
    Q(1, 7) = 55;
    Q(2, 0) = 14;
    Q(2, 1) = 13;
    Q(2, 2) = 16;
    Q(2, 3) = 24;
    Q(2, 4) = 40;
    Q(2, 5) = 57;
    Q(2, 6) = 69;
    Q(2, 7) = 56;
    Q(3, 0) = 14;
    Q(3, 1) = 17;
    Q(3, 2) = 22;
    Q(3, 3) = 29;
    Q(3, 4) = 51;
    Q(3, 5) = 87;
    Q(3, 6) = 80;
    Q(3, 7) = 62;
    Q(4, 0) = 18;
    Q(4, 1) = 22;
    Q(4, 2) = 37;
    Q(4, 3) = 56;
    Q(4, 4) = 68;
    Q(4, 5) = 109;
    Q(4, 6) = 103;
    Q(4, 7) = 77;
    Q(5, 0) = 24;
    Q(5, 1) = 35;
    Q(5, 2) = 55;
    Q(5, 3) = 64;
    Q(5, 4) = 81;
    Q(5, 5) = 104;
    Q(5, 6) = 113;
    Q(5, 7) = 92;
    Q(6, 0) = 49;
    Q(6, 1) = 64;
    Q(6, 2) = 78;
    Q(6, 3) = 87;
    Q(6, 4) = 103;
    Q(6, 5) = 121;
    Q(6, 6) = 120;
    Q(6, 7) = 101;
    Q(7, 0) = 72;
    Q(7, 1) = 92;
    Q(7, 2) = 95;
    Q(7, 3) = 98;
    Q(7, 4) = 112;
    Q(7, 5) = 100;
    Q(7, 6) = 103;
    Q(7, 7) = 99;
    Q *= quality;

    return Q;
}

// --- Displaying an image
void display_image(CImg<double> image, const char *title)
{
    CImgDisplay tmp(image, title);
    while (!tmp.is_closed())
        tmp.wait();
}

// --- Apply DCT on one pixel
double DCTOnOnePx(int i, int j, CImg<double> sub_image)
{
    double pixel = 0;

    for (int a = 0; a < SIZE; ++a)
        for (int b = 0; b < SIZE; ++b)
            pixel += sub_image(a, b) * cos(((2 * a + 1) * i * M_PI) / (2 * SIZE)) * cos(((2 * b + 1) * j * M_PI) / (2 * SIZE));

    return pixel * (2.0 / SIZE) * (i == 0 ? (1.0 / sqrt(2.0)) : 1.0) * (j == 0 ? (1.0 / sqrt(2.0)) : 1.0);
}

// --- Apply Reverse DCT on one pixel
double ReverseDCTOnPx(int i, int j, CImg<double> sub_image)
{
    double pixel = 0;

    for (int a = 0; a < SIZE; ++a)
        for (int b = 0; b < SIZE; ++b)
            pixel += sub_image[b * SIZE + a] * cos(((2 * i + 1) * a * M_PI) / (2 * SIZE)) * cos(((2 * j + 1) * b * M_PI) / (2 * SIZE)) * (a == 0 ? (1.0 / sqrt(2.0)) : 1.0) * (b == 0 ? (1.0 / sqrt(2.0)) : 1.0);

    return pixel * (2.0 / SIZE);
}

// --- Apply DCT on one bloc of the image
void DCTOnOneBloc(int i, int j, CImg<double> &original_image, CImg<double> &compressed_image)
{
    CImg<double> sub_image = original_image.get_crop(i, j, i + (SIZE - 1), j + (SIZE - 1));

    for (int x = 0; x < SIZE; ++x)
        for (int y = 0; y < SIZE; ++y)
            compressed_image(x + i, y + j) = DCTOnOnePx(x, y, sub_image);
}

// --- Apply reverse DCT on one bloc of the image
void ReverseDCTOnOneBloc(int i, int j, CImg<double> &compressed_image, CImg<double> &original_image)
{
    CImg<double> sub_image = compressed_image.get_crop(i, j, i + (SIZE - 1), j + (SIZE - 1));

    for (int x = 0; x < SIZE; ++x)
        for (int y = 0; y < SIZE; ++y)
            original_image(x + i, y + j) = ReverseDCTOnPx(x, y, sub_image);
}

// --- Encoding method
CImg<double> JPEGEncoder(CImg<double> original_image, float quality)
{
    CImg<> Q = getQuantizationMatrix(quality);
    CImg<double> compressed_image(original_image.width(), original_image.height(), 1, 1, 0);
    compressed_image = original_image;

    // Level shifting the image to prepare for the DCT
    for (int i = 0; i < original_image.width(); ++i)
        for (int j = 0; j < original_image.height(); ++j)
            original_image(i, j) -= 128;

    // Applying the DCT by blocs of 8 by 8 pixels
    for (int i = 0; i <= original_image.width() - SIZE; i += SIZE)
        for (int j = 0; j <= original_image.height() - SIZE; j += SIZE)
            DCTOnOneBloc(i, j, original_image, compressed_image);

    // Quantization step
    for (int i = 0; i < original_image.width(); i++)
        for (int j = 0; j < original_image.height(); j++)
            compressed_image(i, j) = ROUND(compressed_image(i, j) / Q(i % SIZE, j % SIZE));

    return compressed_image;
}

// --- Decoding method
CImg<double> JPEGDecoder(CImg<double> compressed_image, float quality)
{
    CImg<> Q = getQuantizationMatrix(quality);
    CImg<double> original_image(compressed_image.width(), compressed_image.height(), 1, 1, 0);
    original_image = compressed_image;

    // Reverse Quantization
    for (int i = 0; i < compressed_image.width(); i++)
        for (int j = 0; j < compressed_image.height(); j++)
            original_image(i, j) = ROUND(original_image(i, j) * Q(i % SIZE, j % SIZE));

    // Applying the inverse DCT by blocs of 8 by 8 pixels
    for (int i = 0; i <= compressed_image.width() - SIZE; i += SIZE)
        for (int j = 0; j <= compressed_image.height() - SIZE; j += SIZE)
            ReverseDCTOnOneBloc(i, j, compressed_image, original_image);

    // Level shifting back the image
    for (int i = 0; i < compressed_image.width(); i++)
        for (int j = 0; j < compressed_image.height(); j++)
            original_image(i, j) += 128;

    return original_image;
}

// --- Distortion testings
void DistortionTests(double quality_start, double quality_end, double quality_step, CImg<double> original_image)
{
    CImg<double> compressed_image;

    vector<double> distortions;
    double distortion = 0;

    for (double quality = quality_start; quality <= quality_end; quality += quality_step)
    {
        compressed_image = JPEGEncoder(original_image, quality);
        distortion = 0;
        for (int i = 0; i < original_image.width(); ++i)
        {
            for (int j = 0; j < original_image.height(); ++j)
            {
                distortion += pow(original_image(i, j) - compressed_image(i, j), 2);
            }
        }
        distortions.push_back(distortion / (double)(original_image.width() * original_image.height()));
    }

    CImg<double> values(1, distortions.size(), 1, 1, 0);

    for (int i = 0; i < distortions.size(); ++i)
        values(0, i) = distortions[i];

    values.display_graph(NULL, 1, 1, "Quality factor", quality_start, quality_end, "Distortion rate");
}

void display_menu()
{

    cout << "1. Apply the DCT and display the result. You must provide a quality factor." << endl;
    cout << "2. Apply the DCT then the reverse DCT and display the decompressed image. You must provide a quality factor." << endl;
    cout << "3. Calculate the distortion rate according to the quality factor. You must provide a minimum quality factor, a maximum quality factor, and a pitch." << endl;
    cout << "You have chosen : ";
}

int main()
{
    int choice = 0;
    float quality = 1.;
    double quality_start, quality_end, quality_step;
    // Read the image "lena.bmp"
    CImg<double> my_image("lena.bmp");
    // Take the luminance information
    my_image.channel(0);

    display_menu();

    cin >> choice;
    switch (choice)
    {
    case 1:
        cout << "Quality factor : ";
        cin >> quality;
        display_image(JPEGEncoder(my_image, quality), "Compressed Image");
        break;
    case 2:
        cout << "Quality factor : ";
        cin >> quality;
        display_image(JPEGDecoder(JPEGEncoder(my_image, quality), quality), "Decompressed Image");
        break;
    case 3:
        cout << "Quality factor minimum : ";
        cin >> quality_start;
        cout << "Quality factor maximum : ";
        cin >> quality_end;
        cout << "Quality step : ";
        cin >> quality_step;
        DistortionTests(quality_start, quality_end, quality_step, my_image);
        break;
    default:
        break;
    }
}
