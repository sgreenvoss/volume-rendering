#include <iostream>
#include "DataStructs.h"
#include <cmath>

using namespace std;

TransferFunction
SetupTransferFunction(void)
{
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3 * 256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    for (i = 0; i < 256; i++)
        rv.opacities[i] = charOpacity[i] / 255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints * 3] = {
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0,
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76
    };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0; i < numControlPoints - 1; i++)
    {
        int start = controlPointPositions[i] * rv.numBins;
        int end = controlPointPositions[i + 1] * rv.numBins + 1;
        cerr << "Working on " << i << "/" << i + 1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins - 1;
        for (int j = start; j <= end; j++)
        {
            double proportion = (j / (rv.numBins - 1.0) - controlPointPositions[i]) / (controlPointPositions[i + 1] - controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0; k < 3; k++)
                rv.colors[3 * j + k] = proportion * (controlPointColors[3 * (i + 1) + k] - controlPointColors[3 * i + k])
                + controlPointColors[3 * i + k];
        }
    }

    return rv;
}

Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus.Set(0, 0, 0);
    rv.up.Set(0, -1, 0);
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position.Set(-8.25e+7, -3.45e+7, 3.35e+7);

    return rv;
}



/*
int main()
{
    TransferFunction tf = SetupTransferFunction();
    for (int i = 0 ; i < tf.numBins ; i++)
    {
        cerr << i << ": " << (int) tf.colors[3*i] << ", " << (int) tf.colors[3*i+1] << ", " << (int) tf.colors[3*i+2] << ", " << tf.opacities[i] << endl;
    }
}
 */