#include "MyVoronoiShader.h"
#include <iostream>


GPixel cToP(const GColor& color) {
    float alpha = color.a * 255.0f;
    int red = GRoundToInt(color.r * alpha);
    int green = GRoundToInt(color.g * alpha);
    int blue = GRoundToInt(color.b * alpha);

    GPixel pix = GPixel_PackARGB(GRoundToInt(alpha), red, green, blue);
    return pix;
}

bool MyVoronoiShader::isOpaque() {
    return false;
}

bool MyVoronoiShader::setContext(const GMatrix& ctm) {
    std::optional<GMatrix> ctmInverse = ctm.invert();
    if (!ctmInverse) {return false;}
    vInverse = ctmInverse.value();
    return true;
}

void MyVoronoiShader::shadeRow(int x, int y, int count, GPixel row[]) {
    // for each incoming point:
    // map through inverse matrix
    // iterate through points array to find closet point
    // insert the color of the closest point
    for (int i = 0; i < count; i++) {
        GPoint pt = {x + i + 0.5f, y + 0.5f};
        vInverse.GMatrix::mapPoints(&pt, 1);

        float smallestDist = std::numeric_limits<float>::max();
        int smallestInd = -1;
        for (int j = 0; j < points.size(); j++) {
            float dist = (points[j] - pt).length();
            if (dist < smallestDist) {
                smallestDist = dist;
                smallestInd = j;
            }
        }

        GPixel pix = cToP(colors[smallestInd]);
        row[i] = pix;
    }
}

