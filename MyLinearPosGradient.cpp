#include "MyLinearPosGradient.h"
#include <iostream>


GPixel getPix(const GColor& color) {
    float alpha = color.a * 255.0f;
    int red = GRoundToInt(color.r * alpha);
    int green = GRoundToInt(color.g * alpha);
    int blue = GRoundToInt(color.b * alpha);

    GPixel pix = GPixel_PackARGB(GRoundToInt(alpha), red, green, blue);
    return pix;
}

bool MyLinearPosGradientShader::isOpaque() {return false;}

bool MyLinearPosGradientShader::setContext(const GMatrix& ctm) {
    // difference vector between p1 and p0
    GVector parallel = {p1.x-p0.x, p1.y-p0.y};
    // perpendicular to parallel vector
    GVector ppend = {-p1.y+p0.y, p1.x-p0.x};
    GVector origin = {p0.x, p0.y};
    // matrix that maps from unit length to original line
    GMatrix mapLine = GMatrix(parallel, ppend, origin);

    std::optional<GMatrix> ctmInverse = ctm.invert();
    std::optional<GMatrix> lineInverse = mapLine.invert();

    if (!(ctmInverse && lineInverse)) {return false;}
    
    posInverse = GMatrix::Concat(lineInverse.value(), ctmInverse.value());
    return true;
}


// Returns a new type of linear gradient. In this variant, the "count" colors are
//      *  positioned along the line p0...p1 not "evenly", but according to the pos[] array.
//      *
//      *  pos[] holds "count" values, each 0...1, which specify the percentage along the
//      *  line where the each color lies.
//      *
//      *  e.g. pos[] = {0, 0.25, 1} would mean 3 colors positioned as follows:
//      *
//      *      color[0] ..... color[1] ..... ..... ..... color[2]
//      *
//      *  color[i] is positioned by computing (1 - pos[i])*p0 + pos[i]*p1
void MyLinearPosGradientShader::shadeRow(int x, int y, int count, GPixel row[]) {
    // float xi = gradInverse[0]*x + gradInverse[2]*y + gradInverse[4];
    // float delta_x = gradInverse[0];

    for (int i = 0; i < count; i++) {
        GPoint pt = {x + i + 0.5f, y + 0.5f};
        posInverse.GMatrix::mapPoints(&pt, 1);

        float xi_clamped = std::max(0.0f, pt.x);
        xi_clamped = std::min(xi_clamped, 1.0f);
        
        //xi_clamped *= scalingFactor;
        int leftColorInd = -2;
        int rightColorInd = -1;
        float distToRight = 0.0f;
        float distToLeft = 0.0f;
        
        if (xi_clamped == 0) {
            //std::cout << "j";
            leftColorInd = 0;
            rightColorInd = 1;
            distToRight = 1.0f;
            distToLeft = 0.0f;
        }
        else if (xi_clamped == 1) {
            //std::cout << "h";
            leftColorInd = positions.size() - 2;
            rightColorInd = positions.size() - 1;
            distToRight = 0.0f;
            distToLeft = 1.0f;
        }
        else {
            //std::cout << "i";
            for (int j = 0; j < positions.size(); j++) {
                if (positions[j] >= xi_clamped) {
                    leftColorInd = j-1;
                    rightColorInd = j;
                    distToLeft = xi_clamped - positions[j];
                    distToRight = positions[j] - xi_clamped;
                    break;
                }
            }
        }

        // determine color from two nearest colors
        float leftNorm = distToRight / (distToLeft + distToRight);
        float rightNorm = distToLeft / (distToLeft + distToRight);
        GColor color = leftNorm*colors[leftColorInd] + rightNorm*colors[rightColorInd];
        
        row[i] = getPix(color);

    }
}



