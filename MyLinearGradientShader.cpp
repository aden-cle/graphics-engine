
#include "MyLinearGradientShader.h"

#include <iostream>


GPixel colortoPix(const GColor& color) {
    float alpha = color.a * 255.0f;
    int red = GRoundToInt(color.r * alpha);
    int green = GRoundToInt(color.g * alpha);
    int blue = GRoundToInt(color.b * alpha);

    GPixel pix = GPixel_PackARGB(GRoundToInt(alpha), red, green, blue);
    return pix;
}


bool MyLinearGradientShader::isOpaque() {
    for (int i = 0; i < colors.size(); i++) {
        if (colors[i].a != 1.0f) {
            return false;
        }
    }
    return true;
}


bool MyLinearGradientShader::setContext(const GMatrix& ctm) {
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
    
    gradInverse = GMatrix::Concat(lineInverse.value(), ctmInverse.value());
    return true;
}


void MyLinearGradientShader::shadeRow(int x, int y, int count, GPixel row[]) {
    // equivalent to 3D vector(x,y,1)
    float xi = gradInverse[0]*x + gradInverse[2]*y + gradInverse[4];
    float delta_x = gradInverse[0];
    
    assert(colors.size() > 1);

    for (int i = 0; i < count; i++) {
        float xi_clamped;
        switch (tileMode) {
            case GTileMode::kClamp:
                // clamp
                xi_clamped = std::max(0.0f, xi);
                xi_clamped = std::min(xi_clamped, 1.0f);
                break;
            case GTileMode::kRepeat:
                // repeat
                xi_clamped = xi - floorf(xi);
                break;
            case GTileMode::kMirror:
                // mirror
                xi_clamped = xi - 2*floorf(xi*0.5f);
                xi_clamped = 1 - fabs(xi_clamped - 1);
                
                break;
        }
        
        xi_clamped *= scalingFactor;
        // get color index
        int k = GFloorToInt(xi_clamped);
        // get distance from color to the left
        float t = xi_clamped - k;
        // determine color from two nearest colors
        GColor color = (1-t)*colors[k] + t*colors[k+1];
        row[i] = colortoPix(color);
        // update along basis vector
        xi += delta_x;
    }
}


std::unique_ptr<GShader> GCreateLinearGradient(
    GPoint p0, GPoint p1,
    const GColor colors[], int count,
    GTileMode tileMode
) {
    // vector constructor creates array with count elements
    std::vector<GColor> vecColors;
    for(int i = 0; i < count; i++) {
        vecColors.push_back(colors[i].pinToUnit());
    }
    // extra element for weighted average calc in shadeRow
    vecColors.push_back({1, 0, 0, 1});
    int scalingFactor = count-1;
    return std::unique_ptr<GShader>(new MyLinearGradientShader(
        p0,
        p1,
        vecColors, 
        scalingFactor, 
        tileMode
    ));
}
