#include "MyCanvas.h"


GPixel colToPix(const GColor& color) {
    float alpha = color.a * 255.0f;
    int red = GRoundToInt(color.r * alpha);
    int green = GRoundToInt(color.g * alpha);
    int blue = GRoundToInt(color.b * alpha);

    GPixel pix = GPixel_PackARGB(GRoundToInt(alpha), red, green, blue);
    return pix;
}

// inefficient option for composeShader multiplication
GColor pixToCol(const GPixel& pix) {
    GColor c = {
        (float)GPixel_GetR(pix), // 0-255
        (float)GPixel_GetG(pix),
        (float)GPixel_GetB(pix),
        (float)GPixel_GetA(pix),
    };
    float invAlpha255 = 1.0f / c.a;
    c.r *= invAlpha255;
    c.g *= invAlpha255;
    c.b *= invAlpha255;
    c.a /= 255.0f;
    
    return c;
}

unsigned divBy255(unsigned value) {
    return (value + 128) * 257 >> 16;
}


bool MyCanvas::ProxyShader::isOpaque() {
    return baseShader.isOpaque();
}
bool MyCanvas::ProxyShader::setContext(const GMatrix& ctm) {
    return baseShader.setContext(GMatrix::Concat(ctm, drawingToTexsMat));
}
void MyCanvas::ProxyShader::shadeRow(int x, int y, int count, GPixel row[]) {
    baseShader.shadeRow(x, y, count, row);
}




bool MyCanvas::BarycentricShader::isOpaque() {
    for (int i = 0; i < 3; i++) {
        if (colors[i].a != 1) {
            return false;
        }
    }
    return true;
}

bool MyCanvas::BarycentricShader::setContext(const GMatrix& ctm) {
    
    std::optional<GMatrix> ctmInverse = ctm.invert();
    std::optional<GMatrix> drawingToUnitMatInverse = drawingToUnitMat.invert();

    if (!(ctmInverse && drawingToUnitMatInverse)) {return false;}
    
    //inv = ( CTM * M )-1
    storedInverse = GMatrix::Concat(drawingToUnitMatInverse.value(), ctmInverse.value());
    return true;

}

void MyCanvas::BarycentricShader::shadeRow(int x, int y, int count, GPixel row[]) {
    GColor diff1 = colors[1] - colors[0];
    GColor diff2 = colors[2] - colors[0];

    GPoint pt = {x + 0.5f, y + 0.5f};
    storedInverse.mapPoints(&pt, 1);
    GColor baryColor = pt.x*diff1 + pt.y*diff2 + colors[0];
    GColor update = storedInverse[0]*diff1 + storedInverse[1]*diff2;
    
    for (int i = 0; i < count; ++i) {
        baryColor = {
            std::max(0.0f, std::min(1.0f, baryColor.r)),
            std::max(0.0f, std::min(1.0f, baryColor.g)),
            std::max(0.0f, std::min(1.0f, baryColor.b)),
            std::max(0.0f, std::min(1.0f, baryColor.a))
        };
        row[i] = colToPix(baryColor);
        baryColor = baryColor + update;
    }
}




bool MyCanvas::ComposeShader::isOpaque() {
    return shader1.isOpaque() && shader2.isOpaque();
}

bool MyCanvas::ComposeShader::setContext(const GMatrix& ctm) {
    bool context1 = shader1.setContext(ctm); 
    bool context2 = shader2.setContext(ctm);
    return context1 && context2;
}

void MyCanvas::ComposeShader::shadeRow(int x, int y, int count, GPixel row[]) {
    GPixel* row1 = new GPixel[count];
    GPixel* row2 = new GPixel[count];

    shader1.shadeRow(x, y, count, row1);
    shader2.shadeRow(x, y, count, row2);

    for (int i = 0; i < count; i++) {
        // TODO: figure out how to efficiently mutliple color vals
        row[i] = colToPix(pixToCol(row1[i]) * pixToCol(row2[i]));
    }
    
    delete[] row1;
    delete[] row2;
}


