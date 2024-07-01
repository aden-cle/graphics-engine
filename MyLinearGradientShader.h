
#include "include/GShader.h"
#include "include/GMatrix.h"

class MyLinearGradientShader : public GShader {
public:
    MyLinearGradientShader(GPoint p0, GPoint p1, std::vector<GColor> colors, int scalingFactor, GTileMode tileMode) :
        p0(p0),
        p1(p1),
        gradInverse(GMatrix()),
        colors(colors), 
        scalingFactor(scalingFactor),
        tileMode(tileMode) {}
        

    bool isOpaque() override;
    bool setContext(const GMatrix& ctm) override;
    void shadeRow(int x, int y, int count, GPixel row[]) override;
    

private:
    GPoint p0;
    GPoint p1;
    GMatrix gradInverse;
    std::vector<GColor> colors;
    int scalingFactor;
    GTileMode tileMode;
};

