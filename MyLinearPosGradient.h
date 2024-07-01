#include "include/GShader.h"
#include "include/GMatrix.h"

class MyLinearPosGradientShader : public GShader {
public:
    MyLinearPosGradientShader(GPoint p0, GPoint p1, std::vector<float> positions, std::vector<GColor> colors, int count) :
        p0(p0),
        p1(p1),
        posInverse(GMatrix()),
        positions(positions),
        colors(colors) {}
        

    bool isOpaque() override;
    bool setContext(const GMatrix& ctm) override;
    void shadeRow(int x, int y, int count, GPixel row[]) override;
    

private:
    GPoint p0;
    GPoint p1;
    GMatrix posInverse;
    std::vector<float> positions;
    std::vector<GColor> colors;
};