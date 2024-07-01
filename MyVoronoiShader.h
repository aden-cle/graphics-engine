#include "include/GShader.h"
#include "include/GMatrix.h"


class MyVoronoiShader : public GShader {
public:
    MyVoronoiShader(std::vector<GPoint> points, std::vector<GColor> colors, int count) :
        points(points), 
        colors(colors),
        vInverse(GMatrix()) {}
        
    bool isOpaque() override;
    bool setContext(const GMatrix& ctm) override;
    void shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    std::vector<GPoint> points;
    std::vector<GColor> colors;
    GMatrix vInverse;
};