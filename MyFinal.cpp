#include "MyFinal.h"


std::unique_ptr<GShader> MyFinal::createVoronoiShader(const GPoint points[],
                                                         const GColor colors[],
                                                         int count) 
{
    std::vector<GPoint> vecPoints;
    for (int i = 0; i < count; i++) {
        vecPoints.push_back(points[i]);
    }

    std::vector<GColor> vecColors;
    for (int i = 0; i < count; i++) {
        vecColors.push_back(colors[i]);
    }

    return std::unique_ptr<GShader>(new MyVoronoiShader(
        vecPoints,
        vecColors,
        count
    ));
}


// Returns a new type of linear gradient. In this variant, the "count" colors are
// positioned along the line p0...p1 not "evenly", but according to the pos[] array.
//      *
//      *  pos[] holds "count" values, each 0...1, which specify the percentage along the
//      *  line where the each color lies.
//      *
//      *  e.g. pos[] = {0, 0.25, 1} would mean 3 colors positioned as follows:
//      *
//      *      color[0] ..... color[1] ..... ..... ..... color[2]
//      *
//      *  color[i] is positioned by computing (1 - pos[i])*p0 + pos[i]*p1
std::unique_ptr<GShader> MyFinal::createLinearPosGradient(GPoint p0, GPoint p1,
                                                             const GColor colors[],
                                                             const float pos[],
                                                             int count) 
{
    std::vector<float> vecPositions;
    for (int i = 0; i < count; i++) {
        vecPositions.push_back(pos[i]);
    }

    std::vector<GColor> vecColors;
    for (int i = 0; i < count; i++) {
        vecColors.push_back(colors[i]);
    }
    vecColors.push_back({0, 0, 0, 1});

    return std::unique_ptr<GShader>(new MyLinearPosGradientShader(
        p0,
        p1,
        vecPositions,
        vecColors,
        count
    ));;
}

std::unique_ptr<GFinal> GCreateFinal() {
    return std::unique_ptr<GFinal>(new MyFinal);
}
