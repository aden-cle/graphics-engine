#include "include/GFinal.h"
#include "MyVoronoiShader.h"
#include "MyLinearPosGradient.h"

class MyFinal : public GFinal {
public:

    std::unique_ptr<GShader> createVoronoiShader(const GPoint points[],
                                                         const GColor colors[],
                                                         int count) override;
    
    std::unique_ptr<GShader> createLinearPosGradient(GPoint p0, GPoint p1,
                                                             const GColor colors[],
                                                             const float pos[],
                                                             int count) override;

};
