#ifndef __EDGE_H
#define __EDGE_H

#include "include/GPoint.h"


class Edge {
    
public:
    Edge(const GPoint& start, const GPoint& end) :
        start(start),
        end(end)
        {
            if (start.y < end.y) {
                lowerPoint = start;
                higherPoint = end;
                winding = 1;
            }
            else {
                lowerPoint = end;
                higherPoint = start;
                winding = -1;
            }
            recipSlope = (end.x - start.x) / (end.y - start.y);
        }
    
    GPoint start;
    GPoint end;
    // lower in value, not on the screen
    GPoint lowerPoint;
    GPoint higherPoint;
    float recipSlope;
    int winding;
    
    /* Return the x-coord on this line that corresponds to the given y-coord */
    float computeX(float yCoord) const {
        return (yCoord-lowerPoint.y) * recipSlope + lowerPoint.x;
    }

    bool isValid(int y) const {
        return y >= GRoundToInt(lowerPoint.y) && y < GRoundToInt(higherPoint.y);
    }
};


#endif
