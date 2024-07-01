#include "include/GPath.h"
#include <iostream>


void GPath::addRect(const GRect& r, Direction dir) {
    GPoint pts[4];
    pts[0] = {r.left, r.top};
    pts[1] = {r.right, r.top};
    pts[2] = {r.right, r.bottom};
    pts[3] = {r.left, r.bottom};

    moveTo(pts[0]);

    if (dir == Direction::kCW_Direction) {
        for (int i = 1; i < 4; i++) {
            lineTo(pts[i]);
        }
    } 
    else {
        for (int i = 3; i > 0; i--) {
            lineTo(pts[i]);
        }
    }
}


void GPath::addPolygon(const GPoint pts[], int count) {
    moveTo(pts[0]);
    for (int i = 1; i < count; i++) {
        lineTo(pts[i]);
    }
}


/**
    *  Append a new contour respecting the Direction. The contour should be an approximate
    *  circle (8 quadratic curves will suffice) with the specified center and radius.
    */
void GPath::addCircle(GPoint center, float radius, Direction dir) {    
    float sin45 = std::sqrt(2.0f) / 2.0f;

    GPoint at0deg { radius, 0 };
    GPoint at45deg { radius*sin45, radius*sin45 };
    // intersection of y - sin45 = -1(x - sin45) with x=1
    // y = sin45 + -1 + sin45
    GPoint control0to45 { radius * 1, radius * (2*sin45-1) };

    GPoint quadCCWEighth[] = {control0to45, at45deg};
    GPoint quadCWEighth[] = {control0to45, at0deg};

    GMatrix rot45CCW = GMatrix::Rotate(gFloatPI / 4.0f);
    GMatrix rot45CW = GMatrix::Rotate(-gFloatPI / 4.0f);
    
    if (dir != kCCW_Direction) {
        moveTo(at0deg + center);

        for (int i = 0; i < 8; i++) {
            quadTo(quadCCWEighth[0] + center, quadCCWEighth[1] + center);
            rot45CCW.GMatrix::mapPoints(quadCCWEighth, 2);
        }
    }
    else {
        moveTo(at45deg + center);
        
        for (int i = 0; i < 8; i++) {
            quadTo(quadCWEighth[0] + center, quadCWEighth[1] + center);
            rot45CW.GMatrix::mapPoints(quadCWEighth, 2);
        }
    }
}


GRect GPath::bounds() const {
    if (fPts.size() == 0) {return GRect::LTRB(0.0f, 0.0f, 0.0f, 0.0f);}

    float leftBound = std::numeric_limits<float>::max();
    float topBound = std::numeric_limits<float>::max();
    float rightBound = std::numeric_limits<float>::min();
    float bottomBound = std::numeric_limits<float>::min();

    for (int i = 0; i < fPts.size(); i++) {
        float p_x = fPts[i].x;
        float p_y = fPts[i].y;
        leftBound = std::min(leftBound, p_x);
        topBound = std::min(topBound, p_y);
        rightBound = std::max(rightBound, p_x);
        bottomBound = std::max(bottomBound, p_y);
    }

    return GRect::LTRB(leftBound, topBound, rightBound, bottomBound);
}

/*
*  Transform the path in-place by the specified matrix.
*/
void GPath::transform(const GMatrix& m) {    
    GPoint* ptr_to_first = fPts.data();
    m.mapPoints(ptr_to_first, fPts.size());
}


/**
    *  Given 0 < t < 1, subdivide the src[] quadratic bezier at t into two new quadratics in dst[]
    *  such that
    *  0...t is stored in dst[0..2]
    *  t...1 is stored in dst[2..4]
    */
void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t) {
    GPoint A = src[0], B = src[1], C = src[2];
    
    // make indexes more readable
    GPoint& d1A = dst[0];
    GPoint& d1B = dst[1];
    GPoint& d1C = dst[2];
    // GPoint& d2A
    GPoint& d2B = dst[3];
    GPoint& d2C = dst[4];
    
    assert(0 < t && t < 1);
    
    d1A = A;
    d2C = C;
    
    // evaluate src at t
    float q = 1-t;
    d1C = q*q*A + 2*t*q*B + t*t*C;
    
    // control points
    d1B = A + (B - A) * t;
    d2B = B + (C - B) * t;
    
    
    // https://stackoverflow.com/questions/37082744/split-one-quadratic-bezier-curve-into-two
    // Answered May 7, 2016 at 4:24
    // User Blindman67
}

/**
    *  Given 0 < t < 1, subdivide the src[] cubic bezier at t into two new cubics in dst[]
    *  such that
    *  0...t is stored in dst[0..3]
    *  t...1 is stored in dst[3..6]
    */
void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t) {
    assert(0 < t && t < 1);
    
    GPoint A = src[0], B = src[1], C = src[2], D = src[3];
    
    GPoint& d1A = dst[0];
    GPoint& d1B = dst[1];
    GPoint& d1C = dst[2];
    GPoint& d1D = dst[3]; // same as d2A
    GPoint& d2B = dst[4];
    GPoint& d2C = dst[5];
    GPoint& d2D = dst[6];
    
    d1A = A;
    d2D = D;
    
    float q = 1-t;
    
    // first tangents in forward direction and reverse
    d1B = A + (B-A)*t;
    d2C = D + (C-D)*q;
    
    // second tangents like p123 in forward direction and reverse
    GPoint BC = B + (C-B)*t;
    GPoint CB = C + (B-C)*q;
    d1C = d1B + (BC-d1B)*t;
    d2B = d2C + (CB-d2C)*q;
    
    // evaluate cubic ABCD at t
    d1D = q*q*q*A + 3*q*q*t*B + 3*q*t*t*C + t*t*t*D;
    
    
    // https://stackoverflow.com/questions/8369488/splitting-a-bezier-curve
    // Answered Dec 6, 2011 at 19:40
    // User: Jonathan
}

