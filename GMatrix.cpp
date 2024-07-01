#include "include/GMatrix.h"
#include <optional>
#include <cmath>


GMatrix::GMatrix() {
    fMat[0] = 1.0;    fMat[2] = 0.0;    fMat[4] = 0.0;
    fMat[1] = 0.0;    fMat[3] = 1.0;    fMat[5] = 0.0;
}


GMatrix GMatrix::Translate(float tx, float ty) {
    return GMatrix(1.0f, 0, tx, 0, 1.0f, ty);
}

GMatrix GMatrix::Scale(float sx, float sy) {
    return GMatrix(sx, 0, 0, 0, sy, 0);
}

GMatrix GMatrix::Rotate(float radians) {
    float sint = std::sin(radians);
    float cost = std::cos(radians);
    return GMatrix(cost, -sint, 0, sint, cost, 0);
}

GMatrix GMatrix::Concat(const GMatrix& a, const GMatrix& b) {
    GMatrix prod = GMatrix(0, 0, 0, 0, 0, 0);

    prod[0] = a[0]*b[0] + a[2]*b[1];
    prod[1] = a[1]*b[0] + a[3]*b[1];
    prod[2] = a[0]*b[2] + a[2]*b[3];
    prod[3] = a[1]*b[2] + a[3]*b[3];
    prod[4] = a[0]*b[4] + a[2]*b[5] + a[4];
    prod[5] = a[1]*b[4] + a[3]*b[5] + a[5];

    return prod;
}

std::optional<GMatrix> GMatrix::invert() const {
    float denom = fMat[0]*fMat[3] - fMat[1]*fMat[2];
    
    // if inverse does not exist, return nothing
    if (denom == 0) {
        return {};
    }
    
    GMatrix inverse = GMatrix(fMat[3], -fMat[2], fMat[2]*fMat[5] - fMat[3]*fMat[4], -fMat[1], fMat[0], -fMat[0]*fMat[5] + fMat[1]*fMat[4]);
    float coeff = 1.0f / denom;
    
    for (int i = 0; i < 6; i++) {
        inverse[i] *= coeff;
    }

    return inverse;
}

void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const {
    for (int i = 0; i < count; i++) {
        float sx = src[i].x;
        float sy = src[i].y;
        float xi = fMat[0]*sx + fMat[2]*sy + fMat[4];
        float yi = fMat[1]*sx + fMat[3]*sy + fMat[5];
        dst[i] = {xi, yi};
    }
}