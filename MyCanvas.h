/*
 *  Copyright 2024 Aden Clemente
 */

#ifndef _g_starter_canvas_h_
#define _g_starter_canvas_h_

#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GBitmap.h"
#include "include/GTypes.h"
#include "include/GShader.h"
#include "include/GPath.h"
#include "Edge.h"
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <unordered_map>
#include <functional>
#include <numeric>
#include <cmath>

class MyCanvas : public GCanvas {
public:
    MyCanvas(const GBitmap& device) : fDevice(device) {
        matrixStack.push_back(GMatrix());
    }

    void drawPath(const GPath& path, const GPaint& paint) override;
    void save() override;
    void restore() override;
    void concat(const GMatrix& matrix) override;
    void clear(const GColor& color) override;
    void drawRect(const GRect& rect, const GPaint& paint) override;
    void drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) override;
    void drawMesh(
        const GPoint verts[], const GColor colors[], const GPoint texs[],
        int count, const int indices[], const GPaint&
    ) override;
    void drawQuad(
        const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
        int level, const GPaint&
    ) override;
    void fillRect(const GRect& rect, const GColor& color);

private:
    // We store a copy of the bitmap
    const GBitmap fDevice;  
    // Store our stack of transformation matrices
    std::vector<GMatrix> matrixStack;
    // Add whatever other fields you need
    
    
    class ProxyShader : public GShader {
    public:
        ProxyShader(GShader& baseShader, const GMatrix& drawingToTexsMat) : baseShader(baseShader), drawingToTexsMat(drawingToTexsMat) {}
        
        bool isOpaque() override;
        bool setContext(const GMatrix& ctm) override;
        void shadeRow(int x, int y, int count, GPixel row[]) override;
    private:
        GShader& baseShader;
        GMatrix drawingToTexsMat;
    };

    class BarycentricShader : public GShader {
    public:
        BarycentricShader(const GMatrix& drawingToUnitMat, const GColor colors[])
            : drawingToUnitMat(drawingToUnitMat), storedInverse(GMatrix()), colors(colors) {}
        
        bool isOpaque() override;
        bool setContext(const GMatrix& ctm) override;
        void shadeRow(int x, int y, int count, GPixel row[]) override;
        
    private:
        const GMatrix drawingToUnitMat;
        GMatrix storedInverse;
        // colors is a pointer - the array must outlive BarycentricShader
        const GColor* colors;
    };

    class ComposeShader : public GShader {
    public:
        ComposeShader(GShader& shader1, GShader& shader2) : shader1(shader1), shader2(shader2) {}
        
        bool isOpaque() override;
        bool setContext(const GMatrix& ctm) override;
        void shadeRow(int x, int y, int count, GPixel row[]) override;
        
    private:
        GShader& shader1;
        GShader& shader2;
    };
    
};

std::string GDrawSomething(GCanvas* canvas, GISize dim);

class MyShader : public GShader {
public:
    MyShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode tileMode)
        : bm(bitmap), lMatrix(localMatrix), inverse(GMatrix()), tileMode(tileMode) {}

    bool isOpaque() override;
    bool setContext(const GMatrix& ctm) override;
    void shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    const GBitmap bm;
    const GMatrix lMatrix;
    GMatrix inverse;
    GTileMode tileMode;
};

#endif
