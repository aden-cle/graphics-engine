/*
 *  Copyright 2024 Aden Clemente
 */

#include "MyCanvas.h"




GPixel colortoPixel(const GColor& color) {
    float alpha = color.a * 255.0f;
    int red = GRoundToInt(color.r * alpha);
    int green = GRoundToInt(color.g * alpha);
    int blue = GRoundToInt(color.b * alpha);

    GPixel pix = GPixel_PackARGB(GRoundToInt(alpha), red, green, blue);
    return pix;
}


GPixel clearPixel(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* Return a clear pixel */
    return GPixel_PackARGB(0, 0, 0, 0);
}


GPixel kSrc(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    return src;
}


GPixel kDst(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    return *dst;
}


unsigned div255(unsigned value) {
    return (value + 128) * 257 >> 16;
}


GPixel srcOver(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* Given a background pixel, dst, a foreground pixel, src,
    in premul form, and a foreground alpha, in the range [0-1],
    return the blended pixel */
    /* S + (1 - Sa)*D */

    int multiplier = 255 - src_a; 

    int alpha = div255(multiplier * GPixel_GetA(*dst));
    int red = div255(multiplier * GPixel_GetR(*dst));
    int green = div255(multiplier * GPixel_GetG(*dst));
    int blue = div255(multiplier * GPixel_GetB(*dst));

    GPixel altered_dst = GPixel_PackARGB(alpha, red, green, blue);
    return src + altered_dst;
}


GPixel dstOver(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* D + (1 - Da)*S
    */

    if (dst_a == 255) {return *dst;}

    int multiplier = 255 - dst_a; 

    int alpha = div255(multiplier * GPixel_GetA(src));
    int red = div255(multiplier * GPixel_GetR(src));
    int green = div255(multiplier * GPixel_GetG(src));
    int blue = div255(multiplier * GPixel_GetB(src));

    GPixel altered_src = GPixel_PackARGB(alpha, red, green, blue);
    return *dst + altered_src;
}


GPixel srcIn(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* Da * S
    */

    if (dst_a == 255) {return src;}
    if (dst_a == 0) {return GPixel_PackARGB(0, 0, 0, 0);}

    int alpha = div255(dst_a * GPixel_GetA(src));
    int red = div255(dst_a * GPixel_GetR(src));
    int green = div255(dst_a * GPixel_GetG(src));
    int blue = div255(dst_a * GPixel_GetB(src));

    GPixel altered_src = GPixel_PackARGB(alpha, red, green, blue);
    return altered_src;
}


GPixel dstIn(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* Sa * D
    */

    int alpha = div255(src_a * GPixel_GetA(*dst));
    int red = div255(src_a * GPixel_GetR(*dst));
    int green = div255(src_a * GPixel_GetG(*dst));
    int blue = div255(src_a * GPixel_GetB(*dst));

    GPixel altered_dst = GPixel_PackARGB(alpha, red, green, blue);
    return altered_dst;
}


GPixel srcOut(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* (1 - Da)*S
    */

    if (dst_a == 255) {return GPixel_PackARGB(0, 0, 0, 0);}
    if (dst_a == 0) {return src;}

    int multiplier = 255 - dst_a; 

    int alpha = div255(multiplier * GPixel_GetA(src));
    int red = div255(multiplier * GPixel_GetR(src));
    int green = div255(multiplier * GPixel_GetG(src));
    int blue = div255(multiplier * GPixel_GetB(src));

    GPixel altered_dst = GPixel_PackARGB(alpha, red, green, blue);
    return altered_dst;
}


GPixel dstOut(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* (1 - Sa)*D
    */

    int multiplier = 255 - src_a; 

    int alpha = div255(multiplier * GPixel_GetA(*dst));
    int red = div255(multiplier * GPixel_GetR(*dst));
    int green = div255(multiplier * GPixel_GetG(*dst));
    int blue = div255(multiplier * GPixel_GetB(*dst));

    GPixel altered_dst = GPixel_PackARGB(alpha, red, green, blue);
    return altered_dst;
}


GPixel srcATop(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* Da*S + (1 - Sa)*D 
    */ 

    if (dst_a == 0) {return dstOut(src, dst, src_a, dst_a);}

    int mult_a = 255 - src_a; 
    int mult_b = dst_a;

    int alpha_a = div255(mult_a * GPixel_GetA(*dst));
    int red_a = div255(mult_a * GPixel_GetR(*dst));
    int green_a = div255(mult_a * GPixel_GetG(*dst));
    int blue_a = div255(mult_a * GPixel_GetB(*dst));

    int alpha_b = div255(mult_b * GPixel_GetA(src));
    int red_b = div255(mult_b * GPixel_GetR(src));
    int green_b = div255(mult_b * GPixel_GetG(src));
    int blue_b = div255(mult_b * GPixel_GetB(src));

    GPixel altered_dst = GPixel_PackARGB(alpha_a, red_a, green_a, blue_a);
    GPixel altered_src = GPixel_PackARGB(alpha_b, red_b, green_b, blue_b);
    return altered_src + altered_dst;
}


GPixel dstATop(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* Sa*D + (1 - Da)*S
    */

    if (dst_a == 255) {return dstIn(src, dst, src_a, dst_a);}

    int mult_a = src_a; 
    int mult_b = 255 - dst_a;

    int alpha_a = div255(mult_a * GPixel_GetA(*dst));
    int red_a = div255(mult_a * GPixel_GetR(*dst));
    int green_a = div255(mult_a * GPixel_GetG(*dst));
    int blue_a = div255(mult_a * GPixel_GetB(*dst));

    int alpha_b = div255(mult_b * GPixel_GetA(src));
    int red_b = div255(mult_b * GPixel_GetR(src));
    int green_b = div255(mult_b * GPixel_GetG(src));
    int blue_b = div255(mult_b * GPixel_GetB(src));

    GPixel altered_dst = GPixel_PackARGB(alpha_a, red_a, green_a, blue_a);
    GPixel altered_src = GPixel_PackARGB(alpha_b, red_b, green_b, blue_b);
    return altered_src + altered_dst;
}


GPixel blendXor(const GPixel src, GPixel* dst, int src_a, int dst_a) {
    /* (1 - Sa)*D + (1 - Da)*S
    */

    if (dst_a == 255) {return dstOut(src, dst, src_a, dst_a);}

    int mult_a = 255 - src_a; 
    int mult_b = 255 - dst_a;

    int alpha_a = div255(mult_a * GPixel_GetA(*dst));
    int red_a = div255(mult_a * GPixel_GetR(*dst));
    int green_a = div255(mult_a * GPixel_GetG(*dst));
    int blue_a = div255(mult_a * GPixel_GetB(*dst));

    int alpha_b = div255(mult_b * GPixel_GetA(src));
    int red_b = div255(mult_b * GPixel_GetR(src));
    int green_b = div255(mult_b * GPixel_GetG(src));
    int blue_b = div255(mult_b * GPixel_GetB(src));

    GPixel altered_dst = GPixel_PackARGB(alpha_a, red_a, green_a, blue_a);
    GPixel altered_src = GPixel_PackARGB(alpha_b, red_b, green_b, blue_b);
    return altered_src + altered_dst;
}


void MyCanvas::clear(const GColor& color) {
    GPixel pix = colortoPixel(color);

    visit_pixels(fDevice, [&](int x, int y, GPixel* p){
        *p = pix;
    });
}


typedef GPixel (*BlendFunc)(const GPixel, GPixel*, int, int);

BlendFunc getBlendFunc(GBlendMode blendMode) {
    static const std::unordered_map<GBlendMode, BlendFunc> blendFuncMap = {
        {GBlendMode::kClear, &clearPixel},
        {GBlendMode::kSrc, &kSrc},
        {GBlendMode::kDst, &kDst},
        {GBlendMode::kSrcOver, &srcOver},
        {GBlendMode::kDstOver, &dstOver},
        {GBlendMode::kSrcIn, &srcIn},
        {GBlendMode::kDstIn, &dstIn},
        {GBlendMode::kSrcOut, &srcOut},
        {GBlendMode::kDstOut, &dstOut},
        {GBlendMode::kSrcATop, &srcATop},
        {GBlendMode::kDstATop, &dstATop},
        {GBlendMode::kXor, &blendXor}
    };
    auto it = blendFuncMap.find(blendMode);
    if (it != blendFuncMap.end()) {
        return it->second;
    }
    else {
        std::cerr << "Error: Invalid blend mode" << std::endl;
        return &srcOver;
    }  
}


GBlendMode optimizeMode(GBlendMode mode, int src_a) {
    if (mode == GBlendMode::kSrcOver && src_a == 255) {
        return GBlendMode::kSrc;
    } else if (mode == GBlendMode::kDstIn && src_a == 0) {
        return GBlendMode::kClear;
    } else if (mode == GBlendMode::kDstIn && src_a == 255) {
        return GBlendMode::kDst;
    } else if (mode == GBlendMode::kDstOut && src_a == 0) {
        return GBlendMode::kDst;
    } else if (mode == GBlendMode::kDstOut && src_a == 255) {
        return GBlendMode::kClear;
    } else if (mode == GBlendMode::kSrcATop && src_a == 255) {
        return GBlendMode::kSrcIn;
    } else if (mode == GBlendMode::kDstATop && src_a == 0) {
        return GBlendMode::kSrcOut;
    } else if (mode == GBlendMode::kXor && src_a == 255) {
        return GBlendMode::kSrcOut;
    } else {
        return mode;
    }
}


void MyCanvas::save() {
    GMatrix ctm = matrixStack.back();
    matrixStack.push_back(ctm);
}


void MyCanvas::restore() {
    matrixStack.pop_back();
}


void MyCanvas::concat(const GMatrix& matrix) {
    GMatrix ctm = matrixStack.back();
    GMatrix new_tm = GMatrix::Concat(ctm, matrix);

    matrixStack[matrixStack.size()-1] = new_tm;
}



std::ostream &operator<<(std::ostream &os, Edge const &e) { 
    os << "lower xy: " << e.lowerPoint.x << ", " << e.lowerPoint.y;
    os << "; higher xy: " << e.higherPoint.x << ", " << e.higherPoint.y;
    os << "; recipslope: " << e.recipSlope << ", winding: " << e.winding << std::endl;
    return os;
}

bool validPoint(const GPoint& point, int deviceWidth, int deviceHeight) {
    return point.y >= 0 && point.y < deviceHeight;
}


void MyCanvas::drawPath(const GPath& path, const GPaint& paint) {
    // Access transformation matrix from top of stack 
    GMatrix ctm = matrixStack.back();
    // Run path through transformation matrix 
    GPath pathCopy = path;
    pathCopy.transform(ctm);

    GColor color = paint.getColor();
    GPixel src = colortoPixel(color); 
    int src_a = GPixel_GetA(src);

    GBlendMode blendMode = paint.getBlendMode();

    if (paint.getShader() == nullptr && (src_a == 0 || src_a == 255)) {
        blendMode = optimizeMode(blendMode, src_a);
    }

    BlendFunc blendFunc = getBlendFunc(blendMode);
    
    int dh = fDevice.height();
    int dw = fDevice.width();

    // vector to store and sort edges
    std::vector<Edge> edges;
    
    // loop through verbs and points in path using GPath::Edger
    GPath::Edger edger(pathCopy);
    GPoint pts[GPath::kMaxNextPoints];

    float tolerance = 0.25f;
    
    while (auto v = edger.next(pts)) {
        switch (v.value()) {
            case GPath::kLine:
                // pts[0..1]
                // don't consider horizontal edges
                // TODO: don't consider completely out of bounds edges
                if (pts[0].y != pts[1].y) {
                    // Creates an Edge inside the vector at last index
                    // location in memory, avoiding a copy operation
                    // Edge constructor adds winding value 
                    edges.emplace_back(pts[0], pts[1]);
                }
                break;
            case GPath::kQuad: {
                // pts[0..2]
                if (!validPoint(pts[0], dw, dh) && !validPoint(pts[1], dw, dh) && !validPoint(pts[2], dw, dh)) {
                    break;
                }
                GPoint firstTermVec = pts[0] - 2*pts[1] + pts[2];
                GPoint E = firstTermVec * 0.25f;
                int N = GCeilToInt(std::sqrt(E.length() / tolerance));
                GPoint secondTermVec = -2*pts[0] + 2*pts[1];
                
                GPoint p0 = pts[0];        // this is equivalent to curve.evalute(0)
                for (int i = 1; i < N; i++) {
                    float t = i * 1.0f / N;       // 0 < t <= 1
                    GPoint p1;
                    p1 = t*t*firstTermVec + t*secondTermVec + pts[0];

                    if (p0.y != p1.y) {
                        // Creates an Edge inside the vector at last index
                        // location in memory, avoiding a copy operation
                        // Edge constructor adds winding value
                        edges.emplace_back(p0, p1);
                    }
                    p0 = p1;
                }
                edges.emplace_back(p0, pts[2]);
                break;
            }
            case GPath::kCubic: {
                // pts[0..3]
                // E0 = A - 2B + C
                // E1 = B - 2C + D
                // E.x = max(abs(E0.x), abs(E1.x))
                // E.y = max(abs(E0.y), abs(E1.y))
                // int num_segs = ceil(sqrt((3*|E|)/(4*tolerance)));
                if (!validPoint(pts[0], dw, dh) && !validPoint(pts[1], dw, dh) && !validPoint(pts[2], dw, dh) && !validPoint(pts[3], dw, dh)) {
                    break;
                }
                GPoint E0 = pts[0] - 2*pts[1] + pts[2];
                GPoint E1 = pts[1] - 2*pts[2] + pts[3];
                GPoint E {
                    std::max(std::abs(E0.x), std::abs(E1.x)),
                    std::max(std::abs(E0.y), std::abs(E1.y))
                };
                int N = GCeilToInt(std::sqrt(3*E.length() / (4*tolerance)));
                
                GPoint firstTermVec = -1*pts[0] + 3*pts[1] - 3*pts[2] + pts[3];
                GPoint secondTermVec = 3*pts[0] - 6*pts[1] + 3*pts[2];
                GPoint thirdTermVec = -3*pts[0] + 3*pts[1];
                
                GPoint p0 = pts[0];        // this is equivalent to curve.evalute(0)
                for (int i = 1; i < N; i++) {
                    float t = i * 1.0f / N;       // 0 < t <= 1
                    GPoint p1;
                    //p1 = (1-t)*(1-t)*pts[0] + 2*t*(1-t)*pts[1] + t*t*pts[2];
                    p1 = t*t*t * firstTermVec + t*t*secondTermVec + t*thirdTermVec + pts[0];
    
                    if (p0.y != p1.y) {
                        // Creates an Edge inside the vector at last index
                        // location in memory, avoiding a copy operation.
                        // Edge constructor adds winding value 
                        edges.emplace_back(p0, p1);
                    }
                    p0 = p1;
                }
                edges.emplace_back(p0, pts[3]);
                break;
            }
            default:
                // do nothing if v is a move
                break;
        }
    }
    
    GRect bounds = pathCopy.bounds();
    int yMin = std::max(0, GRoundToInt(bounds.top));
    int yMax = std::min(fDevice.height(), GRoundToInt(bounds.bottom));

    // different from pseudocode: sort with respect to x at yMin
    std::sort(
        edges.begin(), edges.end(),
        [yMin](const Edge& A, const Edge& B) {
            return A.computeX(yMin+0.5f) < B.computeX(yMin+0.5f);
        }
    );
    
    // std::stable_sort(
    //     edges.begin(), edges.end(),
    //     [](const Edge& A, const Edge& B) {
    //         return A.lowerPoint.y < B.lowerPoint.y;
    //     }
    // );
    // for (int i =0; i < edges.size(); i++){
    //     std::cout << edges[i];
    // }

    // loop through all y’s containing edges
    for (int y = yMin; y < yMax; y++) {
        size_t i = 0;
        int w = 0;
        int L;

        // loop through active edges for this Y value
        while (i < edges.size()) {
            if (!edges[i].isValid(y)) {
                i++;
                continue;
            }

            // TODO: extract from loop
            int x = GRoundToInt(edges[i].computeX(y+0.5f));

            if (w == 0) {
                L = x;
            }

            w += edges[i].winding;  // +1 or -1
            // int test_x = std::max(0, x);
            // test_x = std::min(test_x, fDevice.width()-1);
            // GPixel* dst_p = fDevice.getAddr(test_x, y);
            // *dst_p = src;
            
            if (w == 0) {
                // blit(L, y, R - L);
                int R = x;

                if (paint.getShader() == nullptr) {
                    for (int x = std::max(0, L); x < std::min(fDevice.width(), R); x++) {
                        GPixel* dst_ptr = fDevice.getAddr(x, y);
                
                        int dst_a = GPixel_GetA(*dst_ptr);
                        GPixel blended_pix = blendFunc(src, dst_ptr, src_a, dst_a);
                        *dst_ptr = blended_pix;
                    }
                }
                else {
                    GShader* shader = paint.getShader();
                    
                    if (shader->isOpaque()) {
                        blendMode = optimizeMode(blendMode, 255);
                    }
                    BlendFunc blendFunc = getBlendFunc(blendMode);

                    shader->setContext(ctm); 
                    // assert(R >= L);
                    // GPixel storage[R - L];
                    std::vector<GPixel> storage(R - L, 0);
                    assert(storage.size() == R - L);
                    shader->shadeRow(L, y, R - L, storage.data());
                    
                    for (int x = std::max(0, L); x < std::min(fDevice.width(), R); x++) {
                        GPixel* dst_ptr = fDevice.getAddr(x, y);
                        GPixel shaderSrc = storage[x - L];
                            
                        int dst_a = GPixel_GetA(*dst_ptr);
                        int src_a = GPixel_GetA(shaderSrc);

                        GPixel blended_pix = blendFunc(shaderSrc, dst_ptr, src_a, dst_a);
                        *dst_ptr = blended_pix;
                    }
                }
            }

            if (edges[i].isValid(y+1)) {
                i += 1;
            } 
            else {
                edges.erase(edges.begin() + i);
            }
        }

        // if (w != 0) {
        //     std::cout << w << " ";
        //     std::cout << "Y val: " << y << " ";
        //     std::cout << "height: "<< fDevice.height() << std::endl;
        // }
        
        assert(w == 0);
        
        // now i is the number of remaining valid edges, so
        // account for any new edges that will be valid for next y
        //while (i < edges.size() && edges[i].isValid(y+1)) {
        //     i += 1;
        //}
        // now i also includes the number of edges that will be valid
        
        assert(edges.size() == i);

        // sort_edges( [0…i) based on computed X for y+1 )
        // currently sorts entire vector
        std::sort(
            edges.begin(), edges.end(),
            [y](const Edge& A, const Edge& B) {
                return A.computeX(y+1.5f) < B.computeX(y+1.5f);
            }
        );
    }
}


// TODO: Deprecate and replace with getX member function of Edge
float getXFromY(float yCoord, float x, float y, float recipSlope) {
        /* Given a y-coordinate, the reciprocal of the line's slope,
        and a point which lies on that line, return the x-coord
        that corresponds to the y-coord */
        return (yCoord - y) * recipSlope + x;
}


GMatrix computeBasis(GPoint p0, GPoint p1, GPoint p2) {
    GVector basis01 = {p1.x-p0.x, p1.y-p0.y};
    GVector basis02 = {p2.x-p0.x, p2.y-p0.y};
    GVector basisOrigin = {p0.x, p0.y};
    return GMatrix(basis01, basis02, basisOrigin);
}


void MyCanvas::drawMesh (
    const GPoint verts[], const GColor colors[], const GPoint texs[],
    int count, const int indices[], const GPaint& paint) {
    
    assert(!(colors == nullptr && texs == nullptr));
    // if texs == null and colors != null    // BarycentricShader(using 3 colors)
    // if texs != null and colors == null    // ProxyShader(paint.shader, 3 texs)
    // if neither null                       // ComposeShader(TextureProxyShader(paint.shader), barycentric(using 3 colors)

    GShader& ogShader = *paint.getShader();

    GPoint points[3];
    GColor triColors[3];
    GPoint triTexs[3];
    
    int n = 0;
    for (int i = 0; i < count; ++i) {
        points[0] = verts[indices[n+0]];
        points[1] = verts[indices[n+1]];
        points[2] = verts[indices[n+2]];

        // matrix that maps unit triangle to this triangle
        GMatrix unitToDraw = computeBasis(points[0], points[1], points[2]);
        
        // will hold the new paint that can be set depending on which shaders are set
        GPaint newPaint;
        
        BarycentricShader* optBary = nullptr;
        if (colors != nullptr) {
            triColors[0] = colors[indices[n+0]];
            triColors[1] = colors[indices[n+1]];
            triColors[2] = colors[indices[n+2]];
            
            optBary = new BarycentricShader(unitToDraw, triColors);
            newPaint = GPaint(optBary);
        }

        ProxyShader* optProxy = nullptr;
        if (texs != nullptr) {
            triTexs[0] = texs[indices[n+0]];
            triTexs[1] = texs[indices[n+1]];
            triTexs[2] = texs[indices[n+2]];
            
            GMatrix unitToTex = computeBasis(triTexs[0], triTexs[1], triTexs[2]);
            std::optional<GMatrix> textToUnit = unitToTex.invert();
            if (!textToUnit) {std::cerr << "Error: Inverse of texture matrix DNE" << std::endl;}
            GMatrix textToDraw = GMatrix::Concat(unitToDraw, textToUnit.value());
            
            optProxy = new ProxyShader(ogShader, textToDraw);
            newPaint = GPaint(optProxy);
        }
        
        ComposeShader* optCompose = nullptr;
        if (optProxy != nullptr && optBary != nullptr) {
            optCompose = new ComposeShader(*optProxy, *optBary);
            newPaint = GPaint(optCompose);
        }
        
        drawConvexPolygon(points, 3, newPaint);
        
        delete optBary;
        delete optProxy;
        delete optCompose;
        
        n += 3;
    }
}


void MyCanvas::drawQuad(
    const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
    int level, const GPaint& paint) {
    
    // Level is just the number of lines between segments:
    int numSegments = level + 1;
    
    // if numSegments is 4, there's 16 parallelograms and 32 triangles
    // 32 triangles -> 32*3 indices
    // 4 segments -> (4+1)^2 unique vertex positions/texs/colors
    int numTris = numSegments*numSegments*2;
    std::vector<int> indices;
    indices.reserve(numTris*3);
    std::vector<GPoint> meshVerts((numSegments+1)*(numSegments+1));
    std::vector<GColor> meshColors((numSegments+1)*(numSegments+1));
    std::vector<GPoint> meshTexs((numSegments+1)*(numSegments+1));
    
    for (int i = 0; i < numSegments; ++i) {
        for (int j = 0; j < numSegments; ++j) {
            // for each little quadrilateral in the numSegments*numSegments grid,
            // add 2 triangles AKA 6 indices
            // length = numSegments*numSegments*6
            
            // first find the indices of the 4 unique corners
            int leftI = i, rightI = i+1;
            int botJ = j, topJ = j+1;
            int iStride = numSegments+1; // like "numCorners"
            int topLeftIdx = topJ + (iStride*leftI);
            int topRigtIdx = topJ + (iStride*rightI);
            int botLeftIdx = botJ + (iStride*leftI);
            int botRigtIdx = botJ + (iStride*rightI);
            
            // first triangle is botLeftIdx,botRigtIdx,topRigtIdx
            // second triangle is botLeftIdx,topRigtIdx,topLeftIdx
            indices.push_back(botLeftIdx);
            indices.push_back(botRigtIdx);
            indices.push_back(topRigtIdx);
            indices.push_back(botLeftIdx);
            indices.push_back(topRigtIdx);
            indices.push_back(topLeftIdx);
        }
    }
    
    
    for (int i = 0; i < numSegments+1; ++i) {
        for (int j = 0; j < numSegments+1; ++j) {
            // one vertex per loop iter
            // length = (ns+1)*(ns+1)
            
            // store in the appropriate index
            int iStride = numSegments+1;
            int vertexIdx = j + (iStride * i);
            
            // interpolate
            float tI = (float)i / numSegments;
            float tJ = (float)j / numSegments;
            
            GPoint vertTop = (1-tI)*verts[0] + tI*verts[1];
            GPoint vertBot = (1-tI)*verts[3] + tI*verts[2];
            GPoint vert = (1-tJ)*vertBot + tJ*vertTop;
            meshVerts[vertexIdx] = vert;
            
            if (colors) {
                // interpolate
                GColor colTop = (1-tI)*colors[0] + tI*colors[1];
                GColor colBot = (1-tI)*colors[3] + tI*colors[2];
                GColor col = (1-tJ)*colBot + tJ*colTop;
                meshColors[vertexIdx] = col;
            }
            if (texs) {
                // interpolate
                GPoint texTop = (1-tI)*texs[0] + tI*texs[1];
                GPoint texBot = (1-tI)*texs[3] + tI*texs[2];
                GPoint tex = (1-tJ)*texBot + tJ*texTop;
                meshTexs[vertexIdx] = tex;
            }
        }
    }
    
    drawMesh(
        meshVerts.data(),
        colors ? meshColors.data() : nullptr,
        texs ? meshTexs.data() : nullptr,
        numTris,
        indices.data(),
        paint
    );
    
}



void MyCanvas::drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) {
    if (count <= 2) {return;}

    // Access transformation matrix from top of stack 
    GMatrix ctm = matrixStack.back();
    // Run points through transformation matrix 
    GPoint mappedPoints[count];
    ctm.mapPoints(mappedPoints, points, count);

    GColor color = paint.getColor();
    GPixel new_pix = colortoPixel(color); 
    int src_a = GPixel_GetA(new_pix);

    GBlendMode blendMode = paint.getBlendMode();

    if (paint.getShader() == nullptr && (src_a == 0 || src_a == 255)) {
        blendMode = optimizeMode(blendMode, src_a);
    }

    BlendFunc blendFunc = getBlendFunc(blendMode);

    // create array of line segments, where each element is a line segment,
    // represented by x1, y1, x2, y2
    std::vector<std::vector<float>> segments; 
    for (int i = 0; i < count; i++) {
        int lo = i, hi = i+1;
        if (hi == count) {
            hi = 0;
        }

        if (mappedPoints[lo].y < 0 && mappedPoints[hi].y < 0 || mappedPoints[lo].y >= fDevice.height() && mappedPoints[hi].y >= fDevice.height()) {
            continue;
        }

        int lower = lo, higher = hi;
        if (mappedPoints[lo].y > mappedPoints[hi].y) {
            lower = hi, higher = lo;
        }
        
        std::vector<float> seg{mappedPoints[lower].x, 
            mappedPoints[lower].y,
            mappedPoints[higher].x,
            mappedPoints[higher].y
        };

        segments.push_back(seg);
    }

    int validCount = segments.size();

    // sort segments by reverse order of lower y-coord
    std::sort(segments.begin(), segments.end(),
              [](const std::vector<float>& a, const std::vector<float>& b) {
                  return std::min(a[1], a[3]) < std::min(b[1], b[3]);
              });

    // record reciprocal slope for each segment
    std::vector<float> recipSlopes(validCount);
    for (int i = 0; i < validCount; i++) { 
        
        // designate segments that have a horizontal slope
        if (segments[i][3] == segments[i][1]) {
            recipSlopes[i] = std::numeric_limits<float>::max();
        } 
        else {
            float recipSlope = (segments[i][2] - segments[i][0]) / (segments[i][3] - segments[i][1]);
            recipSlopes[i] = recipSlope; 
        }
    }

    // iterate through line segments, keeping track of which two 
    // we are currently using to fill pixels
    int segA = 0;
    int segB = 1;
    while (segA < validCount && segB < validCount) {
        
        // skip over line segments with a flat slope
        if (recipSlopes[segA] == std::numeric_limits<float>::max()) {
            segA = std::max(segA, segB) + 1;
            continue;
        }
        if (recipSlopes[segB] == std::numeric_limits<float>::max()) {
            segB = std::max(segA, segB) + 1;
            continue;
        }

        // shoot x-rays and fill all pixels which hit these two lines
        int y_start = std::max(0, GRoundToInt(std::max(segments[segA][1], segments[segB][1])));
        int y_stop = std::min(fDevice.height(), GRoundToInt(std::min(segments[segA][3], segments[segB][3])));
        
        bool last_iter = false;
        if ((segA == validCount-1 || segB == validCount-1) && y_stop < fDevice.height()) {
            y_stop += 1;
            last_iter = true;
        }

        float A_x = getXFromY(y_start+0.5f, segments[segA][0], segments[segA][1], recipSlopes[segA]);
        float B_x = getXFromY(y_start+0.5f, segments[segB][0], segments[segB][1], recipSlopes[segB]);
        float delta_a = recipSlopes[segA];
        float delta_b = recipSlopes[segB];
        
        //y0                /
        //y0.25  s         / s
        //y0.5         p  /               /
        //y0.75  s       /   s
        //y1            /
        
        
        for (int y = y_start; y < y_stop; y++) {
            // entry point
            int x_start = std::max(0, GRoundToInt(std::min(A_x, B_x)));
            // exit point
            int x_end = std::min(fDevice.width(), GRoundToInt(std::max(A_x, B_x)));
            
            
            // where is the line? A_x is line coord at y=0.5 (mid)
            // float A_top_x = A_x - 0.25*delta_a; // at y=0.25
            // float A_bot_x = A_x + 0.25*delta_a; // at y=0.75
            // float B_top_x = B_x - 0.25*delta_b;
            // float B_bot_x = B_x + 0.25*delta_b;

            if (paint.getShader() == nullptr) {
                
                for (int x = x_start; x < x_end; x++) {
                    GPixel* prev_pix = fDevice.getAddr(x, y);
                        
                    int big_dst_a = GPixel_GetA(*prev_pix);
                    
                    GPixel blended_pix = blendFunc(new_pix, prev_pix, src_a, big_dst_a);
                    
                    // if (x == x_start || x == x_end - 1) {
                    //     //Test
                    //     int samplesTotal = 4;
                        
                    //     // the actual pixel is at (x, y+0.5)
                    //     float x_left = x - 0.25;
                    //     float x_right = x + 0.25;
                        
                    //     // the 4 samples are at (x_left, y_top), etc
                    //     // we don't care about calculating y_top, just comparing x to A_top_x and B_top_x
                    //     int samplesIn = 0;
                    //     float min_top_x = std::min(A_top_x, B_top_x);
                    //     float max_top_x = std::max(A_top_x, B_top_x);
                    //     float min_bot_x = std::min(A_bot_x, B_bot_x);
                    //     float max_bot_x = std::max(A_bot_x, B_bot_x);
                    //     // check top left sample:
                    //     if (min_top_x <= x_left && x_left <= max_top_x) samplesIn++;
                    //     // top right:
                    //     if (min_top_x <= x_right && x_right <= max_top_x) samplesIn++;
                    //     if (min_bot_x <= x_left && x_left <= max_bot_x) samplesIn++;
                    //     if (min_bot_x <= x_right && x_right <= max_bot_x) samplesIn++;
                        
                    //     float samplingAlpha = samplesIn / samplesTotal;
                        
                    //     // blended_pix = alpha*blended_pix + (1-a)*prev_pix;
                    //     // dst_a is ignored in srcOver
                    //     blended_pix = srcOver(blended_pix, prev_pix, (int)255*samplingAlpha, 255);
                    // }

                    *prev_pix = blended_pix;
                }
            }
            else {
                GShader* shader = paint.getShader();
                
                if (shader->isOpaque()) {
                    blendMode = optimizeMode(blendMode, 255);
                }
                BlendFunc blendFunc = getBlendFunc(blendMode);

                shader->setContext(ctm); 

                GPixel storage[x_end - x_start];
                shader->shadeRow(x_start, y, x_end - x_start, storage);
                
                for (int x = x_start; x < x_end; x++) {
                    GPixel* prev_pix = fDevice.getAddr(x, y);
                        
                    int big_dst_a = GPixel_GetA(*prev_pix);
                    int big_src_a = GPixel_GetA(storage[x - x_start]);
                    GPixel blended_pix = blendFunc(storage[x - x_start], prev_pix, big_src_a, big_dst_a);
                    // antialias
                    
                    
                    *prev_pix = blended_pix;
                }
            }

            A_x += delta_a;
            B_x += delta_b;
        }

        if (y_stop == fDevice.height() || last_iter) {
            break;
        }
        // increment one or both of the line segments if their higher y-coord
        // equals the y coord at the end of the interval
        if (GRoundToInt(segments[segA][3]) == y_stop) {
            segA = std::max(segA, segB) + 1;
        }
        if (GRoundToInt(segments[segB][3]) == y_stop) {
            segB = std::max(segA, segB) + 1;
        }
    }

    return;
}


void MyCanvas::drawRect(const GRect& rect, const GPaint& paint) {
    GPoint rect_points[4];
    rect_points[0] = {rect.left, rect.top};
    rect_points[1] = {rect.right, rect.top};
    rect_points[2] = {rect.right, rect.bottom};
    rect_points[3] = {rect.left, rect.bottom};

    MyCanvas::drawConvexPolygon(rect_points, 4, paint);
    return;
}


std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    return std::unique_ptr<GCanvas>(new MyCanvas(device));
}


std::string GDrawSomething(GCanvas* canvas, GISize dim) {
    GColor red {1.0, 0.0, 0.0, 1.0};
    canvas->clear(red);
    
    GColor green {0.0, 1.0, 0.0, 1.0};
    GPaint green_paint(green);

    //GPoint rectPoints[] = {{10, 10}, {10, 20}, {20, 20}, {20, 10}};
    GPath path;
    path = GPath();
    path.addCircle({100, 100}, 40, GPath::kCCW_Direction);
    //canvas->drawPath(path, green_paint); 
    
    // template <typename T> struct GSize {
    //     T width, height;
    // };
    // typedef GSize<int> GISize;
    
    // draw a letter D with a path
    // .---c
    // |
    // |        
    // |
    // .   c 
    //GPath path1;
    path = GPath();
    GPoint top = {dim.width * 0.1f, dim.height * 0.1f};
    GPoint bottom = {dim.width * 0.1f, dim.height * 0.8f};
    GPoint c1 = bottom + GPoint{dim.width * 0.3, 0.0f};
    GPoint c2 = GPoint{dim.width * 0.3, 0.0f} + top;
    path.moveTo(top);
    path.lineTo(bottom);
    path.cubicTo(c1, c2, top);
    // canvas->drawPath(path, green_paint);

    // draw a letter v
    //
    //      .            .
    //    .                .
    //
    //             
    //             .
    float topRow = dim.height * 0.2;
    float midRow = dim.height * 0.3;
    float bottomRow = dim.height * 0.7;

    const GPoint verts[] = {
        {dim.width * 0.1f, midRow},
        {dim.width * 0.15f, topRow},
        {dim.width * 0.9f, midRow}, 
        {dim.width * 0.85f, topRow},
        {dim.width * 0.5f, bottomRow},
    };
    const GColor colors[] = {
        {0, 0, 1, 1},
        {0, 0.2, 0.8, 1},
        {1, 1, 0, 1},
        {1, 0.8, 0, 1},
        {1, 1, 1, 0.4},
    };
    // const GPoint texs[],
    // int count;
    const int indices[] = {0, 1, 4, 2, 3, 4};
    const GPaint paint;
    canvas->drawMesh (verts, colors, nullptr, 2, indices, paint);

    return "storm king";
}
