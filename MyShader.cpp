#include "MyCanvas.h"

bool MyShader::isOpaque() {
    return bm.isOpaque();
}


bool MyShader::setContext(const GMatrix& ctm) {
    std::optional<GMatrix> ctmInverse = ctm.invert();
    std::optional<GMatrix> lmInverse = lMatrix.invert();

    if (!(ctmInverse && lmInverse)) {return false;}

    // (ctm * lm)^-1 = lm^-1 * ctm^-1
    inverse = GMatrix::Concat(lmInverse.value(), ctmInverse.value());
    return true;
}
  

void MyShader::shadeRow(int x, int y, int count, GPixel row[]) {
    GPoint pt = {x + 0.5f, y + 0.5f};
    inverse.GMatrix::mapPoints(&pt, 1);
    float xi = pt.x;
    float yi = pt.y;
    float delta_x = inverse[0];
    float delta_y = inverse[1];
    
    float invBmWidth = 1.0f / (bm.width());
    float invBmHeight = 1.0f / (bm.height());
    
    for (int i = 0; i < count; i++) {
        // floor
        int xi_floored = GFloorToInt(xi);
        int yi_floored = GFloorToInt(yi);
        
        switch (tileMode) {
            case GTileMode::kClamp:
                // we will clamp later
                break;
            case GTileMode::kRepeat: {
                float xi_frac = xi_floored * invBmWidth;
                // xi_frac will be [0-1)
                xi_frac = xi_frac - floorf(xi_frac);
                // xi_frac * width will be [0-width) so we need to clamp later to avoid width
                xi_floored = GRoundToInt(xi_frac * bm.width());
                
                float yi_frac = yi_floored * invBmHeight;
                yi_frac = yi_frac - floorf(yi_frac);
                yi_floored = GRoundToInt(yi_frac * bm.height());
                break;
            }
            case GTileMode::kMirror: {
                // x-2*floor(x/2)
                float xi_frac = xi_floored * invBmWidth * 0.5f;
                // xi_frac will be [0-1)
                xi_frac = 2*(xi_frac - floorf(xi_frac));
                //1 - |x - 1|
                xi_frac = 1 - fabs(xi_frac - 1);
                // xi_frac * width will be [0-width) so we need to clamp later to avoid width
                xi_floored = GRoundToInt(xi_frac * bm.width());
                
                float yi_frac = yi_floored * invBmHeight * 0.5f;
                yi_frac = 2*(yi_frac - floorf(yi_frac));
                yi_frac = 1 - fabs(yi_frac - 1);
                yi_floored = GRoundToInt(yi_frac * bm.height());
                break;
            }
        }
        
        // clamp no matter which mode for safety with getAddr
        xi_floored = std::max(0, xi_floored);
        yi_floored = std::max(0, yi_floored);
        xi_floored = std::min(bm.width()-1, xi_floored);
        yi_floored = std::min(bm.height()-1, yi_floored);
        
        // get pixel from bitmap
        row[i] = *(bm.getAddr(xi_floored, yi_floored));
        // update along basis vectors
        xi += delta_x;
        yi += delta_y;
    }
}


std::unique_ptr<GShader> GCreateBitmapShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode tileMode) {
    return std::unique_ptr<GShader>(new MyShader(bitmap, localMatrix, tileMode));
}
