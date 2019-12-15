/*
 * File: fractals.cpp
 * --------------------------
 * Name:
 * Section leader:
 * This file contains fractal problems for CS106B.
 */

#include "fractals.h"
#include <cmath> //For floor function (overload for doubles)
#include <algorithm> //for max function
#include "random.h"

using namespace std;

const int LEAF_COLOR = 0x2e8b57;   /* Color of all leaves of recursive tree (level 1) */
const int BRANCH_COLOR = 0x8b7765; /* Color of all branches of recursive tree (level >=2) */
const double PI = 3.14159265;
const double LOGTWO = log(2.0);

////NOTE: monotone cubic interpolate makes for smoother gradient
///**
// * @brief linearInterpolate
// * @param y0
// * @param y1
// * @param x
// * @return
// */
//int linearInterpolate(int y0, int y1, double x) {
//    return round(y0 + x * (y1 - y0));
//}

GPoint newTestPoint(GPoint& testPoint, GPoint ctrlPoint) {
    double newX = testPoint.getX() + 0.5*(ctrlPoint.getX()-testPoint.getX());
    double newY = testPoint.getY() + 0.5*(ctrlPoint.getY()-testPoint.getY());
    return GPoint(newX, newY);
}

/**
 * IMPLEMENTATION NOTES ~ drawSierpinskiTriangle ~
 *
 * Conducts error checking for invalid values for size, x, y, and order. Checks for null case
 * by doing nothing. Follows a reursive strategy by implementing the base case to draw an order
 * one Sierpinski triangle (unfilled) and then implements the recursive step by drawing
 * three triangles of the previous generation in the correct orientation.
 *
 * NOTE: ODDLY BUGGY: WILL WORK SOMETIMES BUT NOT OTHERS. SYSTEM TROUBLES? ISSUE
 * WITH INITIAL POINTS? WHO KNOWS!
 *
 */
void drawSierpinskiTriangle(GWindow& gw, double x, double y, double size, int order, int color) {
    GPoint p0(x, y);
    GPoint p1(x+size/2, y+size*cos(PI/6.0));
    GPoint p2(x+size, y);
    Vector<GPoint> controlPoints;
    controlPoints.add(p0);
    controlPoints.add(p1);
    controlPoints.add(p2);
    //Test for initial point in the middle of a non-Sierpinski point
    GPoint testPoint(x+size*0.5, y+size*cos(PI/6.0)*0.5);
    //GPoint testPoint(randomReal(x, x+size), randomReal(y, y+size*cos(PI/6.0)));
    gw.setColor(color);
    gw.drawLine(p0, p1);
    gw.drawLine(p1, p2);
    gw.drawLine(p2, p0);
    double slopeP1P0 = (p1.getY()-p0.getY()) / (p1.getX()-p0.getX());
    double slopeP2P1 = (p2.getY()-p1.getY()) / (p2.getX()-p1.getX());
    int i = 0;
    while (i < 15000) {
        if (testPoint.getX() < (x+size / 2)) {
            if ((testPoint.getY()-y) < (slopeP1P0*(testPoint.getX()-x))) {
                gw.drawPixel(testPoint.getX(), testPoint.getY());
                int indexOfRandomPoint = randomInteger(0, controlPoints.size()-1);
                testPoint = newTestPoint(testPoint, controlPoints[indexOfRandomPoint]);
                ++i;
            }
        } else {
            if ((testPoint.getY()-y) < (slopeP2P1*(testPoint.getX()-x)+2*size*cos(PI/6.0))) {
                gw.drawPixel(testPoint.getX(), testPoint.getY());
                int indexOfRandomPoint = randomInteger(0, controlPoints.size()-1);
                testPoint = newTestPoint(testPoint, controlPoints[indexOfRandomPoint]);
                ++i;
            }
        }
    }
    // TODO: Implement different colors for different layers
    // int N = 2048;
    // bool smoothColors = false
    //Vector<int> palette = setPalette(smoothColors, N);
//    if (size < 0 || x < 0 || y < 0 || order < 0) cerr << "insert error message here" << endl;
//    if (order == 0) {} //Do nothing if the order is zero
//    else if (order == 1) //Base case for drawing 1st order triangle
//    {
////        if (color != 0) {
////            gw.setColor(color); //not utilizing the color palette: one color for every iteration
////        }
////        else {
////            //TODO: add numIterations parameter which is passed by reference everywhere
////            // s.t. for each generation we can update its value
////            gw.setColor(palette[numIterations % palette.size()]);
////        }
//        gw.setColor(color);
//        gw.drawLine(x, y, x + size / 2, y + (size * sin(PI / 3)));
//        gw.drawLine(x + size / 2, y + (size * sin(PI / 3)), x + size, y);
//        gw.drawLine(x + size, y, x, y);
//    }
//    else
//    {
//        //TODO: add setColor function and use color button to change color of lines drawm
//        //s.t. each layer of the recursion is drawn with a specific color
//        //make three calls of the next lowest order, recursively
//        drawSierpinskiTriangle(gw, x, y, size / 2, order - 1, color);
//        drawSierpinskiTriangle(gw, x + size / 4, y + (size * sin(PI / 3)) / 2, size / 2, order - 1, color);
//        drawSierpinskiTriangle(gw, x + size / 2, y, size / 2, order - 1, color);
//    }

}

/**
 * IMPLEMENTATION NOTES: ~ drawSquare ~
 *
 * Rather self-explanatory.
 *
 */
void drawSquare(GWindow& gw, double x, double y, double size) {
    gw.fillRect(x, y, size, size);
}

/**
 * IMPLEMENTATION NOTES: ~ drawMengerSponge ~
 *
 * Works very similarly to drawSierpinskiTriangle: utilizes a recursive strategy
 * and the recursive step draws nine squares of the previous generation.
 *
 */
void drawMengerSponge(GWindow& gw, double x, double y, double size, int order, int color) {
    if (size < 0 || x < 0 || y < 0 || order < 0) cerr << "insert error message here" << endl;
    if (order == 0) {} //Do nothing if the order is zero
    else if (order == 1) { //The base case
        gw.setColor(color); //Trial color with black
        drawSquare(gw, x, y, size);
    } else { //The recursive step: Must draw 8 squares along border overall
             //square by dividing into thirds
        drawMengerSponge(gw, x, y, size / 3, order - 1, color);
        drawMengerSponge(gw, x + size / 3, y, size / 3, order - 1, color);
        drawMengerSponge(gw, x + 2 * size / 3, y, size / 3, order - 1, color);
        drawMengerSponge(gw, x + 2 * size / 3, y + size / 3, size / 3, order - 1, color);
        drawMengerSponge(gw, x + 2 * size / 3, y + 2 * size / 3, size / 3, order - 1, color);
        drawMengerSponge(gw, x + size / 3, y + 2 * size / 3, size / 3, order - 1, color);
        drawMengerSponge(gw, x, y + 2 * size / 3, size / 3, order - 1, color);
        drawMengerSponge(gw, x, y + size / 3, size / 3, order - 1, color);
    }
}

/**
 * IMPLEMENTATION NOTES: ~ drawTree ~
 *
 * This version of the drawTree function does the recursion. In this case the
 * "base case" is of order 0, so we do all the drawing in the recursive step essentially.
 * Namely, each generation requires a "base" trunk, and the 7 child branches are drawn
 * at the tip of the parent trunk for this generation.
 *
 */
void drawTree(GWindow& gw, double x, double y, double size, int order, double angle)
{
    if (x < 0 || y < 0 || order < 0 || size < 0) cerr << "Insert error message here" << endl;
    else if (order != 0)
    {
        if (order == 1)
        {
            gw.setColor(LEAF_COLOR);
        }
        else
        {
            gw.setColor(BRANCH_COLOR);
        }
        gw.drawPolarLine(x, y, size / 2, angle);
        for (int theta = -45; theta <= 45; theta += 15)
        {
            drawTree(gw, x + cos(PI * angle / 180) * size / 2, y - sin(PI * angle / 180) * size / 2, size / 2, order - 1, angle + theta);
        }
    }
}

/**
 * IMPLEMENTATION NOTES: ~ drawTree ~
 *
 * This is simply a wrapper function which defines where the base of each trunk
 * should be drawn as well as the angle of base branch / trunk.
 *
 */
void drawTree(GWindow& gw, double x, double y, double size, int order) {
    drawTree(gw, x + size / 2, y + size, size, order, 90);
}

/**
 * IMPLEMENTATION NOTES ~ mandleBrotSet ~
 *
 * Creates a GBufferedImage object which is added to the graphics window and is
 * converted to a grid for ease of access to individual pixel "elements." Pixel
 * drawing is done via two loops and within inner loop there exist two cases:
 * the argument value for color is not black, meaning we solid fill the MandelBrot
 * Set with the color argument and leave all pixels in the complement uncolored
 * (in this case they appear white). The second case involves explicit use of the
 * color palette. All pixels determined to be inside the Mandelbrot set are colored
 * black and the number of iterations for all other pixels is used to determine
 * the color of that pixel. A final if else statement is used to differentiate whether
 * pixels should be drawn with continuous or discrete color gradients / banding.
 *
 */
void mandelbrotSet(GWindow& gw, double minX, double incX,
                   double minY, double incY, int maxIterations, int color, bool smoothGradient) {

    // Creates palette of colors
    // To use palette:
    // pixels[r][c] = palette[numIterations % palette.size()];
    int N = 2048;
    //bool smoothColors = true;
    //Initialize and populate the color palette
    Vector<int> palette = setPalette(smoothGradient, N);

    int width = gw.getCanvasWidth();
    int height = gw.getCanvasHeight();
    GBufferedImage image(width,height,0xffffff);
    gw.add(&image);
    Grid<int> pixels = image.toGrid(); // Convert image to grid

    for (int i = 0; i < pixels.numRows(); i++)
    {
        for (int j = 0; j < pixels.numCols(); j++)
        {
            //implementing this function to work with remaining iterations and not
            //number of iterations makes for harder reading...
            //NOTE: maxIterations - remainingIterations == numIterations
            //z is c within recursive def z_n+1 = z_n^2 + c
            Complex z = Complex(minX + j * incX, minY + i * incY), zn(0.0, 0.0);
            int remainingIterations = mandelbrotSetIterations(z, zn, maxIterations);

            if (color != 0) //Not using palette, only color pixels in set
            {
                if (remainingIterations == 0)
                {
                    pixels[i][j] = color;
                }
            } //TODO: Serious debugging on color gradient code
            else //utilizing palette: obtain "contours" of same colour by colouring
                 //according to how many iterations are needed to determine whether
                 //point in set or not (or more generally, regardless whether the point
                 //is in the set or not)
            {
                // Color points within the set black
                if (remainingIterations == 0)
                {
                    pixels[i][j] = color;
                }
                else
                {
                    if (smoothGradient) {
                        //Implementation of smooth color gradient
                        int numIterations = maxIterations - remainingIterations;
                        int pixelColor = gradientColorMapping(numIterations, zn, palette);
                        pixels[i][j] = pixelColor;
                    } else {
                        //Below is implementation of "plotting" of pixels w/o continuous gradient
                        pixels[i][j] = palette[(maxIterations - remainingIterations) % palette.size()];
                    }

                }
            }
        }
    }

    image.fromGrid(pixels); // Converts and puts the grid back into the image
}

/**
 * IMPLEMENTATION NOTES ~ mandelBrotSetIterations ~
 *
 * Serves as a wrapper function to the function which does the recursion. This
 * function merely initializes z -> (0,0) to initiate the recursion.
 *
 */
int mandelbrotSetIterations(Complex c, Complex& zn, int maxIterations) {
    return mandelbrotSetIterations(Complex(0, 0), c, zn, maxIterations);
}

//int mandelbrotSetIterations(Complex z, Complex c, int remainingIterations, int &maxIterations) {
//    if ((z * z + c).abs() >= 4 || remainingIterations == 0)
//    {
//        return maxIterations - remainingIterations;
//    }
//    else
//    {
//        return mandelbrotSetIterations(z * z + c, c, remainingIterations - 1, maxIterations);
//    }
//}

/**
 * IMPLEMENTATION NOTES: ~ mandelBrotSetIterations ~
 *
 * This function actually conducts the recursion. It sets zn to the current
 * value of z under the recursive definition for the mandelbrot set. It then sets up
 * a recursive strategy by dividing control flow between cases in which we've either
 * diverged or determined that the point c is in the mandelbrot set, or continue the
 * recursion if neither of these conditions have been met. The parameter zn is called
 * by reference s.t. the caller has access to the changed value.
 *
 */
int mandelbrotSetIterations(Complex z, Complex c, Complex& zn, int remainingIterations) {
    zn = z * z + c;
    if ((zn).abs() >= 4 || remainingIterations == 0)
    {
        return remainingIterations;
    }
    else
    {
        return mandelbrotSetIterations(zn, c, zn, remainingIterations - 1);
    }
}

/**
 * IMPLEMENTATION NOTES ~ initMap ~
 *
 * I mean, this one's pretty intuitive chief...
 *
 */
Map<double, int> initMap() {
    Map<double, int> gradientMap;
    gradientMap.put(0.0, 1892);
    gradientMap.put(0.16, 2124747);
    gradientMap.put(0.42, 15597567);
    gradientMap.put(0.6425, 16755200);
    gradientMap.put(0.8575, 512);
    return gradientMap;
}

/**
 * IMPLEMENTATION NOTES: ~ gradientColorMapping ~
 *
 * Implements the normalized iteration count algorithmn, which utilizes the
 * escape time algorithm to produce the number of iterations and diverged
 * value of z to compute the color of the pixel at that z, creating a smooth
 * gradient of colors for all elements in the complement of the Mandelbrot set.
 * As of 12/20/2018, this algorithm needs work.
 *
 */
int gradientColorMapping(int numIterations, Complex z, Vector<int> colorPalette) {
    //if (numIterations <= maxIterations) { // Have a precheck for this so somewhat redundant?
        double logZn = log(z.abs());
        double nu = log(logZn / LOGTWO) / LOGTWO;
        numIterations = (int)(sqrt(numIterations + 1 - nu) * 256) % colorPalette.size();
        return colorPalette[numIterations];
    //}
}

/**
 * IMPLEMENTATION NOTES: ~ linspace ~
 *
 * This algorithm produces the step value h based off of input parameters a, b, and N.
 * It then assigns a to a double val which is placed into the vector xs at position i
 * and is updated within the loop until all values in xs are evenly spaced between a, and b.
 *
 * EFFICIENCY NOTE: Could probably leave xs uninitialized (remove an O(N) operation) and
 * instead conduct N O(1) operations by calling member method add. s.t. we cut the
 * number of operations conducted by a factor of two
 *
 */
Vector<double> linspace(double a, double b, int N) {
    double h = (b - a) / (N - 1);
    Vector<double> xs;
    double val = a;
    for (int i = 0; i < N; ++i) {
        xs.add(val);
        val += h;
    }
    return xs;
}

/**
 * IMPLEMENTATION NOTES ~ monotonicCubicInterpolate ~
 *
 * Honestly this one's a doosey. Essentially, the algorithm initializes the map
 * if it is empty (which occurs on the first and only the first call to this
 * function) and then places the keys in the xs vector and the values in the ys
 * vector. From here, it calculates the slopes at each point, the c0, c1, and c2
 * coefficients for the cubic polynomial interpolation, and then determines which
 * onterval d is in and returns the result of f(d) where f is the interpolating
 * cubic polynomial.
 *
 */
double monotonicCubicInterpolate(double d, Map<double,int>& inMap) {
    // Precondition check
    if (inMap.isEmpty()) {
        inMap = initMap();
    }
    // Define a tolerance for comparing doubles
    double tol = 1.0e-8;
    Vector<double> xs, ys, dxs, dys, ms;

    // Loading keys and values in inMap into xs and ys vectors for ordered access..
    // TODO: Create a more efficient implementation
    for(double key : inMap) {
        xs.push_back(key);
        ys.push_back(inMap.get(key)); //Implicit typecast from int -> double
    }

    // Get consecutive differences and slopes
    for(int i = 0; i < xs.size() - 1; i++) {
        double dx = xs[i + 1] - xs[i];
        double dy = ys[i + 1] - ys[i];
        dxs.push_back(dx);
        dys.push_back(dy);
        ms.push_back(dy / dx);
    }

    // Get degree - 1 coeff
    Vector<double> c1s;
    // Edge case: obtain first element
    c1s.push_back(ms[0]);
    for (int i = 0; i < dxs.size() - 1; i++) {
        double m = ms[i], mNext = ms[i + 1];
        if (m * mNext <= 0) {
            c1s.push_back(0.0);
        } else {
            double dx_ = dxs[i];
            double dxNext = dxs[i + 1];
            double common = dx_ + dxNext;
            c1s.push_back(3 * common / ((common + dxNext) / m + (common + dx_) / mNext));
        }
    }
    // edge case: obtain last element
    c1s.push_back(ms[ms.size() - 1]);

    // Get degree - 2 and degree - 3 coeff
    Vector <double> c2s, c3s;
    for (int i = 0; i < c1s.size() - 1; i++) {
        double c1 = c1s[i], m_ = ms[i], invDx = 1 / dxs[i], common_ = c1 + c1s[i + 1] - m_ - m_;
        c2s.push_back((m_ - c1 - common_) * invDx);
        c3s.push_back(common_ * invDx * invDx);
    }

    // Now interpolate at d and return f(d)
    if (abs(d - xs[xs.size() - 1]) < tol) return ys[xs.size() - 1];

    int low = 0, mid, high = c3s.size() - 1;
    while (low <= high) {
        mid = floor(0.5 * (low + high)); //implicit downcast? double -> int from floor return??
        double xHere = xs[mid];
        if (xHere < d) low = mid + 1;
        else if (xHere > d) high = mid - 1;
        else if (abs(xHere - d) < tol) return ys[mid];
    }

    int i = max(0, high); //where is max function defined??
    double diff = d - xs[i], diffSq = diff * diff;
    return ys[i] + c1s[i] * diff + c2s[i] * diffSq + c3s[i] * diff * diffSq;
}

/**
 * IMPLEMENTATION NOTES ~ setPalette ~
 *
 * Contains two control flow cases: either we are drawing a smooth gradient or
 * discrete bands of colors. For the continuous gradient, we assign a vector
 * xs the result of linspace using the correct interval [a,b] and input parameter
 * N, create a map object (should probably change whether that map is a constant and
 * whether we should return by reference...but all that is slightly beyond what
 * you've mastered so far) and then generates the interpolated values of the xs,
 * placing them in a vector of ints representing the interpolated colors of the
 * gradient. The other control route is the "stock" setPalette function, which
 * has a designated color palette hardcoded as a string object and then that
 * string is parsed into the individual colors which are placed in the colors vector
 * (after properly converting the string representation of those colors to rgb values)
 * and that vector is returned.
 *
 */
Vector<int> setPalette(bool smoothColors, int N) {
    //allocate space in memory for a vector of integers representing colors
    Vector<int> colors;
    // Follow this control path if we want "continuous" colors
    if (smoothColors) {
//        int N = 2048;
        Vector<double> xs;
        // TODO create implementation that is independent of the largest key in gradientMap
        xs = linspace(0.0, 0.8575, N);
        Map<double, int> gradMap;

        for (int i = 0; i < N; ++i) {
            colors.add((int)monotonicCubicInterpolate(xs[i], gradMap));
        }
    }
    // Otherwise setPalette is simply the unmodified version
    else {
        //https://www.colourlovers.com/palette/848743/(_%E2%80%9D_)
        //string colorSt = "#490A3D,#BD1550,#E97F02,#F8CA00,#8A9B0F,#490A3D,#BD1550,#E97F02";
        // Feel free to replace with any palette.
        // You can find palettes at:
        // http://www.colourlovers.com/palettes

        //My custom gradient 1
        //string colorSt = "#903E07,#960B0B,#0F0B1B,#030927,#050F42,#071559,#06166A,#152CA4,#1330C6,#4040F4,#8686E2,#F7E0D6,#F5C16E,#F1AD1A,#C78506,#8C6B06";

        //My custom gradient 2
        string colorSt = "#421E0F,#19071A,#09012F,#040449,#000764,#0C2C8A,#1852B1,#397DD1,#86B5E5,#D3ECF8,#F1E9BF,#F8C95F,#FFAA00,#CC8000,#995700,#6A3403";

        //My custom gradient 3
        //string colorSt = "#000764,#206BCB,#EDFFFF,#FFAA00,#000200";
        // Example palettes:
        // http://www.colourlovers.com/palette/4480793/in_the_middle
        //string colorSt = "#A0B965,#908F84,#BF3C43,#9D8E70,#C9BE91,#A0B965,#908F84,#BF3C43";

        // http://www.colourlovers.com/palette/4480786/Classy_Glass
        // string colorSt = "#9AB0E9,#C47624,#25269A,#B72202,#00002E,#9AB0E9,#C47624,#25269A";

        // The following is the "Hope" palette:
        // http://www.colourlovers.com/palette/524048/Hope
        //string colorSt =  "#04182B,#5A8C8C,#F2D99D,#738585,#AB1111,#04182B,#5A8C8C,#F2D99D";
        Vector<string>colorsStrVec = stringSplit(colorSt,",");
        for (string color : colorsStrVec) {
            colors.add(convertColorToRGB(trim(color)));
            //cout << convertColorToRGB(trim(color)) << endl;
        }
    }
    //Either way, want to return palette of colors
    return colors;
}
