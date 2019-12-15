// my project

#include <iostream>
#include <string>
#include "console.h"
#include "gwindow.h" // for GWindow
#include "simpio.h"  // for getLine
#include "vector.h"  // for Vector
#include "map.h"
#include <cmath>
#include <algorithm>
using namespace std;

//function prototypes

/**
 * @brief initMap - Creates a map object from doubles in the unit interval to
 *        integers which encode the rgb value of each of five control points
 *        for a custom gradient
 * @return returns a map from doubles to ints
 */
Map<double, int> initMap();

/**
 * @brief cubicInterpolate - generates a vector object storing the resulting
 *        gradoent from monotonic cubic interpolation of values from inMap
 * @param inMap - the generated map of control points for the gradient
 * @return a vector object of say, 2048 ints generated from cubic interpolation
 *         of control points
 */
double monotonicCubicInterpolate(double d, Map<double,int>& inMap);

Vector<int> setPalette(bool smoothColors, int N);

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

int main() {
    int N = 2048;
//    Vector<double> xis, yis(N, 0.0);
    Vector<double> xis = linspace(0.0, 0.8575, N);

////    for (int i = 0; i < 12; ++i) {
////        cout << xis[i] << " " << yis[i] << endl;
////    }
//    Map<double, int> gradMap;
//    for (int i = 0; i < N; ++i) {
//        yis[i] = monotonicCubicInterpolate(xis[i], gradMap);
//    }
    bool smoothColors = true;
    Vector<int> yis = setPalette(smoothColors, N);
    for (int i = 0; i < N; ++i) {
        cout << xis[i] << " " << yis[i] << endl;
    }

    string myStr = getLine("Press any button to quit...");
    return 0;
}

//Helper function to initialize the gradient map
Map<double, int> initMap() {
    Map<double, int> gradientMap;
    gradientMap.put(0.0, 1892);
    gradientMap.put(0.16, 2124747);
    gradientMap.put(0.42, 15597567);
    gradientMap.put(0.6425, 16755200);
    gradientMap.put(0.8575, 512);
    return gradientMap;
}

//A helper function which calculates a cubic interpolation of values
//representing a cubic interpolation of colors from inMap
// TODO: Return Vector<ints> and typecast element-wise to double as caller,
// or return Vector<double>, letting the callee typecast (which has already been done)
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

Vector<int> setPalette(bool smoothColors, int N) {
    Vector<int> colors;
    if (smoothColors) {
        //int N = 2048;
        Vector<double> xs;
        // TODO create implementation that is independent of the largest key in gradientMap
        xs = linspace(0.0, 0.8575, N);
        Map<double, int> gradMap;

        for (int i = 0; i < N; ++i) {
            // WARNING: Truncation warnings
            colors.add((int)monotonicCubicInterpolate(xs[i], gradMap));
        }
    }
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

