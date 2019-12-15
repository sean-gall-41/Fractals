/*
 * CS 106B/X Fractal Recursion Problems
 * This file declares necessary function prototypes, so that other files
 * are able to call the recursive functions you will write.
 * 
 * !!! DO NOT MODIFY THIS FILE !!!
 * !!! DO NOT MODIFY THIS FILE !!!
 * !!! DO NOT MODIFY THIS FILE !!!
 *
 * Please do not modify this provided file. Your turned-in files should work
 * with an unmodified version of all provided code files.
 *
 * (If you want to declare your own helper function prototypes,
 *  do that near the top of your .cpp file, not here.)
 */

#ifndef _recursionproblems_h
#define _recursionproblems_h

#include <iostream>
#include <string>
#include "gbufferedimage.h"
#include "gwindow.h"
#include "vector.h"
#include "map.h"
#include "complex.h"

/**
 * Creates a map object from doubles in the unit interval to
 * integers which encode the rgb value of each of five control points
 * for a custom gradient
 * @usage Map<double, int> myMap = initMap();
 * @return returns a map from doubles to ints
 */
Map<double, int> initMap();

/**
 * Uses the current iteration number to map a real number in the interval [0,1]
 * to an appropriate value within the gradient
 *
 * @usage int myColor = gradientColorMapping(n, z, palette);
 * @param n
 * @param maxIterations
 * @return
 */
int gradientColorMapping(int numIterations, Complex z, Vector<int> colorPalette);

/**
 * a C++ equivalent of numpy's linspace
 *
 * @usage Vector<double> xs = linspace(0.0, 1.0, int N = 2048);
 * @param a - first value in interval [a, b]
 * @param b - second value in interval [a, b]
 * @param N - the number of evenly spaced values to generate within [a, b]
 * @return a vector of N evenly spaced values in the interval [a, b]
 */
Vector<double> linspace(double a, double b, int N);

/**
 * generates an interpolated value f(d) using points in inMap as the reference points
 * Will be used to generate a vector of interpolated points
 *
 * @usage double y = monotonicCubicInterpolate(x, colorMapping);
 * @param inMap - the generated map of control points for the gradient
 * @return a vector object of say, 2048 ints generated from cubic interpolation
 *         of control points
 */
double monotonicCubicInterpolate(double d, Map<double,int>& inMap);

/**
 * Draws a Sierpinski triangle of the specified size and order, placing its
 * top-left corner at position (x, y).
 *
 * This will be called by fractalgui.cpp.
 *
 * @param gw - The window in which to draw the Sierpinski triangle.
 * @param x - The x coordinate of the top-left corner of the triangle.
 * @param y - The y coordinate of the top-left corner of the triangle.
 * @param size - The length of one side of the triangle.
 * @param order - The order of the fractal.
 */
void drawSierpinskiTriangle(GWindow& window, double x, double y, double size, int order, int color);

/**
 * A helper function which draws a single square with the specified size side
 * length, placing its top-left corner at position (x, y).
 *
 * @param gw - The window in which to draw the square.
 * @param x - The x coordinate of the top-left corner of the square.
 * @param y - The y coordinate of the top-left corner of the square.
 */
void drawSquare(GWindow& gw, double x, double y, double size);

/**
 * Draws a Menger sponge of the specified size and order, placing its
 * top-left corner at position (x, y).
 *
 * This will be caled fractalgui.cpp.
 *
 * @param gw - The window in which to draw the Menger Sponge.
 * @param x - The x coordinate of the top-left corner of the square.
 * @param y - The y coordinate of the top-left corner of the square.
 * @param size - The length of one side of the square.
 * @param order - The order of the fractal.
 */
void drawMengerSponge(GWindow& window, double x, double y, double size, int order, int color);

/**
 * Draws a recursive tree fractal of the specified size and order, placing the
 * bounding box's top-left corner at position (x, y).
 *
 * Used within wrapper function of the same name, with signature only different
 * in lack of angle parameter
 *
 * @param gw - The window in which to draw the Sierpinski triangle.
 * @param x - The x coordinate of the top-left corner of the triangle.
 * @param y - The y coordinate of the top-left corner of the triangle.
 * @param size - The length of one side of the triangle.
 * @param order - The order of the fractal.
 * @param angle - the angle which each branch is drawn at w.r.t. the x-axis
 */
void drawTree(GWindow& gw, double x, double y, double size, int order, double angle);

/**
 * Serves as a wrapper function for abover drawTree method which sets the angle object
 * initially to 90 degrees as well as the initial coordinates (x, y) for each "level"
 * of the recursion
 *
 * @param gw - The window in which to draw the recursive tree.
 * @param x - The x coordinate of the top-left corner of the bounding box.
 * @param y - The y coordinate of the top-left corner of the bounding box.
 * @param size - The length of one side of the bounding box.
 * @param order - The order of the fractal.
 */
void drawTree(GWindow& window, double x, double y, double size, int order);

/**
 * Draws a Mandelbrot Set in the graphical window give, with maxIterations
 * (size in GUI) and in a given color (zero for palette)
 *
 * This will be called by fractalgui.cpp.
 *
 * @param gw - The window in which to draw the Mandelbrot set.
 * @param minX - left-most column of grid
 * @param incX - increment value of columns of grid
 * @param minY - top-most row of grid
 * @param incY - increment value of rows of grid
 * @param maxIterations - The maximum number of iterations to run recursive step
 * @param color - The color of the fractal; zero if palette is to be used
 */
void mandelbrotSet(GWindow& window, double minX, double incX,
                   double minY, double incY, int maxIterations, int inColor, bool smoothGradient);

/**
 * Serves as a wrapper function for the overloaded version of this function (which
 * actually does the recursion) This version initializes z -> (0,0) for the recursion
 *
 * @param c - the complex value we would like to determine is in the set or not
 * @param zn - the current complex value in the recursive formula z_(n+1) := z^2 + c
 * @param maxIterations - the maximum number of iterations to run recursion until
 *        "bailout" thus determining this c as belonging to the mandelbrot set
 * @return the number of iterations necessary to reach divergence or bailout
 */
int mandelbrotSetIterations(Complex c, Complex& zn, int maxIterations);

/**
 * Runs the Mandelbrot Set recursive formula on complex number c a maximum
 * of maxIterations times.
 *
 * This will be called by you. Think about how this fits with the other two functions.
 *
 * @param z - Complex number serving as zn of this generation
 * @param c - Complex number to use for recursive formula - the complex coordinate at which
 *        this pixel shall be drawn in mandelBrot set
 * @param zn - the current value of z in the recursive formula
 * @param maxIterations - The maximum number of iterations to run recursive step
 * @return number of iterations needed to determine if c is unbounded (i.e. belongs
 *         to Mandelbrot set or not)
 */
int mandelbrotSetIterations(Complex z, Complex c, Complex& zn, int remainingIterations);

//int floodFill(GWindow& window, int x, int y, int color);

/**
 * initializes the palette used to color all points in the complement of the Mandelbrot set
 *
 * @usage double myPalette = setPalette(true, int N = 2048)
 * @param smoothColors - bool value used to determine if drawing smooth gradient
 *        of colors or not
 * @param N - the desired length of the smoothed palette
 * @return the palette of colors used to color elements in the complement of the
 *         Mandelbrot set
 */
Vector<int> setPalette(bool smoothColors, int N);

#endif

/*
 * !!! DO NOT MODIFY THIS FILE !!!
 * !!! DO NOT MODIFY THIS FILE !!!
 * !!! DO NOT MODIFY THIS FILE !!!
 */
