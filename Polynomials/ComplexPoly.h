#ifndef _COMPLEX_POLY_H
#define _COMPLEX_POLY_H 1
/************************************************************************
 *                                                                      *
 * This subclass of object implements a polynomial object.  The         *
 * coefficients of the polynomial are Complex variables.                *
 * The usual operations are defined for this class, +, x, -, /, +=, etc.*
 * This class should be useful for digital signal processing.           *
 *                                                                      *
 * File:ComplexPoly.h                                                   *
 *                                                                      *
 * Note: for p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]              *
 *       p(x) has order n, but the length of the array is n+1.          *
 *                                                                      *
 ************************************************************************/
#ifdef USE_IO_STREAM
#include <fstream.h>
#endif

class   Complex;                                    // Class prototype
class   DoublePoly;
#ifndef NULL
#define NULL    0
#endif
#ifndef LOGICAL
#define LOGICAL    char
#endif


#define MAX_SHIFT           1024                        // Maximum shift right or left
#define MAX_STRING_SIZE     2048                        // Polynomial string size

class ComplexPoly
{
private:
  int       _order;                                     // order of the polynomial
  int       _size;                                      // actual length of the array
  Complex   *_coefficients;                             // array containing the polynomial coefficients
  char      *_polyString;                               // E.g., "0.5*x^n + -.2*x^n-1 ..."

  Complex   *_polyRoots;                                // Complex roots of a polynomial
//
// Private methods:
//
  void updateMemory (int newOrder);                     // Adjust memory size
  double maxCoeffAbs ();
  void updateOrder();                                   // Check that x[order+1] != 0, else update

public:
//
// Public methods:
//
/*********************************
 * Constructors, destructors:    *
 *********************************/
  ComplexPoly(int pOrder=0, const Complex *coeff=NULL);         // Initializes a polynomial
  ComplexPoly(const ComplexPoly &poly);
  ComplexPoly(Complex x);
  ComplexPoly(double x);
  ~ComplexPoly();

/**********************
 * Set parameters:    *
 **********************/
  void assign(int pOrder, const Complex *coeff);                // Sets a new polynomial
  void setRoots(int pOrder, const Complex *roots);              // Sets a new polynomial from roots

/**********************
 * Get parameters:    *
 **********************/
  const Complex *getCoefficients() const  {return _coefficients;} // Return the polynomial coefficients
  int getOrder() const                    {return _order;}      // Return the order of the polynomial

/**********************
 * Calculating:       *
 **********************/
  void reverse();                                           // Reverse the order of coefficients of A
  Complex evaluateAtComplexPoint(Complex &point) const;     // Evaluate the A polynomial at a complex point
  Complex *roots();                                         // Evaluate and return the roots of A

/****************************
 * The following functions  *
 * print the polynomial     *
 ****************************/
  const char *getPolyString();

#ifdef USE_IO_STREAM
  void printPoly(ofstream &fout);
#endif

/****************************
 * The following operators  *
 * modify the polynomial    *
 * by a second polynomial   *
 ****************************/
  ComplexPoly&  operator =  (const ComplexPoly& y);
  ComplexPoly&  operator =  (const DoublePoly& y);
  ComplexPoly&  operator =  (Complex y);
  ComplexPoly&  operator =  (double y);
  ComplexPoly&  operator += (const ComplexPoly& y);
  ComplexPoly&  operator -= (const ComplexPoly& y);
  ComplexPoly&  operator *= (const ComplexPoly& y);
  ComplexPoly&  operator /= (const ComplexPoly& y);
  ComplexPoly&  operator %= (const ComplexPoly& y);
  ComplexPoly&  operator >>= (int numberShifts);
  ComplexPoly&  operator <<= (int numberShifts);

/********************************
 * The following operators      *
 * modify the polynomial        *
 * by a Complex.                *
 ********************************/
  ComplexPoly&  operator += (Complex y);
  ComplexPoly&  operator -= (Complex y);
  ComplexPoly&  operator *= (Complex y);
  ComplexPoly&  operator /= (Complex y);

/********************************
 * The following operators      *
 * modify the polynomial        *
 * by a double.                 *
 ********************************/
  ComplexPoly&  operator += (double y);
  ComplexPoly&  operator -= (double y);
  ComplexPoly&  operator *= (double y);
  ComplexPoly&  operator /= (double y);

};
//
// The following are nonmember functions taking two arguments
//
/****************************
 * The following operators  *
 * make comparisons between *
 * two polynomials.         *
 ****************************/
  LOGICAL  operator == (const ComplexPoly& x, const ComplexPoly& y);           // Checks if two polynomials are equal
  LOGICAL  operator != (const ComplexPoly& x, const ComplexPoly& y);           // Checks if two polynomials are not equal
  LOGICAL  operator >  (const ComplexPoly& x, const ComplexPoly& y);           // Checks if x > y
  LOGICAL  operator >= (const ComplexPoly& x, const ComplexPoly& y);           // Checks if x >= y
  LOGICAL  operator <  (const ComplexPoly& x, const ComplexPoly& y);           // Checks if x < y
  LOGICAL  operator <= (const ComplexPoly& x, const ComplexPoly& y);           // Checks if x <= y

/****************************
 * The following operators  *
 * take 1 or 2 polynomials  *
 * and return a third       *
 ****************************/
  ComplexPoly   operator - (const ComplexPoly& x);                          // Negation operator
  ComplexPoly   operator << (const ComplexPoly& x, int numberShift);        // Shift left, 0 fill
  ComplexPoly   operator >> (const ComplexPoly& x, int numberShift);        // Shift right, 0 fill
  ComplexPoly   operator + (const ComplexPoly& x, const ComplexPoly& y);    // Adds two polynomials
  ComplexPoly   operator - (const ComplexPoly& x, const ComplexPoly& y);    // Subtracts two polynomials
  ComplexPoly   operator * (const ComplexPoly& x, const ComplexPoly& y);    // Multiplies two polynomials
  ComplexPoly   operator / (const ComplexPoly& x, const ComplexPoly& y);    // Divides two polynomials
  ComplexPoly   operator % (const ComplexPoly& x, const ComplexPoly& y);    // Remainder between two polynomials

/********************************
 * The following operators      *
 * modify a Complex by a        *
 * poly, and return a poly.     *
 ********************************/
  ComplexPoly   operator + (Complex x, const ComplexPoly &y);
  ComplexPoly   operator - (Complex x, const ComplexPoly &y);
  ComplexPoly   operator * (Complex x, const ComplexPoly &y);
  ComplexPoly   operator / (Complex x, const ComplexPoly &y);

  ComplexPoly   operator + (const ComplexPoly &x, Complex y);
  ComplexPoly   operator - (const ComplexPoly &x, Complex y);
  ComplexPoly   operator * (const ComplexPoly &x, Complex y);
  ComplexPoly   operator / (const ComplexPoly &x, Complex y);

/********************************
 * The following operators      *
 * modify a double by a         *
 * poly, and return a poly.     *
 ********************************/
  ComplexPoly   operator + (double x, const ComplexPoly &y);
  ComplexPoly   operator - (double x, const ComplexPoly &y);
  ComplexPoly   operator * (double x, const ComplexPoly &y);
  ComplexPoly   operator / (double x, const ComplexPoly &y);

  ComplexPoly   operator + (const ComplexPoly &x, double y);
  ComplexPoly   operator - (const ComplexPoly &x, double y);
  ComplexPoly   operator * (const ComplexPoly &x, double y);
  ComplexPoly   operator / (const ComplexPoly &x, double y);

/****************************************
 * The following operators perform      *
 * special functions on a polynomial.   *
 ****************************************/
  ComplexPoly pow(const ComplexPoly &x, int n);

#endif
