#ifndef _DOUBLE_POLY_H
#define _DOUBLE_POLY_H   1
/************************************************************************
 *                                                                      *
 * This subclass of object implements a polynomial object.  The         *
 * coefficients of the polynomial are double precision variables.       *
 * The usual operations are defined for this class, +, x, -, /, +=, etc.*
 * This class should be useful for digital signal processing.           *
 *                                                                      *
 * File:DoublePoly.h                                                    *
 *                                                                      *
 * Note: for p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]              *
 *       p(x) has order n, but the length of the array is n+1.          *
 *                                                                      *
 ************************************************************************/

#ifdef USE_IO_STREAM
#include <fstream.h>
#endif

#include <stdio.h>

class   Complex;                                    // Class prototype

#ifndef LOGICAL
#define LOGICAL char
#endif


#define MAX_SHIFT           1024                    // Maximum shift right or left
#define MAX_STRING_SIZE     2048                    // Polynomial string size

class DoublePoly
{
private:
  int       _order;                                  // orders of the polynomial
  int       _size;                                   // actual length of the array
  double    *_coefficients;                          // array containing the polynomial coefficients
  char      *_polyString;                            // E.g., "0.5*x^n + -.2*x^n-1 ..."

  Complex   *_polyRoots;                             // Complex roots of a polynomial
//
// Private methods:
//
  void updateMemory (int newOrder);                     // Adjust memory size
  double maxCoeffNorm ();
  void updateOrder();                                   // Check that x[order+1] != 0, else update

public:
//
// Public methods:
//
/*********************************
 * Constructors, destructors:    *
 *********************************/
  DoublePoly(int pOrder=0, const double *coeff=NULL);           // Initializes a polynomial
  DoublePoly(const DoublePoly &poly);
  DoublePoly(double x);
  ~DoublePoly();

/**********************
 * Set parameters:    *
 **********************/
  void assign(int pOrder, const double *coeff);                 // Sets a new polynomial
  void assignRoots(int pOrder, const Complex *roots);           // Sets a polynomial from roots

/**********************
 * Get parameters:    *
 **********************/
  const double *getCoefficients() const {return _coefficients;}         // Return the polynomial coefficients
  int getOrder() const                  {return _order;}                // Return the order of the polynomial

/**********************
 * Calculating:       *
 **********************/
  void reverse();                                           // Reverse the order of coefficients of A
  double  evaluateAtRealPoint(double point);
  Complex evaluateAtComplexPoint(Complex &point) const;     // Evaluate the polynomial at a complex point
  Complex *roots();                                         // Evaluate and return the roots of the polynomial

/****************************
 * The following operators  *
 * are for printing the     *
 * polynomial string.       *
 ****************************/
  const char *getPolyString();

/****************************
 * The following operators  *
 * modify the polynomial    *
 * by a second polynomial   *
 ****************************/
  DoublePoly&   operator =  (const DoublePoly& y);
  DoublePoly&   operator =  (double y);
  DoublePoly&   operator += (const DoublePoly& y);
  DoublePoly&   operator -= (const DoublePoly& y);
  DoublePoly&   operator *= (const DoublePoly& y);
  DoublePoly&   operator /= (const DoublePoly& y);
  DoublePoly&   operator %= (const DoublePoly& y);
  DoublePoly&   operator >>= (int numberShifts);
  DoublePoly&   operator <<= (int numberShifts);

/****************************
 * The following operators  *
 * modify the polynomial    *
 * by a double.             *
 ****************************/
  DoublePoly&   operator += (double y);
  DoublePoly&   operator -= (double y);
  DoublePoly&   operator *= (double y);
  DoublePoly&   operator /= (double y);

};

/****************************
 * The following operators  *
 * are for printing the     *
 * polynomial string.       *
 ****************************/
#ifdef USE_IO_STREAM
  ostream&  operator << (ostream& s, const DoublePoly& x);        // Outputs a polynomial
#endif
//
// The following are nonmember functions taking two arguments
//

/****************************
 * The following operators  *
 * make comparisons between *
 * two polynomials.         *
 ****************************/
  LOGICAL operator == (const DoublePoly& x, const DoublePoly& y);       // Checks if two polynomials are equal
  LOGICAL operator != (const DoublePoly& x, const DoublePoly& y);       // Checks if two polynomials are not equal
  LOGICAL operator >  (const DoublePoly& x, const DoublePoly& y);       // Checks if x > y
  LOGICAL operator >= (const DoublePoly& x, const DoublePoly& y);       // Checks if x >= y
  LOGICAL operator <  (const DoublePoly& x, const DoublePoly& y);       // Checks if x < y
  LOGICAL operator <= (const DoublePoly& x, const DoublePoly& y);       // Checks if x <= y

/****************************
 * The following operators  *
 * take 1 or 2 polynomials  *
 * and return a third       *
 ****************************/
  DoublePoly    operator - (const DoublePoly& x);                           // Negation operator
  DoublePoly    operator >> (const DoublePoly& x, int numberShift);         // Shift right, 0 fill
  DoublePoly    operator << (const DoublePoly& x, int numberShift);         // Shift left, 0 fill
  DoublePoly    operator + (const DoublePoly& x, const DoublePoly& y);      // Adds two polynomials
  DoublePoly    operator - (const DoublePoly& x, const DoublePoly& y);      // Subtracts two polynomials
  DoublePoly    operator * (const DoublePoly& x, const DoublePoly& y);      // Multiplies two polynomials
  DoublePoly    operator / (const DoublePoly& x, const DoublePoly& y);      // Divides two polynomials
  DoublePoly    operator % (const DoublePoly& x, const DoublePoly& y);      // Remainder between two polynomials


/****************************
 * The following operators  *
 * modify a double by a     *
 * polynomial, and return a *
 * polynomial.              *
 ****************************/
  DoublePoly    operator + (double x, const DoublePoly &y);
  DoublePoly    operator - (double x, const DoublePoly &y);
  DoublePoly    operator * (double x, const DoublePoly &y);
  DoublePoly    operator / (double x, const DoublePoly &y);

  DoublePoly    operator + (const DoublePoly &x, double y);
  DoublePoly    operator - (const DoublePoly &x, double y);
  DoublePoly    operator * (const DoublePoly &x, double y);
  DoublePoly    operator / (const DoublePoly &x, double y);

/********************************
 * The following operators      *
 * perform special function     *
 * on a polynomial.             *
 ********************************/
  DoublePoly    pow (const DoublePoly &x, int n);


#endif

