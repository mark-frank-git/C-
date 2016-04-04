#ifndef _POLYNOMIAL_H
#define _POLYNOMIAL_H   1
/************************************************************************
 *                                                                      *
 * This subclass of object implements a polynomial object.  The         *
 * coefficients of the polynomial are double precision variables.       *
 * The usual operations are defined for this class, +, x, -, /, +=, etc.*
 * This class should be useful for digital signal processing.           *
 *                                                                      *
 * File:Polynomial.h                                                    *
 *                                                                      *
 * Note: for p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]              *
 *       p(x) has order n, but the length of the array is n+1.          *
 *                                                                      *
 * Revision history:                                                    *
 *  3. 02/18/97  - Started                                              *
 ************************************************************************/
/********************** Won't compile in NT ***************************************
#include <fstream.h>
 **********************************************************************************/
#include <stdio.h>

class   Complex;                                    // Class prototype

#ifndef LOGICAL
#define LOGICAL char
#endif


#define MAX_SHIFT           1024                    // Maximum shift right or left
#define MAX_STRING_SIZE     2048                    // Polynomial string size

class Polynomial
{
private:
  int       order;                                  // orders of the polynomial
  int       size;                                   // actual length of the array
  double    *coefficients;                          // array containing the polynomial coefficients
  char      *polyString;                            // E.g., "0.5*x^n + -.2*x^n-1 ..."

  Complex   *polyRoots;                             // Complex roots of a polynomial
//
// Private methods:
//
  void updateMemory (int newOrder);                     // Adjust memory size
  void updateOrder();                                   // Check that x[order+1] != 0, else update

public:
//
// Public methods:
//
/*********************************
 * Constructors, destructors:    *
 *********************************/
  Polynomial(int pOrder=0, const double *coeff=NULL);           // Initializes a polynomial
  Polynomial(const Polynomial &poly);
  Polynomial(double x);
  ~Polynomial();

/**********************
 * Set parameters:    *
 **********************/
  void assign(int pOrder, const double *coeff);                 // Sets a new polynomial

/**********************
 * Get parameters:    *
 **********************/
  const double *getCoefficients() const {return coefficients;}  // Return the polynomial coefficients
  int getOrder() const                  {return order;}         // Return the order of the polynomial

/**********************
 * Calculating:       *
 **********************/
  void reverse();                                           // Reverse the order of coefficients of A
  double  evaluateAtRealPoint(double point);
  Complex evaluateAtComplexPoint(Complex &point) const;     // Evaluate the A polynomial at a complex point
  Complex *roots();                                         // Evaluate and return the roots of A

/****************************
 * The following functions  *
 * aid in printing the      *
 * polynomial               *
 ****************************/
  const char *getPolyString();

/****************************
 * The following operators  *
 * modify the polynomial    *
 * by a second polynomial   *
 ****************************/
  Polynomial&   operator =  (const Polynomial& y);
  Polynomial&   operator =  (double y);
  Polynomial&   operator += (const Polynomial& y);
  Polynomial&   operator -= (const Polynomial& y);
  Polynomial&   operator *= (const Polynomial& y);
  Polynomial&   operator /= (const Polynomial& y);
  Polynomial&   operator %= (const Polynomial& y);
  Polynomial&   operator >>= (int numberShifts);
  Polynomial&   operator <<= (int numberShifts);

/****************************
 * The following operators  *
 * modify the polynomial    *
 * by a double.             *
 ****************************/
  Polynomial&   operator += (double y);
  Polynomial&   operator -= (double y);
  Polynomial&   operator *= (double y);
  Polynomial&   operator /= (double y);

};
//
// The following are nonmember functions taking two arguments
//

/****************************
 * The following operators  *
 * are for printing the     *
 * polynomial string.       *
 ****************************/
/********************** Won't compile in NT ********************************************
  ostream&  operator << (ostream& s, const Polynomial& x);        // Outputs a polynomial
 ***************************************************************************************/

/****************************
 * The following operators  *
 * make comparisons between *
 * two polynomials.         *
 ****************************/
  LOGICAL operator == (const Polynomial& x, const Polynomial& y);       // Checks if two polynomials are equal
  LOGICAL operator != (const Polynomial& x, const Polynomial& y);       // Checks if two polynomials are not equal
  LOGICAL operator >  (const Polynomial& x, const Polynomial& y);       // Checks if x > y
  LOGICAL operator >= (const Polynomial& x, const Polynomial& y);       // Checks if x >= y
  LOGICAL operator <  (const Polynomial& x, const Polynomial& y);       // Checks if x < y
  LOGICAL operator <= (const Polynomial& x, const Polynomial& y);       // Checks if x <= y

/****************************
 * The following operators  *
 * take 1 or 2 polynomials  *
 * and return a third       *
 ****************************/
  Polynomial    operator - (const Polynomial& x);                           // Negation operator
  Polynomial    operator >> (const Polynomial& x, int numberShift);         // Shift right, 0 fill
  Polynomial    operator << (const Polynomial& x, int numberShift);         // Shift left, 0 fill
  Polynomial    operator + (const Polynomial& x, const Polynomial& y);      // Adds two polynomials
  Polynomial    operator - (const Polynomial& x, const Polynomial& y);      // Subtracts two polynomials
  Polynomial    operator * (const Polynomial& x, const Polynomial& y);      // Multiplies two polynomials
  Polynomial    operator / (const Polynomial& x, const Polynomial& y);      // Divides two polynomials
  Polynomial    operator % (const Polynomial& x, const Polynomial& y);      // Remainder between two polynomials


/****************************
 * The following operators  *
 * modify a double by a     *
 * polynomial, and return a *
 * polynomial.              *
 ****************************/
  Polynomial    operator + (double x, const Polynomial &y);
  Polynomial    operator - (double x, const Polynomial &y);
  Polynomial    operator * (double x, const Polynomial &y);
  Polynomial    operator / (double x, const Polynomial &y);

  Polynomial    operator + (const Polynomial &x, double y);
  Polynomial    operator - (const Polynomial &x, double y);
  Polynomial    operator * (const Polynomial &x, double y);
  Polynomial    operator / (const Polynomial &x, double y);

/********************************
 * The following operators      *
 * perform special function     *
 * on a polynomial.             *
 ********************************/
  Polynomial    pow (const Polynomial &x, int n);


#endif

