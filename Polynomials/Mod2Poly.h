#ifndef _MOD2POLY_H
#define _MOD2POLY_H     1
/************************************************************************
 *                                                                      *
 * This subclass of object implements a polynomial object.  Operations  *
 * on the polynomial are taken over the binary number field, GF(2).     *
 * The usual operations are defined for this class, +, x, -, /, +=, etc.*
 * This class should be useful for analyzing linear block codes.        *
 *                                                                      *
 * File:Mod2Poly.h                                                      *
 *                                                                      *
 * Note: for p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]              *
 *       p(x) has order n, but the size of the array is n+1.            *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 12/22/95 - Started                                               *
 *  2. 09/18/96 - Added ^ and & operators.                              *
 *  3. 12/11/96 - Add getPolyString().                                  *
 *  4. 03/10/97 - Reverse order of polynomial.                          *
 ************************************************************************/
#ifdef USE_IO_STREAM
#include <fstream.h>
#endif

#ifndef NULL
#define NULL                            0
#endif
#ifndef LOGICAL
#define LOGICAL                         char
#endif
#ifndef YES
#define YES                             1
#define NO                              0
#endif
#ifndef BIT
#define BIT                             short
#endif
#define MOD_2_ADD                       ^                       // exclusive or for modulo 2 additions
#define MOD_2_PLUS_EQUAL                ^=
#define MOD_2_MULT                      &                       // and for modulo 2 multiplications
#define MOD_2_TIMES_EQUAL               &=
#define MAX_SHIFT                       1024                    // Maximum shift right or left
#define MAX_STRING_SIZE                 2048                    // Polynomial string size

class Mod2Poly
{
protected:
  int   _order;                                         // order of the polynomial
  int   _size;                                          // size of the allocated memory for coefficients
  BIT   *_coefficients;                                 // arrays containing the polynomial coefficients
  char  *_polyString;                                   // E.g., "1 +  + x^2 + x^3 +  + x^5"

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
  Mod2Poly(int pOrder=0, const BIT *coeff=NULL);        // Initializes a polynomial
  Mod2Poly(const Mod2Poly& poly);
  Mod2Poly(BIT input);
  ~Mod2Poly();

/********************************
 * The following functions      *
 * set parameters.              *
 ********************************/
  void assign(int pOrder, const BIT *coeff);            // Sets a new polynomial

/********************************
 * The following functions      *
 * get parameters.              *
 ********************************/
  const BIT *getCoefficients()  const {return _coefficients;}   // Returns the coefficients
  int getOrder()                const {return _order;}          // Return the order of the polynomial

/********************************
 * Calculating:                 *
 ********************************/
  int hammingWeight();                                          // Find the number of non zero coefficients
  void reverse();

/********************************
 * The following functions      *
 * print the polynomial         *
 ********************************/
  const char *getPolyString();

/********************************
 * The following operators      *
 * modify the polynomial        *
 * by a second polynomial       *
 ********************************/
  Mod2Poly&     operator =  (const Mod2Poly& y);
  Mod2Poly&     operator =  (BIT y);
  Mod2Poly&     operator += (const Mod2Poly& y);
  Mod2Poly&     operator -= (const Mod2Poly& y);
  Mod2Poly&     operator *= (const Mod2Poly& y);
  Mod2Poly&     operator /= (const Mod2Poly& y);
  Mod2Poly&     operator %= (const Mod2Poly& y);
  Mod2Poly&     operator ^= (const Mod2Poly& y);                                                // Same as +=
  Mod2Poly&     operator &= (const Mod2Poly& y);
  Mod2Poly&     operator >>= (int numberShifts);
  Mod2Poly&     operator <<= (int numberShifts);

};
//
// The following are nonmember functions taking two arguments
//

/****************************
 * The following operators      *
 * are for printing the         *
 * polynomial string.           *
 ****************************/
#ifdef IO_STREAM
  ostream&  operator << (ostream& s, Mod2Poly& x);                      // Outputs a polynomial
#endif

/****************************
 * The following operators      *
 * make comparisons between     *
 * two polynomials.             *
 ****************************/
  LOGICAL       operator == (const Mod2Poly& x, const Mod2Poly& y);     // Checks if two polynomials are equal
  LOGICAL       operator != (const Mod2Poly& x, const Mod2Poly& y);     // Checks if two polynomials are not equal
  LOGICAL       operator >  (const Mod2Poly& x, const Mod2Poly& y);     // Checks if x > y
  LOGICAL       operator >= (const Mod2Poly& x, const Mod2Poly& y);     // Checks if x >= y
  LOGICAL       operator <  (const Mod2Poly& x, const Mod2Poly& y);     // Checks if x < y
  LOGICAL       operator <= (const Mod2Poly& x, const Mod2Poly& y);     // Checks if x <= y
  int   hammingDistance(const Mod2Poly& x, const Mod2Poly& y);          // Returns Hamming distance btwn 2 polys

/****************************
 * The following operators      *
 * take 1 or 2 polynomials      *
 * and return a third           *
 ****************************/
  Mod2Poly      operator - (const Mod2Poly& x);                         // Negation operator
  Mod2Poly      operator << (const Mod2Poly& x, int numberShift);       // Shift left, 0 fill
  Mod2Poly      operator >> (const Mod2Poly& x, int numberShift);       // Shift right, 0 fill
  Mod2Poly      operator + (const Mod2Poly& x, const Mod2Poly& y);      // Adds two polynomials
  Mod2Poly      operator - (const Mod2Poly& x, const Mod2Poly& y);      // Subtracts two polynomials
  Mod2Poly      operator * (const Mod2Poly& x, const Mod2Poly& y);      // Multiplies two polynomials
  Mod2Poly      operator / (const Mod2Poly& x, const Mod2Poly& y);      // Divides two polynomials
  Mod2Poly      operator % (const Mod2Poly& x, const Mod2Poly& y);      // Remainder between two polynomials
  Mod2Poly      operator ^ (const Mod2Poly& x, const Mod2Poly& y);      // Bit wise exclusive or same as +
  Mod2Poly      operator & (const Mod2Poly& x, const Mod2Poly& y);      // Bit wise and of two polys

#endif