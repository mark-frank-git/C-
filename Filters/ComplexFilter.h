#ifndef _COMPLEX_FILTER_H
#define _COMPLEX_FILTER_H 1
/************************************************************************
 *                                                                      *
 * This subclass of AbstractFilter is for filters with Complex          *
 * coefficients.                                                        *
 *                                                                      *
 * File:ComplexFilter.h                                                 *
 *                                                                      *
 * The filter is stored in two forms:                                   *
 *                                                                      *
 *          b[n] + b[n-1]s + ... + b[0]s**(n)                           *
 *    H(s) = ------------------------------------                       *
 *          1    + a[n-1]s + ... + a[0]s**(n)                           *
 *                                                                      *
 *           (s-zero[0]) * (s-zero[1]) ... (s-zero[n_zero])             *
 *    H(s) = ----------------------------------------------             *
 *           (s-pole[0]) * (s-pole[1]) ... (s-pole[n_pole])             *
 *                                                                      *
 * Or, equivalently for digital filters:                                *
 *                                                                      *
 *           b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                       *
 *  H(z)   = ------------------------------------                       *
 *           1    + a[1]z^(-1) + ... + a[n]z^-(n)                       *
 *                                                                      *
 *  and:                                                                *
 *           (z-zero[0]) * (z-zero[1]) ... (z-zero[n_zero])             *
 *    H(z) = ----------------------------------------------             *
 *           (z-pole[0]) * (z-pole[1]) ... (z-pole[n_pole])             *
 *                                                                      *
 * NOTE: The order of the filter is n, but the number of coefficients   *
 *       is equal to n+1.                                               *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 10/08/99  - Started as subclass of AbstractFilter                *
 ************************************************************************/
#include "AbstractFilter.h"
#include <stdio.h>

#ifndef LOGICAL
#define LOGICAL char
#endif
#ifndef YES
#define YES     1
#define NO      0
#endif

class   ComplexPoly;                            // Class prototypes
class   Complex;
class   TestData;

class ComplexFilter: public AbstractFilter
{
protected:
  ComplexPoly   *_aPolyObject;                  // Denominator polynomial
  ComplexPoly   *_bPolyObject;                  // Numerator polynomial
  TestData      *_fileReader;                   // Reads coefficients from file
  int           _numberFileData;                // Number of read in file coefficients
  int           _oldNumberData;                 // Old number of read in coeffs

  Complex       *_fileData;                     // Read in file coefficients/poles/zeros

//
// Private methods:
//
  virtual       Complex transferResponseAt(Complex &point);
                                                // Return transfer function response at a point
  LOGICAL       readDataFile(const char *fileName, int fileType, LOGICAL realCoeff);
  float         denominatorCenterOfMass();      // Calculate center of mass of denominator coefficients
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  ComplexFilter(Complex *aCoeffs, Complex *bCoeffs, int order);
  ComplexFilter(int type=LOW_PASS, double centerFreq=0., double cutoffFreq=10., int order=2);
  virtual  ~ComplexFilter();

/**************************
 * Setting parameters:    *
 **************************/
  virtual       void setFilterACoeffs(const Complex *a);
  virtual       void setFilterBCoeffs(const Complex *b);
  virtual       void setFilterPoles(int number, const Complex *poles);
  virtual       void setFilterZeros(int number, const Complex *zeros);
  virtual       LOGICAL readACoeffsFromFile(const char *fileName, int fileType, LOGICAL realCoeff=NO);
  virtual       LOGICAL readBCoeffsFromFile(const char *fileName, int fileType, LOGICAL realCoeff=NO);
  virtual       LOGICAL readPolesFromFile(const char *fileName, int fileType);
  virtual       LOGICAL readZerosFromFile(const char *fileName, int fileType);
  virtual       void setFilterOrder(int order);

/**********************
 * Get parameters:    *
 **********************/
  int   aOrder();
  int   bOrder();
  const Complex *aCoeffs();
  const Complex *bCoeffs();

};

#endif
