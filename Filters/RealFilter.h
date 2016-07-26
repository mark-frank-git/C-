#ifndef _REAL_FILTER_H
#define _REAL_FILTER_H 1
/************************************************************************
 *                                                                      *
 * This subclass of AbstractFilter is for filters with real             *
 * real coefficients.                                                   *
 *                                                                      *
 * File:RealFilter.h                                                    *
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
 * or equivalently, for digital filters:                                *
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
 *  1. 10/08/00  - Started as subclass of AbstractFilter                *
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

class   DoublePoly;                             // Class prototypes
class   Complex;
class   TestData;

class RealFilter: public AbstractFilter
{
protected:
  DoublePoly    *_aPolyObject;                  // Denominator polynomial
  DoublePoly    *_bPolyObject;                  // Numerator polynomial
  TestData      *_fileReader;                   // Reads coefficients from file
  int           _numberFileCoeffs;              // Number of read in file coefficients
  int           _oldNumberCoeffs;               // Old number of read in coeffs

  double        *_fileCoeffs;                   // Read in file coefficients

//
// Private methods:
//
  void          transferFromPoles(Complex *poles, Complex *zeros, int numPoles, int numZeros,
                                  LOGICAL realAxisPoles=NO);
                                                // Finds the transfer function from pole-zero form
  virtual       Complex transferResponseAt(Complex &point);
                                                // Return transfer function response at a point
  LOGICAL       readCoeffFile(const char *fileName, int fileType);
  LOGICAL       writeCoeffFile(const char *fileName, int fileType, const double *coeffs, int numberCoeffs);
  LOGICAL       writePoleZeroToFile(const char *fileName, int fileType, const Complex *roots, int numberRoots);

//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  RealFilter(double *aCoeffs, double *bCoeffs, int order);
  RealFilter(int type=LOW_PASS, double centerFreq=0., double cutoffFreq=10., int order=2);
  virtual  ~RealFilter();

/**************************
 * Setting parameters:    *
 **************************/
  virtual       void setFilterACoeffs(const double *a);
  virtual       void setFilterBCoeffs(const double *b);
  virtual       LOGICAL readACoeffsFromFile(const char *fileName, int fileType);
  virtual       LOGICAL readBCoeffsFromFile(const char *fileName, int fileType);
  virtual       LOGICAL writeACoeffsToFile(const char *fileName, int fileType);
  virtual       LOGICAL writeBCoeffsToFile(const char *fileName, int fileType);
  virtual       LOGICAL writePolesToFile(const char *fileName, int fileType);
  virtual       LOGICAL writeZerosToFile(const char *fileName, int fileType);
  virtual       void setFilterOrder(int order);

/**********************
 * Get parameters:    *
 **********************/
  int   aOrder();
  int   bOrder();
  const double *aCoeffs();
  const double *bCoeffs();
  const DoublePoly *aPoly()             {return _aPolyObject;}
  const DoublePoly *bPoly()             {return _bPolyObject;}
  virtual const Complex *filterPoles()  {return NULL;}  // override in subclass
  virtual const Complex *filterZeros()  {return NULL;}  // override in subclass

};

#endif
