#ifndef _POLYNOMIAL_FIT_H
#define _POLYNOMIAL_FIT_H    1
/************************************************************************
 *                                                                      *
 * This subclass of object fits a polynomial to a set of data.          *
 *                                                                      *
 *                                                                      *
 * File:PolynomialFit.h                                                 *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 11/17/00 - Started.                                              *
 ************************************************************************/

#ifndef LOGICAL
#define LOGICAL                 char
#endif

class   DoublePoly;

class PolynomialFit
{
private:
  DoublePoly    *_fitPoly;              // Polynomial fit to data
  int           _polyOrder;             // Order of polynomial to fit
  double        _lambda;                // (A'A + lambda*I), regularization parameter
  double        _matrixDeterminant;     // Determinant of Matrix used to fit data

public:
/**********************
 * Constructors:        *
 **********************/
  PolynomialFit(int order = 3);
  ~PolynomialFit();
        
/********************************
 * PolynomialFit functions      *
 ********************************/
  LOGICAL fitToData(const float *xData, const float *yData, int numberPoints);
  const DoublePoly      *fitPoly()              {return _fitPoly;}
  const double          *fitCoefficients();

 /**********************
 * Set parameters:    *
 **********************/
  void  setPolyOrder(int order);
  void  setRegularization(double lambda)        {_lambda = lambda; return;}

/**********************
 * Get parameters:    *
 **********************/
  int           polyOrder()                     {return _polyOrder;}
  double        getDeterminant()                {return _matrixDeterminant;}


};
#endif
