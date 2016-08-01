#ifndef _DELAYFIR_H
#define _DELAYFIR_H     1
/************************************************************************
 *                                                                      *
 * This subclass of AbstractFIR adds functionality for calculating      *
 * the filter coefficients for implmenting delay FIR filters.           *
 * The code is based on the article, "Splitting the unit delay,"        *
 * T. Laakso, V. Valimaki, M. Karjalainen, and U. Laine, IEEE Sig Proc  *
 * Magazine, Jan. 1996.                                                 *
 *                                                                      *
 * File:DelayFIR.h                                                      *
 *                                                                      *
 * The filter is stored in the forms:                                   *
 *                                                                      *
 *    H(z) =  b[n] + b[n-1]z + ... + b[0]z^n                            *
 *         =  b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                      *
 *                                                                      *
 * Where b polynomials are stored in the polynomial objects.            *
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 ************************************************************************/

#include        "AbstractFIR.h"

#define         MAX_DELAY       1.
#define         MIN_DELAY       -1.

#define         LAGRANGE_DELAY  0               // Delay types

#ifndef YES
#define YES     1
#define NO      0
#endif

#ifndef LOGICAL
#define LOGICAL char
#endif

class Complex;                                  // class prototype

class DelayFIR: public AbstractFIR
{
private:
  float         _normalizedDelay;                       // Delay relative to 1 sampling period

//
// Private functions:
//
  LOGICAL       transferForLagrange();                  // calculates transfer function (tap weights)

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  DelayFIR(int taps = 7, float samplingFreq=100., float delay=0.25);
  ~DelayFIR();

/**********************
 * Set parameters:    *
 **********************/
  void  setDelay(float delay);

/**********************
 * Get parameters:    *
 **********************/
  float delay()                                 {return _normalizedDelay;}
        
/****************************************
 * Filtering float and Complex data     *
 * assuming an FIR structure.  These    *
 * are defined in the superclass.       *
 ****************************************/

/**********************
 * Getting Outputs:   *
 **********************/
  float         findTotalDelay();
  void          findTransferFunction();                         // Find transfer function

};
#endif
