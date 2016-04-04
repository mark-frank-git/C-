/************************************************************************
 *                                                                      *
 * This subclass of AbstractFIR implements a moving average type        *
 * of FIR filter, i.e., all coefficients = 1.                           *
 * filters.                                                             *
 *                                                                      *
 * File:MAFIR.h                                                         *
 *                                                                      *
 * The filter is stored in the forms:                                   *
 *                                                                      *
 *    H(z) =  b[n] + b[n-1]z + ... + b[0]z^n                            *
 *         =  b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                      *
 *                                                                      *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/17/00 - Started.                                              *
 ************************************************************************/

#include "MAFIR.h"                                      // Object prototypes
#include <Polynomials/DoublePoly.h>
#include <math.h>
#include <stdio.h>
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the MAFIR class.
//
// Input:           taps:       Number of filter taps
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
MAFIR::MAFIR(int taps)
              :AbstractFIR()
{
  setNumberTaps(taps);
  
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the MAFIR class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
MAFIR::~MAFIR()
{
  return;
}

// ############################# Public Method ###############################
// setNumberTaps -- Sets a new number of taps for the filter.
// Input:       taps            New number of filter taps
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void MAFIR::setNumberTaps(int taps)
{
  int           i;
  double        *tap_coeffs;
//
// Call super's method
//
  AbstractFIR::setNumberTaps(taps);

  tap_coeffs    = new double[_numberTaps];
  for(i=0; i<_numberTaps; i++)
    tap_coeffs[i]       = 1./_numberTaps;
  _bPolyObject->assign(_numberTaps-1, tap_coeffs);
  delete [] tap_coeffs;
  return;
}

// ############################# Public Method ###############################
// setFilterOrder -- Sets the order of the filter (numberTaps - 1).
// Input:       taps            New number of filter taps
//          
// Output:                      none
// Notes:
// ############################# Public Method ###############################
void MAFIR::setFilterOrder(int order)
{
  setNumberTaps(order+1);
  return;
}

