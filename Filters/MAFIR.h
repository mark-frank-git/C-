#ifndef _MA_FIR_H
#define _MA_FIR_H    1
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
 ************************************************************************/

#include        "AbstractFIR.h"

#ifndef YES
#define YES     1
#define NO      0
#endif

class MAFIR: public AbstractFIR
{
private:

//
// Private functions:
//

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  MAFIR(int taps=2);
  ~MAFIR();

/**********************
 * Set parameters:    *
 **********************/
  virtual void  setNumberTaps(int taps);
  virtual void  setFilterOrder(int order);


};
#endif
