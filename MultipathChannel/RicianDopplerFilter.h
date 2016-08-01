#ifndef _RICIAN_DOPPLER_FILTER_H
#define _RICIAN_DOPPLER_FILTER_H    1
/************************************************************************
 *                                                                      *
 * This class implements a filter which gives the classical Doppler     *
 * spectrum with peaks at +/- the Doppler frequency.                    *
 *                                                                      *
 * File:RicianDopplerFilter.h                                           *
 *                                                                      *
 * The filter is stored in the following form:                          *
 *                                                                      *
 *           b[n] + b[n-1]z + ... + b[0]z^n                             *
 *    H(z) = ------------------------------------                       *
 *           1    + a[n-1]z + ... + a[0]z^n                             *
 *                                                                      *
 * See: D.I. Laurenson, D.G.M. Cruickshank, and G.J.R. Povey, "A        *
 * computationally efficient multipath channel simulator for the COST   *
 * 207 models," IEE Colloquium on Computer Modeling of Communication    *
 * Systems, 1994, pp. 8/1-8/6.                                          *
 *                                                                      *
 ************************************************************************/

#include        "AbstractDopplerFilter.h"


class RicianDopplerFilter: public AbstractDopplerFilter
{
private:
//
// Private functions:
//
//
// Local functions for finding the transfer function:
//
// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  RicianDopplerFilter(float doppler=DEFAULT_DOPPLER_FREQ, float sampling = DEFAULT_SAMPLING_FREQ);
  ~RicianDopplerFilter();
        
/****************************************
 * We override the filterComplexData.   *
 ****************************************/
  virtual       Complex filterComplexData(Complex &input);              // Filter complex data

};
#endif
