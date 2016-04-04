#ifndef _FEEDFORWARDFEEDBACK_H
#define _FEEDFORWARDFEEDBACK_H    1
/************************************************************************
 *                                                                      *
 * This class solves for the feedforward and/or feedback coefficients   *
 * for a general delta sigma modulator.  It is based on the block       *
 * diagram in 5.6.4 in Delta-Sigma Data Converters: Theory, Design,     *
 * and Simulation.                                                      *
 *                                                                      *
 * File:FeedforwardFeedback.h                                           *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/23/00  - Started.                                             *
 ************************************************************************/

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX_DELTA_ORDER 3

#define DELTA_FIRST_ORDER       1
#define DELTA_SECOND_ORDER      2
#define DELTA_THIRD_ORDER       3

class Complex;                                  // class prototype

class FeedforwardFeedback
{
protected:
  int           _filterOrder;                   // Number of poles/zeros of transfer functions
  int           _oldFFCoefficients;             // Old number of FF coefficients
  int           _oldFBCoefficients;             // Old number of FB coefficients

  Complex       *_ntfZeros;                     // input zeros of noise transfer function
  Complex       *_stfZeros;                     // input zeros of signal transfer function
  Complex       *_tfPoles;                      // input poles of noise/signal transfer function

  Complex       *_ffCoefficients;               // Calculated feedforward coefficients
  Complex       *_fbCoefficients;               // Calculated feedback coefficients

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
  FeedforwardFeedback(int order=1, Complex *ntfZeros=NULL, Complex *stfZeros=NULL,
                      Complex *poles=NULL);
  ~FeedforwardFeedback();

/**********************
 * Set parameters:    *
 **********************/
  void  setFilterOrder(int order);
  void  setNTFZeros(Complex *zeros)                     {_ntfZeros      = zeros;                return;}
  void  setSTFZeros(Complex *zeros)                     {_stfZeros      = zeros;                return;}
  void  setTFPoles(Complex *poles)                      {_tfPoles       = poles;                return;}

/**********************
 * Get parameters:    *
 **********************/
  int           filterOrder()           {return _filterOrder;}

/**********************
 * Getting Outputs:   *
 **********************/
  const Complex *calculateFFCoefficients();
  const Complex *calculateFBCoefficients();

};
#endif
