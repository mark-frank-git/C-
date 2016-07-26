/************************************************************************
 *                                                                      *
 * This subclass of AbstractFIR adds functionality for calculating      *
 * the filter coefficients of certain types of FIR                      *
 * filters.                                                             *
 *                                                                      *
 * File:FIRFilter.cc                                                    *
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
 *  1. 01/06/00 - Derived from DigitalFilter.                           *
 ************************************************************************/

#include "FIRFilter.h"                                  // Object prototypes
#include "DataWindow.h"

#if defined(WIN32)
#include <GNU/Complex.h>
#include <Polynomials/DoublePoly.h>
#include <C_Libraries/constants.h>
#include <Specfuns/specfuns.h>
#else
#include "Complex.h"
#include "DoublePoly.h"
#include "constants.h"
#include "specfuns.h"
#endif

#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

#define EDGE_COEFF1     1.045                           // EDGE pulse shape coefficients
#define EDGE_COEFF2     0.218

// ############################# Private Method ###############################
// windowImpulseResponse -- This routine finds the impulse response for either
//                          an ideal brick wall filter or a raised cosine filter.
//
// Input:       i:              Index of impulse response
//              alpha:          (number_taps-1)/2
//              wc:             cutoff frequency in radians
//              w1:             bandpass/stop low frequency in radians
//              w2:             bandpass/stop high frequency in radians
//          
// Output:                      Returns the impulse response at the given index
//
// Notes:
//  1. The impulse response of the raised cosine filter is from Jeruchim et al p. 345
//  2. The impulse response of the root raised cosine filter is from Jeruchim et al p. 346
//  2. The impulse response of the raised cosine for non-low pass are TBD.
// ############################# Private Method ###############################
double FIRFilter::windowImpulseResponse(int i, double alpha, double wc, double w1, double w2)
{
  double output;
  double i_minus_alpha, num, den, beta, rb;

  i_minus_alpha = (double)i - alpha;
  if(_firFilterType == FIR_WINDOW)
  {
    switch(_filterPassType)
    {
      case LOW_PASS:
      default:
        if(ABS(i_minus_alpha) > EPS)
          output    = sin(wc*(i-alpha))/PI/(i-alpha);
        else
          output    = wc/PI;
        break;
      case HIGH_PASS:
        if(ABS(i_minus_alpha) > EPS)
          output    = (sin(PI*(i-alpha)) - sin(wc*(i-alpha)))/PI/(i-alpha);
        else
          output    = (PI-wc)/PI;
        break;
      case BAND_PASS:
        if(ABS(i_minus_alpha) > EPS)
          output    = (sin(w2*(i-alpha)) - sin(w1*(i-alpha)))/PI/(i-alpha);
        else
          output    = (w2-w1)/PI;
        break;
      case BAND_STOP:
        if(ABS(i_minus_alpha) > EPS)
          output    = (sin(w1*(i-alpha)) + sin(PI*(i-alpha)) - sin(w2*(i-alpha)))/PI/(i-alpha);
         else
          output    = (TWOPI - 2*(w2-w1))/TWOPI;
         break;
    }
  }
  else if(_firFilterType == FIR_RAISED_COS)                     // FIR raised cosine
  {
    switch(_filterPassType)
    {
      case LOW_PASS:                                            // Only low pass implemented
      default:
        if(ABS(i_minus_alpha) > EPS)
          output    = sin(PI*_symbolsPerSample*i_minus_alpha)/(PI*_symbolsPerSample*i_minus_alpha);
        else
          output    = 1.;
        num     = cos(PI*_raisedCosineAlpha*_symbolsPerSample*i_minus_alpha);
        den     = 2.*_raisedCosineAlpha*_symbolsPerSample*i_minus_alpha;
        den     = 1. - den*den;
        if(den != 0.)
          output    *= num/den;
 
    }
  }
  else                                                          // FIR root raised cosine
  {
    switch(_filterPassType)
    {
      case LOW_PASS:                                            // Only low pass implemented
      default:
        beta    = _raisedCosineAlpha*_symbolsPerSample/2.;
        rb      = _symbolsPerSample;
        if(ABS(i_minus_alpha) > EPS)
          num   = sin(PI*(rb - 2.*beta)*i_minus_alpha)/(8.*beta*i_minus_alpha);
        else
          num   = (rb - 2.*beta)*PI/8./beta;
        num     += cos(PI*(rb + 2.*beta)*i_minus_alpha);
        den     = PI*( (8.*beta*i_minus_alpha)*(8.*beta*i_minus_alpha) - 1.)/sqrt(rb);
        if(den != 0.)
          output    = num/den;
        else
          output    = 1.;
 
    }
  }
  return output;
}

// ############################# Private Method ###############################
// transferForFIRWindow -- This routine finds the coefficients of the filter
//                         found by multiplying the window function by either
//                         a sinc pulse whose spectrum is an ideal filter, or
//                         the impulse response of a raised cosine filter.
//
// Input:                   None
//          
// Output:                  Return YES if OK, return NO if polynomials don't exist
//
// Notes:
//  1. The instanceVariable, _bPolyObject is modified.
//  2. The number of taps is equal to order + 1.
//  3. The impulse response of the raised cosine filter is from Leon Couch p. 167
// ############################# Private Method ###############################
LOGICAL FIRFilter::transferForFIRWindow()
{
  int       i, n;
  const     double *window;
  double    wo, wc, w1, w2, beta;
  double    *b, b_max;
  
  if(_bPolyObject == NULL)
     return NO;
   
  n             = _filterOrder+1;
  b             = new double[n];
  wo            = TWOPI*_fo/_fs;
  wc            = TWOPI*_fc/_fs;
  w1            = wo - wc;
  w2            = wo + wc;

//
// Calculate the window function, e.g., Hamming
//
  window        = _dataWindow->windowFunction(n);
//
// Now, modulate to frequency of interest
//
  b_max         = 0.;
  beta          = (double)(n-1)/2.;
  for(i=0; i<n; i++)
  {
     b[i]   = windowImpulseResponse(i, beta, wc, w1, w2) * window[i];
     b_max  = MAX(b_max, b[i]);
  }
//
// Normalize coefficients, and save in poly object
//
  if(b_max > 0.)
  {
    for(i=0; i<n; i++)
    {
      b[i]  /= b_max;
    }
  }
  _bPolyObject->assign(_filterOrder, b);
  delete [] b;
//
// Set unity pass band gain:
//
  setPassBandGain();

  return YES;
}

// ############################# Private Method ###############################
// transferForPowerLawFilter -- This routine finds the coefficients of the FIR
//                              filter for implementing a 1/(f^alpha) power law
//                              response.
//
// Input:                   None
//          
// Output:                  Return YES if OK, return NO if polynomials don't exist
//
// Notes:
//  1. The instanceVariables, _bPolyObject are modified.
//  2. The number of taps is equal to order + 1.
//  3. The instance variables, _powerLawAlpha  affect
//     the output, in addition to the normal filter parameters.
//  4. Based on N.J. Kasdin, "Discrete simulation of ...", Proc. IEEE, vol. 83, no.5
// ############################# Private Method ###############################
LOGICAL FIRFilter::transferForPowerLawFilter()
{
  int       n, old_type;
  double    *b;
  
  b             = new double[_numberTaps];

//
// Set up the filter coefficients according to equation (104) in the reference
//
  if(_numberTaps)
  {
    b[0]        = 1.;
    for(n=1; n<_numberTaps; n++)
      b[n]      = (_powerLawAlpha/2. + n - 1)*b[n-1]/n;
//
// Store coefficients in polynomials:
//
    _bPolyObject->assign(_filterOrder, b);
    delete [] b;
//
// Set unity pass band gain:
//
    old_type            = _filterPassType;
    _filterPassType     = LOW_PASS;
    setPassBandGain();
    _filterPassType     = old_type;
  }

  return YES;
}

// ############################# Private Method ###############################
// transferForGMSKFilter -- This routine finds the coefficients of the FIR
//                          filter for implementing a GMSK type of filter, taken
//                          from Brad Badke's paper.
//
// Input:                   None
//          
// Output:                  Return YES if OK, return NO if polynomials don't exist
//
// Notes:
//  1. The instanceVariables, _bPolyObject are modified.
//  2. The number of taps is equal to order + 1.
// ############################# Private Method ###############################
LOGICAL FIRFilter::transferForGMSKFilter()
{
  int           n;
  float         time;
  double        *b;
  
  b             = new double[_numberTaps];

//
// Set up the filter coefficients according to equation (2) in the reference
// The center time (3) should be centered in the _numberTaps:
// Let T be normalized to 1.
//
  time          = GMSK_CENTER - _symbolsPerSample*(_numberTaps-1)/2.;
  for(n=0; n<_numberTaps; n++)
  {
    b[n]        = qx(GMSK_K1*_gmskBT*(-time)) - qx(GMSK_K1*_gmskBT*(1.-time));
    time        += _symbolsPerSample;
  }
//
// Store coefficients in polynomials:
//
  _bPolyObject->assign(_filterOrder, b);
  delete [] b;
  setPassBandGain();
  return YES;
}

// ############################# Private Method ###############################
// transferForEDGEFilter -- This routine finds the coefficients of the FIR
//                          filter for implementing an EDGE type of filter, taken
//                          from Tropian's paper at www.tropian.com/tech/tech_docs/
//                          edge_paper.pdf.
//
// Input:                   None
//          
// Output:                  Return YES if OK, return NO if polynomials don't exist
//
// Notes:
//  1. The instanceVariables, _bPolyObject are modified.
//  2. The number of taps is equal to order + 1.
//  3. The pulse response is an approximation to that described in GSM documents,
//     and is given by, p(t)~exp[-1.045(t/T)^2 - 0.218(t/T)^4]
// ############################# Private Method ###############################
LOGICAL FIRFilter:: transferForEDGEFilter ()
{
  int           n;
  float         time;
  double        *b, arg;
  
  b             = new double[_numberTaps];

//
// Set up the filter coefficients according to equation above
// The center time should be centered in the _numberTaps:
// Let T be normalized to 1.
//
  time          = -_symbolsPerSample*(_numberTaps-1)/2.;
  for(n=0; n<_numberTaps; n++)
  {
    arg         = time*time*(EDGE_COEFF1 + EDGE_COEFF2*time*time);
    b[n]        = exp(-arg);
    time        += _symbolsPerSample;
  }
//
// Store coefficients in polynomials:
//
  _bPolyObject->assign(_filterOrder, b);
  delete [] b;
  setPassBandGain();
  return YES;
}

// ############################# Private Method ###############################
// initInstanceVariables -- This routine initializes the object's instance variables
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
void FIRFilter::initInstanceVariables()
{
  _dataWindow           = NULL;
//
// Add these:
//
  setFIRFilterType(FIR_WINDOW);
  setRaisedCosineAlpha(DEFAULT_ALPHA);
  setPowerLawAlpha(DEFAULT_POWER_ALPHA);
  setSymbolsPerSample(DEFAULT_SYM_SAMP);
  setGMSKBT(DEFAULT_BT);

  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the FIRFilter class.
//
// Input:           pass:           Filter pass type: LOW_PASS, BAND_PASS, etc.
//                  centerFreq:     Filter center frequency in Hertz
//                  cutoffFreq:     Filter cutoff frequency in Hertz
//                  samplingFreq:   Digital sampling frequency in Hertz
//                  order           Filter order
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
FIRFilter::FIRFilter(int pass, double centerFreq, double cutoffFreq,
                                  double samplingFreq, int order)
              :AbstractFIR(pass, centerFreq, cutoffFreq, samplingFreq, order)
{
  initInstanceVariables();
  _dataWindow           = new DataWindow();
  
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the FIRFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
FIRFilter::~FIRFilter()
{
  delete _dataWindow;
  return;
}

// ############################# Public Method ###############################
// setWindowType -- Set a new windowing type for the data window.
//
// Input:   type:           E.g., Hamming
//          
// Output:                  None
//
// Notes:
// ############################# Public Method ###############################
void FIRFilter::setWindowType(int type)
{
  _dataWindow->setWindowType(type);
  return;
}

// ############################# Public Method ###############################
// filterDelay -- Returns the delay of the filter, assuming a linear phase
//                design.
//
// Input:               None
//
// Output:              Filter delay in samples
//          
// Notes:
// ############################# Public Method ###############################
float FIRFilter::filterDelay()
{
  return ((_numberTaps - 1.)/2.);
}

// ############################# Public Method ###############################
// windowType -- Returns the current window type setting.
//
// Input:                  None
//
// Output:   type:         E.g., Hamming
//          
// Notes:
// ############################# Public Method ###############################
int FIRFilter::windowType()
{
  return _dataWindow->windowType();
}

// ############################# Public Method ###############################
// findTransferFunction -- Find the filter's transfer function
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
//  1. The instanceVariables, _bPolyObject are modified.
// ############################# Public Method ###############################
void FIRFilter::findTransferFunction()
{
  switch(_firFilterType)
  {
    case FIR_WINDOW:
    case FIR_RAISED_COS:
    case FIR_ROOT_RAISED_COS:
    default:
      transferForFIRWindow();
      break;
    case FIR_POWER_LAW:
      transferForPowerLawFilter();
      break;
    case FIR_GMSK:
      transferForGMSKFilter();
      break;
    case FIR_EDGE:
      transferForEDGEFilter();
      break;
    case FIR_USER_DEFINED:                                      // Do nothing, user has entered coefficients
      break;
  }
  zeroOutTaps();
  return;
}
