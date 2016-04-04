#ifndef _DELTASIGMAPREDICTOR_H
#define _DELTASIGMAPREDICTOR_H    1
/************************************************************************
 *                                                                      *
 * This class calculates and solves the linear prediction equations     *
 * for a bandpass delta sigma converter.                                *
 *                                                                      *
 * File:DeltaSigmaPredictor.h                                           *
 *                                                                      *
 * The autocorrelation or linear prediction equations take the form:    *
 *                                                                      *
 *    Rxx * A = R                                                       *
 * where                                                                *
 * Rxx = |rxx(0) rxx(1) .. rxx(p-1)|                                    *
 *       |rxx(1) rxx(0) .. rxx(p-2)|                                    *
 *       |                         |                                    *
 *       |rxx(p-1)..        rxx(0) |                                    *
 * is the symmetric autocorrelation matrix.                             *
 * A = [a(0) a(1) .. a(p-1)]T   is the coefficient vector.              *
 * R = [rxx(1) rxx(2) .. rxx(p)] is the one step predictor correlation  *
 *                               vector.                                *
 * The autocorrelation vector is given by:                              *
 * rxx(n) = exp(-j*IF*n)*sin(BW*n)/(2*PI*n).                            *
 * Where the IF and BW frequencies are given in rad/s.  See the January *
 * 2000 issue of the IEEE Signal Processing Magazine for more details.  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/04/00  - Started.                                             *
 ************************************************************************/

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX_PRED_COEFFS 256

class LinearPredictor;                          // class prototype
class Complex;

class DeltaSigmaPredictor
{
protected:
  int           _numberCoeffs;                  // Number of predictor coefficients
  int           _oldNumber;                     // Old number of coefficients

  double        _ifFrequency;                   // IF frequency in rads [0, PI]
  double        _bwFrequency;                   // Two sided bandwidth in rads [0, PI]
  double        _lambda;                        // constant added to diagonal of autocorrelation matrix

  Complex       *_autocorrelationVector;        // Autocorrelation vector
  Complex       *_predictorCoefficients;        // Calculated predictor coefficients

  LinearPredictor       *_linearPredictor;      // Predictor coefficient calculator

//
// Private functions:
//
  void  calculateAutocorrelationVector();       // Calculate the autocorrelation vector

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  DeltaSigmaPredictor(int number=0, double ifFreq = PI/2., double bwFreq = PI/8., double lambda=0.0);
  ~DeltaSigmaPredictor();

/**********************
 * Set parameters:    *
 **********************/
  void  setNumberCoeffs(int number);
  void  setIFFrequency(double ifFreq);
  void  setBWFrequency(double bwFreq);
  void  setLambda(double lambda);

/**********************
 * Get parameters:    *
 **********************/
  int           numberCoeffs()          {return _numberCoeffs;}
  double        ifFrequency()           {return _ifFrequency;}
  double        bwFrequency()           {return _bwFrequency;}
  double        lambda()                {return _lambda;}
  const Complex *autocorrelationVector(){return _autocorrelationVector;}

/**********************
 * Getting Outputs:   *
 **********************/
  const Complex *calculatePredictorCoefficients();

};
#endif
