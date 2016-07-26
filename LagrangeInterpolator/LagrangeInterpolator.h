#ifndef _LAGRANGE_INTERPOLATOR_H
#define _LAGRANGE_INTERPOLATOR_H 1
/************************************************************************
 *                                                                      *
 * This subclass of object implements a class for interpolating an      *
 * input sequence using Lagrange's formula.  There are also methods     *
 * for real time processing of input data samples.                      *
 * See Abramowitz & Stegun, pp. 878-879.                                *
 *                                                                      *
 * File:LagrangeInterpolator.h                                          *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/23/00  - Started.                                             *
 ************************************************************************/

#define MIN_INTERP_POINTS       2               // Linear interpolation
#define MAX_INTERP_POINTS       6
#define DEFAULT_TIME_STEP       1.


class CLagrangeInterpolator
{
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  explicit CLagrangeInterpolator(int numberPoints=3);
  ~CLagrangeInterpolator();

/********************************
 * Initializing interpolator    *
 ********************************/
  void          reset();
  
/**************************
 * Setting parameters:    *
 **************************/
  void          setNumberInterpPoints(int points);

/**********************
 * Get parameters:    *
 **********************/
  int           numberInterpPoints()    {return m_iNumberInterpPoints;}

/********************************
 * The following methods        *
 * are used for interpolating   *
 * arrays of data:              *
 ********************************/
  void          interpolateArray(const float *input, int inputSize, int startIndex,
                         double outputStepSize, double firstStep, double *output, int outputSize);   // Interpolate an evenly spaced array
  void          interpolateUnevenArray(const float *inputYData, const float *inputXData,
                         const float *outputXData, int inputSize, double *output, int outputSize); // Interpolate an unevenly spaced array

protected:
  int           m_iNumberInterpPoints;            // # of points to use in interpolation

//
// Private methods:
//
  double        interpolate(double f_m2, double f_m1, double f_0, double f_p1, double f_p2, double p);
  double        interpolateDouble(double *f, double p);

};

#endif

