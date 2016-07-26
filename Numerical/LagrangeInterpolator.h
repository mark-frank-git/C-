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
 ************************************************************************/

#ifndef LOGICAL
#define LOGICAL                 char
#endif

#ifndef YES
#define YES                     1
#define NO                      0
#endif
#define MIN_INTERP_POINTS       2               // Linear interpolation
#define MAX_INTERP_POINTS       6
#define DEFAULT_TIME_STEP       1.

class   DoubleBuffer;

class LagrangeInterpolator
{
protected:
  int           _numberInterpPoints;            // # of points to use in interpolation
  int           _outputSize;                    // Size of interpolator output
  int           _inputSample;                   // Input sample counter for real time processing
  int           _outputSample;                  // Output sample counter for real time processing

  float         *_interpolatorOutput;           // output interpolated array
  float         _inputTimeStep;                 // Input sample time (Ts) in seconds for real time processing
  float         _outputTimeStep;                // Output sample time (Ts) in seconds for real time processing

  DoubleBuffer  *_inputBuffer;                  // Input buffer for real time processing

//
// Private methods:
//
  double        interpolate(float f_m2, float f_m1, float f_0, float f_p1, float f_p2, float p);
  double        interpolateDouble(double *f, double p);
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  LagrangeInterpolator(int numberPoints=3);
  ~LagrangeInterpolator();

/********************************
 * Initializing interpolator    *
 ********************************/
  void          reset();
  
/**************************
 * Setting parameters:    *
 **************************/
  void          setNumberInterpPoints(int points);
  void          setInputTimeStep(float time);
  void          setOutputTimeStep(float time);

/**********************
 * Get parameters:    *
 **********************/
  int           numberInterpPoints()    {return _numberInterpPoints;}
  int           outputSize()            {return _outputSize;}
  const float   *interpolatorOutput()   {return _interpolatorOutput;}

/********************************
 * The following methods        *
 * are used for interpolating   *
 * arrays of data:              *
 ********************************/
  void          interpolateArray(const float *input, int inputSize, int startIndex,
                         float outputStepSize, float firstStep, int outputSize);   // Interpolate an evenly spaced array
  void          interpolateArray(const float *inputYData, const float *inputXData,
                         const float *outputXData, int inputSize, int outputSize); // Interpolate an unevenly spaced array


/********************************
 * The following methods        *
 * are used for interpolating   *
 * real time data.              *
 ********************************/
  int           numberOfInputPointsNeeded();                                    // Returns # of input points needed for next output
  void          bufferNextInput(double inputData);
  double        getNextOutput();

  
};

#endif
