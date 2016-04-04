#ifndef _SINC_INTERPOLATOR_H
#define _SINC_INTERPOLATOR_H 1
/************************************************************************
 *                                                                      *
 * This subclass of object implements a class for interpolating an      *
 * input sequence using a windowed sinc pulse.                          *
 *                                                                      *
 * File:SincInterpolator.h                                              *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 09/03/00  - Started.                                             *
 ************************************************************************/

#ifndef LOGICAL
#define LOGICAL                 char
#endif

#ifndef YES
#define YES                     1
#define NO                      0
#endif
#define DEFAULT_TIME_STEP       1.
#define DEFAULT_OVER_SAMPLING   10
#define DEFAULT_TRANS_WIDTH     0.05
#define MIN_TRANSITION          0.0001
#define MAX_TRANSITION          0.9
#define DEFAULT_CROSSINGS       3
#define DEFAULT_WINDOW          0               // Rectangular

class   DoubleBuffer;                           // Class prototypes
class   FIRFilter;

class SincInterpolator
{
protected:
  int           _numberCrossings;               // # of (one sided) sinc crossings to use
  int           _oversamplingFactor;            // Oversampling factor used in calculating sinc function
  int           _windowType;                    // Type of window function to use in sinc function
  int           _inputSample;                   // Input sample counter for real time processing
  int           _outputSample;                  // Output sample counter for real time processing
  int           _firstOutputSample;             // Index of first non-zero output sample
  int           _numberFilterTaps;              // Calculated number of sinc taps
  int           _numberBufferedSamples;         // Calculated number of input samples needed

  float         _inputTimeStep;                 // Input sample time (Ts) in seconds for real time processing
  float         _outputTimeStep;                // Output sample time (Ts) in seconds for real time processing
  float         _transitionWidth;               // Width of filter transition width as a fraction [0,1] of fs/2

  double        *_sincFunction;                 // The calculated sinc function
  DoubleBuffer  *_inputBuffer;                  // Input buffer for real time processing
  FIRFilter     *_firFilter;                    // FIR filter for generating sinc function

//
// Private methods:
//
  void          calculateSincFunction();        // Calculates new sinc function
  double        interpolateDouble(double p);
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  SincInterpolator(int numberCrossing=DEFAULT_CROSSINGS);
  ~SincInterpolator();

/********************************
 * Initializing interpolator    *
 ********************************/
  void          reset();
  
/**************************
 * Setting parameters:    *
 **************************/
  void          setNumberCrossings(int crossings);
  void          setOversamplingFactor(int factor);
  void          setWindowType(int type);
  
  void          setInputTimeStep(float time);
  void          setOutputTimeStep(float time);
  void          setTransitionWidth(float width);

/**********************
 * Get parameters:    *
 **********************/
  int           numberCrossings()       {return _numberCrossings;}
  int           overSamplingFactor()    {return _oversamplingFactor;}
  int           windowType()            {return _windowType;}
  int           numberFilterTaps()      {return _numberFilterTaps;}

  float         inputTimeStep()         {return _inputTimeStep;}
  float         outputTimeStep()        {return _outputTimeStep;}
  float         transitionWidth()       {return _transitionWidth;}


/********************************
 * The following methods        *
 * are used for interpolating   *
 * real time data.              *
 ********************************/
  int           numberOfInputPointsNeeded();                            // Returns # of input points needed for next output
  void          bufferNextInput(double inputData);
  double        getNextOutput();

  
};

#endif
