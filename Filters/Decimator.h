#ifndef _DECIMATOR_H
#define _DECIMATOR_H    1
/************************************************************************
 *                                                                      *
 * This class implements a decimator (no filtering).                    *
 *                                                                      *
 * File:Decimator.h                                                     *
 *                                                                      *
 ************************************************************************/

#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif

class   Complex;

class Decimator
{
protected:
  int                   _decimateFactor;                // Factor to decimate data
  int                   _oldDecimatePoints;             // Old size of _decimateOutput;
  int                   _oldComplexPoints;              // Old size of _complexOutput

  float                 *_decimateOutput;               // Output of decimator
  Complex               *_complexOutput;                // Complex output of decimator
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  Decimator(int decimation=1);
  virtual  ~Decimator();

/**************************
 * Setting parameters:    *
 **************************/
  void  setDecimateFactor(int factor)           {_decimateFactor = factor; return;}

/************************
 * Get parameters:      *
 ***********************/
  int           getDecimateFactor()             {return _decimateFactor;}

/************************
 * Decimating and       *
 * interpolating.       *
 ************************/
  float         *decimateInput(float *input, int numberPts);            // Decimate a float array
  Complex       *decimateComplexData(Complex *input, int numberPts);    // Decimate a complex array


};

#endif
