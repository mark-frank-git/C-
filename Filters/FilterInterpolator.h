#ifndef _FILTER_INTERPOLATOR_H
#define _FILTER_INTERPOLATOR_H 1
/************************************************************************
 *                                                                      *
 * This subclass of object implements a class for interpolating an      *
 * input sequence using filter coefficients (e.g, sinc filter).         *
 * See Theory of Bandlimited interpolation by J.O. Smith.               *
 * http://www-ccrma.stanford.edu.                                       *
 *                                                                      *
 * File:FilterInterpolator.h                                            *
 *                                                                      *
 * The interpolator operates as follows:                                *
 * 1. For each output point center the filter at the output point.      *
 * 2. With the filter centered at the output point, 
 *                                                                      *
 ************************************************************************/

#include <stdio.h>

#ifndef LOGICAL
#define LOGICAL                 char
#endif

#ifndef YES
#define YES                     1
#define NO                      0
#endif

class FilterInterpolator
{
protected:
  int           _filterLength;                  // Size of the interpolating filter
  int           _outputSize;                    // Size of interpolator output

  LOGICAL       _filterOdd;                     // YES = odd length filter

  float         *_filterCoefficients;
  float         *_interpolatorOutput;           // output interpolated array

//
// Private methods:
//
  float         interpolateFilter(int index, float p);  // interpolate filter between samples
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  FilterInterpolator(const float *coefficient=NULL, int length=0);
  ~FilterInterpolator();

  
/**************************
 * Setting parameters:    *
 **************************/
  void          setFilterCoefficients(const float *coefficient, int length);

/**********************
 * Get parameters:    *
 **********************/
  int           outputSize()            {return _outputSize;}
  int           filterLength()          {return _filterLength;}
  const float   *filterCoefficients()   {return _filterCoefficients;}
  const float   *interpolatorOutput()   {return _interpolatorOutput;}

/**********************
 * Interpolating:    *
 **********************/
  void          interpolateArray(float *input, int inputSize, int startIndex, float outputStepSize,
                         float firstStep, int outputSize);
  void          filterArray(float *input, int inputSize, int startIndex, float outputStepSize,
                                            float firstStep, int outputSize);
};

#endif
