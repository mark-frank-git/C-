#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H 1
/************************************************************************
 *                                                                      *
 * This class implements an interpolator (no filtering).                *
 *                                                                      *
 * File:Interpolator.h                                                  *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 06/08/99  - Started.                                             *
 ************************************************************************/

#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif


class Interpolator
{
protected:
  int                   _interpolateFactor;             // Factor to interpolate data
  int                   _oldInterpolatePoints;          // Old size of _interpolateOutput;

  float                 *_interpolateOutput;            // Output of decimator
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  Interpolator(int interpolation=1);
  virtual  ~Interpolator();

/**************************
 * Setting parameters:    *
 **************************/
  void  setInterpolateFactor(int factor)                {_interpolateFactor = factor; return;}

/************************
 * Get parameters:      *
 ***********************/
  int   getInterpolateFactor()          {return _interpolateFactor;}

/************************
 * Decimating and       *
 * interpolating.       *
 ************************/
  float         *interpolateInput(float * input, int numberPts);        // Interpolate a float array


};

#endif
