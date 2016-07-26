#ifndef _TWOPOINTINTERP_H
#define _TWOPOINTINTERP_H 1
/************************************************************************
 *                                                                      *
 * This is used for finding the peak location from two points.  The     *
 * algorithm assumes the underlying function is triangle shaped, like   *
 * a PN autocorrelation.                                                *
 *                                                                      *
 * File: /User/frank/C++/Measurements/TwoPointInterp.h                  *
 *                                                                      *
 ************************************************************************/

#ifndef LOGICAL
#define LOGICAL    char
#endif
#define MIN_TRIANGLE_WIDTH      0.1

class TwoPointInterp
{
protected:
  float         _triangleWidth;                 // Width of triangle base in samples

//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  TwoPointInterp();
  virtual ~TwoPointInterp();

/*******************************
 * These methods set parameters*
 *******************************/
  void  setTriangleWidth(float width);

/****************************************
 * These methods find the peak position *
 * from a set of autocorrelation data.  *
 ****************************************/
  float peakPosition(float *autoData, int low_index, int high_index);

};

#endif