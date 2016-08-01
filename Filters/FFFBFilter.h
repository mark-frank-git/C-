#ifndef _FF_FB_FILTER_H
#define _FF_FB_FILTER_H 1
/************************************************************************
 *                                                                      *
 * This class implements a generalized feedforward feedback filter      *
 * for use in delta sigma modulators.  See 5.6.4 in Delta-Sigma Data    *
 * Converters: Theory, Design, and Simulation.                          *
 *                                                                      *
 * File: FFFBFilter.h                                                   *
 *                                                                      *
 * The structure of the filter is as follows:                           *
 *                                                                      *
 * x -----------------------------                                      *
 *    |            |              |                                     *
 *    b0           b1             bn                                    *
 *    |     1      |       1      |     1                               *
 *    +  -------   +    -------   +  -------   -----> Quant() ----> y   *
 *   -|  z - p0   -|    z - p1   -|  z - pn                   |         *
 *    |            |              |                           |         *
 *    a0           a1             an                          |         *
 *    |            |              |                           |         *
 *    ---------------------------------------------------------         *
 *                                                                      *
 * where the pi's set the transfer function poles, the bi's set the     *
 * signal transfer function zeros, and the ci's set the noise transfer  *
 * function zeros. Also, note that the coefficients and poles are       *
 * complex in order to effect a bandpass filter.                        *
 *                                                                      *
 ************************************************************************/
#include <stdio.h>
#include <DataIO/TestData.h>

#define MAX_FILTER_ORDER        1024
#ifndef LOGICAL
#define LOGICAL char
#endif
#ifndef YES
#define YES     1
#define NO      0
#endif

class   Complex;

class FFFBFilter
{
protected:
  int           _filterOrder;                   // Order of the filter (# of poles)
  int           _numberFileData;                // Number of read in file coefficients
  int           _oldNumberData;                 // Old number of read in coeffs
  int           _shiftSize;                     // Number of shift stages allocated
  int           _quantLevels;                   // Number of quantization levels for output

  double        _quantScale;                    // Quantization saturation level
  double        _feedbackGain;                  // The gain in the feedback path

  Complex       *_aCoeffs;                      // Feedback coefficients
  Complex       *_bCoeffs;                      // Feedforward coefficients
  Complex       *_filterPoles;                  // Filter poles
  Complex       *_inputShift;                   // complex shift register
  Complex       *_backShift;                    // complex shift register
  Complex       *_fileData;                     // Read in file coefficients/poles/zeros

  TestData      *_fileReader;                   // Reads coefficients from file

//
// Private methods:
//
  LOGICAL       readDataFile(const char *fileName, int fileType);
  void  checkShiftSize(int size);                       // Check for size of shifts
  inline Complex filterNextDouble(double input);        // Local function for filtering
  inline Complex filterNextComplex(Complex input);      // Local function for filtering
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  FFFBFilter(Complex *aCoeffs=NULL, Complex *bCoeffs=NULL, Complex *poles=NULL, int order=2);
  virtual  ~FFFBFilter();

/**************************
 * Setting parameters:    *
 **************************/
  virtual       void    setQuantLevels(int levels)      {_quantLevels   = levels;       return;}
  virtual       void    setQuantScale(double scale)     {_quantScale    = scale;        return;}
  virtual       void    setFeedbackGain(double gain)    {_feedbackGain  = gain;         return;}
  virtual       void    setFilterACoeffs(const Complex *a);
  virtual       void    setFilterBCoeffs(const Complex *b);
  virtual       void    setFilterPoles(const Complex *poles);
  virtual       LOGICAL readACoeffsFromFile(const char *fileName, int fileType);
  virtual       LOGICAL readBCoeffsFromFile(const char *fileName, int fileType);
  virtual       LOGICAL readPolesFromFile(const char *fileName, int fileType);
  virtual       void    setFilterOrder(int order);
  virtual       void    zeroOutTaps();

/**********************
 * Get parameters:    *
 **********************/
  int           filterOrder()           {return _filterOrder;}
  int           quantLevels()           {return _quantLevels;}
  double        quantScale()            {return _quantScale;}
  double        feedbackGain()          {return _feedbackGain;}
  const Complex *aCoeffs()              {return _aCoeffs;}
  const Complex *bCoeffs()              {return _bCoeffs;}
  const Complex *poles()                {return _filterPoles;}
        
/****************************************
 * Filtering float and Complex data     *
 * assuming an IIR structure.  A sub-   *
 * class may wish to override these     *
 * with more efficient implementations. *
 ****************************************/
  virtual       Complex filterFloatData(float input);                   // Filter float data
  virtual       Complex filterDoubleData(double input);                 // Filter double data
  virtual       LOGICAL filterComplexArray(Complex *x, int numberPts);  // Filter complex array
  virtual       Complex filterComplexData(Complex &input);              // Filter complex data

};

#endif
