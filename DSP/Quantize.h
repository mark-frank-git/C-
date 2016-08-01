#ifndef _QUANTIZER_H
#define _QUANTIZER_H    1
/************************************************************************
 *                                                                      *
 * This subclass of object implements a quantizer.                      *
 *                                                                      *
 * File:Quantize.h                                                      *
 *                                                                      *
 *                                                                      *
 ************************************************************************/

#define DEFAULT_LEVELS  255
#define DEFAULT_SCALE   1.
class Quantize
{
private:
  int           _quantizerLevels;               // # of levels
  int           _quantizerType;                 // Y axis symmetry or X axis
  int           _outputLength;                  // Length of output array
  double        _quantizerFullScale;            // Full scale range
  double        _resolution;                    // 1 bit resolution
  double        _halfResolution;                // 1/2 bit resolution
  double        _maxLevel;                      // Max level for quantizer

  float         *_outputData;                   // Output quantized data
  
  void          initParameters();               // Initialize parameters for quantizing

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  Quantize(int levels=DEFAULT_LEVELS, double scale=DEFAULT_SCALE);
  ~Quantize();

/**********************
 * Set parameters:    *
 **********************/
  void          setQuantizerLevels(int levels);
  void          setQuantizerScale(double scale);

/**********************
 * Get parameters:    *
 **********************/
  int           quantizerLevels()               {return _quantizerLevels;}
  double        quantizerScale()                {return _quantizerFullScale;}
  double        quantizerStepSize();

/************************
 * Quantizing:          *
 ************************/
  float         *quantizeFloatArray(const float *inputArray, int length);
  float         quantizeFloat(float inputFloat);
  double        quantizeDouble(double inputDouble);
  

};
#endif
