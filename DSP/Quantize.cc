/************************************************************************
 *                                                                      *
 * This subclass of object implements a quantizer.                      *
 *                                                                      *
 * File:Quantize.cc                                                     *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 09/18/97  - Started                                              *
 *  2. 11/02/00  - Add quantizerStepSize().                             *
 ************************************************************************/

#include <stdio.h>
#include "Quantize.h"

#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)          ((a) >= 0 ? (a) : (-a))
#define SGN(a)          ( ((a)>=0.) ? 1 : -1 )

/************************************************
 * Quantizer types (depend on number of levels) *
 ************************************************/
#define Y_SYMMETRIC_QUANT               0                       // Symmetric about y axis
#define X_SYMMETRIC_QUANT               1                       // Symmetric about x axis
#define ONE_BIT_QUANT                   2                       // Signum

#define OFFSET_QUANTIZER                0.5
#ifndef EPS
#define EPS     1.e-10
#endif

// ############################# Private Method ###############################
// initParameters --  Initialize the quantizer's parameters
//
// Input:               None
//          
// Output:              None
//
// Notes:
// ############################# Private Method ###############################
void Quantize::initParameters()
{
  if(_quantizerLevels > 2)
  {
    if(_quantizerLevels%2 != 0)                         // Odd number of levels, symmetric about Y axis
    {
      _quantizerType    = Y_SYMMETRIC_QUANT;
      _halfResolution   = _quantizerFullScale/(_quantizerLevels-1);
      _resolution       = 2.*_halfResolution;
      _maxLevel         = _quantizerFullScale;
    }
    else                                                // Even number of levels, symmetric about X axis        
    {
      _quantizerType    = X_SYMMETRIC_QUANT;
      _halfResolution   = _quantizerFullScale/_quantizerLevels;
      _resolution       = 2.*_halfResolution;
      _maxLevel         = _quantizerFullScale - _halfResolution;
    }
  }
  else
  {
     if(_quantizerLevels == 2)
       _quantizerType   = ONE_BIT_QUANT;
     else
     {
       _quantizerType   = Y_SYMMETRIC_QUANT;    // 1 level not really valid?
       _halfResolution  = _resolution    = 1.;
     }
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Quantize class.
//
// Input:       levels:         Number of quantizer levels
//              scale:          Positive full scale value
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
Quantize::Quantize(int levels, double scale)
{
  _quantizerLevels      = levels;
  _quantizerFullScale   = scale;
  _outputData           = NULL;
  _outputLength         = 0;
  initParameters();

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the Quantize class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
Quantize::~Quantize()
{
  delete [] _outputData;
  return;
}


// ############################# Public Method ###############################
// setQuantizerLevels -- Set new number of quantizer levels
//
// Input:       levels:         # of quantizer levels
//          
// Output:                      None
//
// ############################# Public Method ###############################
void Quantize::setQuantizerLevels(int levels)
{
  _quantizerLevels      = MAX(2, levels);
  initParameters();
  return;
}

// ############################# Public Method ###############################
// setQuantizerScale -- Set new quantizer scale
//
// Input:       scale:          new positive full scale
//          
// Output:                      None
//
// ############################# Public Method ###############################
void Quantize::setQuantizerScale(double scale)
{
  _quantizerFullScale   = MAX(EPS, scale);
  initParameters();
  return;
}

// ############################# Public Method ###############################
// quantizerStepSize -- Returns the step size of the quantizer
//
// Input:                       None
//          
// Output:                      Quantizer step size
//
// ############################# Public Method ###############################
double Quantize::quantizerStepSize()
{
  double        denominator, step_size;

  if(_quantizerLevels%2 == 0)           // Even quantizer
    denominator = _quantizerLevels/2.;
  else
    denominator = (_quantizerLevels-1)/2.;
  step_size     = 0.;
  if(denominator > 0.)
    step_size   = _quantizerFullScale/denominator;
  return step_size;
}

// ############################# Public Method ###############################
// quantizeFloatArray -- Quantize an array of floats
//
// Input:       inputArray:             float array to quantizer
//              length:                 length of inputArray
//          
// Output:                              output quantized array
//
// ############################# Public Method ###############################
float *Quantize::quantizeFloatArray(const float *inputArray, int length)
{
  int                   i;
  double                abs_input, multiplier;

//
// Initialize output array:
//
  if(_outputLength < length)
  {
    delete []   _outputData;
    length              = MAX(1, length);
    _outputLength       = length;
    _outputData         = new float[length];
  }
//
// Peform the quantization:
//
  switch(_quantizerType)
  {
    case Y_SYMMETRIC_QUANT:                             // Symmetrical about Y axis
      for(i=0; i<length; i++)                           // Loop over length of input data
      {
        abs_input               = ABS(inputArray[i]);
        multiplier              = (int)(abs_input/_resolution + OFFSET_QUANTIZER);
        _outputData[i]          = MIN((multiplier*_resolution), _maxLevel);
        if(inputArray[i] < 0.)
          _outputData[i]        = -_outputData[i];
      }
      break;
    case X_SYMMETRIC_QUANT:                             // Symmetrical about X axis
      for(i=0; i<length; i++)                           // Loop over length of input data
      {
        abs_input               = ABS(inputArray[i]);
        if(abs_input < _resolution)
          _outputData[i]        = _halfResolution;
        else
        {
          multiplier            = (int)(abs_input/_resolution);
          multiplier            += OFFSET_QUANTIZER;
          _outputData[i]        = MIN((multiplier*_resolution), _maxLevel);
        }
        if(inputArray[i] < 0.)
          _outputData[i]        = -_outputData[i];
      }
      break;
    case ONE_BIT_QUANT:
      for(i=0; i<length; i++)                           // Loop over length of input data
      {
       _outputData[i]           = SGN(inputArray[i]);
      }
      break;
  }

  return _outputData;
}

// ############################# Public Method ###############################
// quantizeFloat -- Quantize a single input float value
//
// Input:       inputFloat:             float data to be quantized
//          
// Output:                              output quantized value
//
// ############################# Public Method ###############################
float Quantize::quantizeFloat(float inputFloat)
{
  float         abs_input, multiplier, output_data;
//
// Peform the quantization:
//
  switch(_quantizerType)
  {
    default:
    case Y_SYMMETRIC_QUANT:                             // Symmetrical about Y axis
      abs_input         = ABS(inputFloat);
      multiplier        = (int)(abs_input/_resolution + OFFSET_QUANTIZER);
      output_data       = MIN((multiplier*_resolution), _maxLevel);
      if(inputFloat < 0.)
        output_data     = -output_data;
      break;
    case X_SYMMETRIC_QUANT:                             // Symmetrical about X axis
      abs_input         = ABS(inputFloat);
      if(abs_input < _resolution)
        output_data     = _halfResolution;
      else
      {
          multiplier    = (int)(abs_input/_resolution);
          multiplier    += OFFSET_QUANTIZER;
          output_data   = MIN((multiplier*_resolution), _maxLevel);
      }
      if(inputFloat < 0.)
          output_data   = -output_data;
      break;
    case ONE_BIT_QUANT:
      output_data       = SGN(inputFloat);
      break;
  }

  return output_data;
}



// ############################# Public Method ###############################
// quantizeDouble -- Quantize a single input double value
//
// Input:       inputDouble:            double variable to be quantized
//          
// Output:                              output quantized value
//
// ############################# Public Method ###############################
double Quantize::quantizeDouble(double inputDouble)
{
  double                abs_input, multiplier, output_data;
//
// Peform the quantization:
//
  switch(_quantizerType)
  {
    default:
    case Y_SYMMETRIC_QUANT:                             // Symmetrical about Y axis
      abs_input         = ABS(inputDouble);
      multiplier        = (int)(abs_input/_resolution + OFFSET_QUANTIZER);
      output_data       = MIN((multiplier*_resolution), _maxLevel);
      if(inputDouble < 0.)
        output_data     = -output_data;
      break;
    case X_SYMMETRIC_QUANT:                             // Symmetrical about X axis
      abs_input         = ABS(inputDouble);
      if(abs_input < _resolution)
        output_data     = _halfResolution;
      else
      {
          multiplier    = (int)(abs_input/_resolution);
          multiplier    += OFFSET_QUANTIZER;
          output_data   = MIN((multiplier*_resolution), _maxLevel);
      }
      if(inputDouble < 0.)
          output_data   = -output_data;
      break;
    case ONE_BIT_QUANT:
      output_data       = SGN(inputDouble);
      break;
  }

  return output_data;
}


