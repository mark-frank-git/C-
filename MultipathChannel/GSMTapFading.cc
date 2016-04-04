/************************************************************************
 *                                                                      *
 * This class implements the fading for one of the taps in multipath    *
 * channel model.  It uses a sum of sinusoids, or Jakes Generator.      *
 *                                                                      *
 * File:GSMTapFading.cc                                                 *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1.  12/19/01  - Started.                                            *
 ************************************************************************/
#include "GSMTapFading.h"                                       // Object prototypes
#include "environment_types.h"
#include <GNU/Complex.h>
#include <Generators/JakesGenerator.h>
#include <C_Libraries/constants.h>
#include <math.h>
#include <stdio.h>

#define NO_FADE_MAG     0.707107

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#define SGN(a)    ( ((a)>=0.) ? 1 : -1 )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the GSMTapFading class.
//
// Input:       filterType:     Doppler filter type, Classical, Rician, etc.
//              doppler:        Doppler frequency in Hertz
//              sampling:       Sampling frequency in Hertz
//              noiseSeed:      Random number generator seed
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
GSMTapFading::GSMTapFading(int filterType, int fadeNumber, int numberSineWaves, int numberFaders,
                           float doppler, float sampling)

{
  _jakesGen     = new JakesGenerator(fadeNumber);
  setFilterType(filterType);
  setFadeNumber(fadeNumber);
  setNumberSineWaves(numberSineWaves);
  setNumberFaders(numberFaders);
  setDopplerFrequency(doppler);
  setSamplingFrequency(sampling);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the GSMTapFading class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
GSMTapFading::~GSMTapFading()
{
  delete _jakesGen;
  return;
}

// ############################# Public Method ###############################
// resetTap -- Resets the tap by initializing filters, and re-running an
//             initial number of samples.
// Input:               None
//          
// Output:              None
// Notes:
// ############################# Public Method ###############################
void GSMTapFading::resetTap()
{
  _jakesGen->initGenerator();
  return;
}

// ############################# Public Method ###############################
// setFilterType -- Sets a new filter type (Classical or Rician).
// Input:       filterType:     New filter type.
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void GSMTapFading::setFilterType(int filterType)
{
  _dopplerFilterType    = filterType;
  return;
}

// ############################# Public Method ###############################
// setFadeNumber -- Sets a new fade number.
// Input:       number:         New fade number for independent faders.
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void GSMTapFading:: setFadeNumber(int number)
{
  _jakesGen->setFadeNumber(number);
  return;
}

// ############################# Public Method ###############################
// setNumberSineWaves -- Sets a new number of sine waves.
// Input:       number:         New number of sinusoids in Jakes generator.
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void GSMTapFading:: setNumberSineWaves (int number)
{
  _jakesGen->setNumberSineWaves(number);
  return;
}

// ############################# Public Method ###############################
// setNumberFaders -- Sets a new number of faders.
// Input:       faders:         New number of faders.
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void GSMTapFading:: setNumberFaders(int faders)
{
  _jakesGen->setNumberFaders(faders);
  return;
}

// ############################# Public Method ###############################
// setDopplerFrequency -- Sets a new Doppler frequency.
// Input:       sampling:       New Doppler frequency in Hertz.
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void GSMTapFading:: setDopplerFrequency(float doppler)
{
  _jakesGen->setDopplerFreq(doppler);
  return;
}

// ############################# Public Method ###############################
// setSamplingFrequency -- Sets a new sampling frequency.
// Input:       sampling:       New sampling frequency in Hertz.
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void GSMTapFading:: setSamplingFrequency(float sampling)
{
  _samplingFrequency    = MAX(MIN_SAMPLING_FREQ, sampling);
  _jakesGen->setSamplingTime(1./_samplingFrequency);
  return;
}

// ############################# Public Method ###############################
// tapFading -- Returns a new value for the tap fading.
// Input:                       None
//          
// Output:                      tap fading
// Notes:
// ############################# Public Method ###############################
Complex  GSMTapFading::tapFading()
{
  Complex       fading_output;
  switch(_dopplerFilterType)
  {
    case CLASSICAL_TYPE:
      return _jakesGen->complexFading();
    case RICIAN_TYPE:
      return _jakesGen->complexFading();
    case GAUSSIAN_TYPE:
      return _jakesGen->complexFading();
    case NO_FADING_TYPE:
    default:
      fading_output     = Complex(NO_FADE_MAG, NO_FADE_MAG);
      return fading_output;
  }
}
