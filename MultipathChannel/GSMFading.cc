/************************************************************************
 *                                                                      *
 * This class implements a multipath fading channel based on the COST   *
 * 207 model.                                                           *
 *                                                                      *
 * File:GSMFading.cc                                                    *
 *                                                                      *
 *                                                                      *
 ************************************************************************/
#include "GSMFading.h"                                  // Object prototypes
#include "GSMTapFading.h"
#include "cost207_channels.h"
#include <GNU/Complex.h>
#include <Buffers/ComplexBuffer.h>
#include <Generators/JakesGenerator.h>
#include <C_Libraries/constants.h>
#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0. ? (a) : (-a))
#define SGN(a)    ( ((a)>=0.) ? 1 : -1 )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )
#define EPS1    0.1

#define NORMALIZE_TAPS  1

// ############################# Private Method ###############################
// allocateObjects -- Allocates the objects based on environmentType.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
void GSMFading::allocateObjects ()
{
  int           i, offset, filter_type, number_taps;
  int           fade_number;
  float         doppler_freq, max_delay_time;
  double        total_power;

  _numberPaths          = number_paths[_environmentType];
  if(_oldNumberPaths < _numberPaths)
  {
    deleteObjects();                            // Delete old objects
    _oldNumberPaths             = _numberPaths;
    _tapLocations               = new int[_numberPaths];
    _tapFading                  = new GSMTapFading * [_numberPaths];
    _tapGains                   = new float[_numberPaths];
  }
  else                                          // Delete _tapFading for news[] below
  {
    for(i=0; i<_oldNumberPaths; i++)
    {
      delete _tapFading[i];
      _tapFading[i]             = NULL;
    }
  }
//
// Get the offset into the environment type arrays:
//
  offset                        = 0;
  for(i=0; i<_environmentType; i++)
  {
    offset                      += number_paths[i];
  }
//
// Allocate the tap fading, the type of
// Doppler filter needs to be selected.
//
  doppler_freq                  = dopplerFrequency();
  for(i=0; i<_numberPaths; i++)
  {
    fade_number                 = _fadeSeed + i;
    filter_type                 = filter_types[offset+i];
    _tapFading[i]               = new GSMTapFading(filter_type, fade_number, DEFAULT_NUMBER_SINE_WAVES, _numberPaths,
                                        doppler_freq, _samplingFrequency);

  }
//
// Allocate the tapped delay line:
//
  max_delay_time                = MICROSECONDS*tap_times[offset+_numberPaths-1];
  number_taps                   = ROUND(max_delay_time*_samplingFrequency);
  delete _delayLine;
  _delayLine                    = new ComplexBuffer(number_taps, 1);
  _delayLine->resetBuffer();
//
// Find the tap locations and gains:
//
  total_power                   = 0.;
  for(i=0; i<_numberPaths; i++)
  {
    _tapLocations[i]            = ROUND((MICROSECONDS*tap_times[offset+i]*_samplingFrequency));
    _tapGains[i]                = tap_powers[offset+i];                 // These are given in dB
    _tapGains[i]                = pow(10., _tapGains[i]/10.);
    total_power                 += _tapGains[i];
  }
#if NORMALIZE_TAPS
  if(total_power != 0.)
  {
    for(i=0; i<_numberPaths; i++)
    {
      _tapGains[i]              = sqrt(_tapGains[i]/total_power);       // Convert power to voltage gain
    }
  }
#else
  for(i=0; i<_numberPaths; i++)
  {
      _tapGains[i]              = sqrt(_tapGains[i]);                   // Convert power to voltage gain
  }
#endif
  
  return;
}

// ############################# Private Method ###############################
// deleteObjects -- Deletes all of the allocated arrays and objects.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
void GSMFading::deleteObjects()
{
  int   i;
  for(i=0; i<_oldNumberPaths; i++)
  {
    delete _tapFading[i];
  }
  delete [] _tapFading;
  delete [] _tapGains;
  delete [] _tapLocations;
  delete _delayLine;
  _delayLine    = NULL;                 // In case of double deletes above

  return;
}


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the GSMFading class.
//
// Input:       fadeSeed:       Seed for independent faders
//              type:           Environment type, RA, etc..
//              velocity:       Mobile velocity in meters/s.
//              carrier:        Carrier frequency in Hertz
//              sampling:       Delay line sampling in Hertz.
//
// Output:                          None
//
// Notes:
// 1. Reset channel is not called from here, so user needs to explicitly call it
//    before running samples.
// ############################# Class Constructor ###############################
GSMFading::GSMFading(int fadeSeed, int type, float velocity,
            float carrier, float sampling)

{
  _numberPaths          = _oldNumberPaths = 0;
  _tapLocations         = NULL;
  _tapFading            = NULL;
  _tapGains             = NULL;
  _delayLine            = NULL;
  _oldEnvironment       = -1;
  _oldVelocity          = _oldCarrier   = _oldSampling  = 0.;
  setFadeSeed(fadeSeed);
  setEnvironmentType(type);
  setMobileVelocity(velocity);
  setCarrierFrequency(carrier);
  setSamplingFrequency(sampling);
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the GSMFading class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
GSMFading::~GSMFading()
{
  deleteObjects();
  return;
}

// ############################# Public Method ###############################
// resetChannel -- Resets the channel's conditions to initial values
// Input:               None
//          
// Output:              None
// Notes:
// ############################# Public Method ###############################
void GSMFading::resetChannel()
{
  int           i;
  Complex       filter_input;
//
// Allocate objects with new parameter settings:
//
  if( _oldEnvironment!=_environmentType || (ABS((_oldVelocity - _mobileVelocity)) > EPS1) ||
      (ABS((_oldCarrier - _carrierFrequency)) > EPS1) || (ABS((_oldSampling - _samplingFrequency)) > EPS1) )
  {
    _oldVelocity        = _mobileVelocity;
    _oldCarrier         = _carrierFrequency;
    _oldSampling        = _samplingFrequency;
    _oldEnvironment     = _environmentType;
    allocateObjects();    
  }
//
// Reset the fading generator:
//
  for(i=0; i<_numberPaths; i++)
  {
    _tapFading[i]->resetTap();
  }
  return;
}

// ############################# Public Method ###############################
// setFadeSeed -- Sets a seed for independent faders.
// Input:       seed:           New fading seed
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void GSMFading::setFadeSeed(int seed)
{
  _fadeSeed             = MAX(0, seed);
  return;
}

// ############################# Public Method ###############################
// setEnvironmentType -- Sets a new type of environment.
// Input:       type:           New environment type, RA, TU, etc.
//          
// Output:                      None
// Notes:
// 1. The fading simulator is not re-loaded with the new environment type until
//    resetChannel() is called.
// ############################# Public Method ###############################
void GSMFading::setEnvironmentType(int type)
{
  _environmentType      = MAX(0, type);
  _environmentType      = MIN(_environmentType, (NUMBER_ENV_TYPES));
  return;
}

// ############################# Public Method ###############################
// setMobileVelocity -- Sets a new mobile velocity in m/s.
// Input:       velocity:       New velocity in meters/second
//          
// Output:                      None
// Notes:
// 1. The fading simulator is not re-loaded with the new velocity until
//    resetChannel() is called.
// ############################# Public Method ###############################
void GSMFading:: setMobileVelocity(float velocity)
{
  _mobileVelocity       = MAX(MIN_VELOCITY, velocity);
}

// ############################# Public Method ###############################
// setCarrierFrequency -- Sets a new carrier frequency in Hertz.
// Input:       carrier:        New carrier frequency in Hertz
//          
// Output:                      None
// Notes:
// 1. The fading simulator is not re-loaded with the new frequency until
//    resetChannel() is called.
// ############################# Public Method ###############################
void GSMFading::setCarrierFrequency(float carrier)
{
  _carrierFrequency     = MAX(MIN_CARRIER, carrier);
}

// ############################# Public Method ###############################
// setSamplingFrequency -- Sets a new sampling frequency in Hertz.
// Input:       carrier:        New sampling frequency in Hertz
//          
// Output:                      None
// Notes:
// 1. The fading simulator is not re-loaded with the new frequency until
//    resetChannel() is called.
// ############################# Public Method ###############################
void GSMFading::setSamplingFrequency(float sampling)
{
  _samplingFrequency    = MAX(MIN_SAMPLING, sampling);
}

// ############################# Public Method ###############################
// dopplerFrequency -- returns the doppler frequency given the current settings
//                     for the mobile velocity and the carrier frequency.
//
// Input:                       None
//
// Output:                      Doppler frequency in Hertz
//
// Notes:
// ############################# Public Method ###############################
float GSMFading::dopplerFrequency()
{
  float doppler;

  doppler       = _carrierFrequency*_mobileVelocity/VELOCITY_LIGHT;
  return doppler;
}

// ############################# Public Method ###############################
// processInput -- Sends input data through the channel, returns output.
// Input:       inputData:      New data to be filtered.
//          
// Output:                      Channel output
// Notes:
// ############################# Public Method ###############################
Complex GSMFading::processInput(const Complex &channelInput)
{
  int           i;
  Complex       channel_output, filter_input, buffer_input[1];
  double        tap_sum;
  const Complex *delay_output;
//
// Get the channel output by summing outputs from delay line:
//
  channel_output        = Complex(0., 0.);
  tap_sum               = 0.;
  for(i=0; i<_numberPaths; i++)
  {
    if(_tapLocations[i] == 0)
    {
      channel_output    += _tapGains[i]*_tapFading[i]->tapFading()*channelInput;
    }
    else
    {
      delay_output      = _delayLine->readOldSampleAt(_tapLocations[i]-1);
      channel_output    += _tapGains[i]*_tapFading[i]->tapFading()*delay_output[0];
    }
    tap_sum             += _tapGains[i];
  }
//
// Add new input to delay line:
//
  buffer_input[0]       = channelInput;
  _delayLine->writeToBuffer(buffer_input);
  return channel_output;
}
