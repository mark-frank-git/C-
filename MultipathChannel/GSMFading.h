#ifndef _GSM_FADING_H
#define _GSM_FADING_H    1
/************************************************************************
 *                                                                      *
 * This class implements a multipath fading channel based on the COST   *
 * 207 model.                                                           *
 *                                                                      *
 * File:GSMFading.h                                                     *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1.  10/23/02  - Started.                                            *
 *  2.  01/28/03  - Fixed normalizing of tap gains.                     *
 *  3.  02/04/03  - Added _fadeSeed to get independent faders, e.g.,    *
 *                  for a cochannel interferer.                         *
 ************************************************************************/
#include "environment_types.h"

class   Complex;
class   ComplexBuffer;
class   GSMTapFading;

#define DEFAULT_SPEED           13.411          // 13.411 m/s ~ 30 MPH
#define DEFAULT_CARRIER_FREQ    900.e6          // 900 MHz
#define DEFAULT_SAMP_FREQ       (18.*1625e3/6)  // 18 * GSM symbol rate

#define MIN_VELOCITY            0.01            // Minimum settable velocity
#define MIN_CARRIER             100.
#define MIN_SAMPLING            100.

class GSMFading
{
private:
  int                   _fadeSeed;                      // Seed to get indpendent GSMFading objects
  int                   _environmentType;               // Fading environment
  int                   _oldEnvironment;                // Saves on calls to allocate

  float                 _mobileVelocity;                // Mobile velocity in meters/second
  float                 _carrierFrequency;              // Carrier frequency in Hertz
  float                 _samplingFrequency;             // Tapped delay line sampling frequency in Hz
  float                 _oldVelocity, _oldCarrier, _oldSampling;        // Saves on calls to allocate
  
  int                   _numberPaths;                   // Number of paths to used in the delay line
  int                   _oldNumberPaths;                // Used for delete and alloc
  int                   *_tapLocations;                 // Array of tap locations corresponding to paths
  float                 *_tapGains;                     // Scalar gain associated with each tap
  GSMTapFading          **_tapFading;                   // Array of tap fading corresponding to paths
  ComplexBuffer         *_delayLine;                    // Tapped delay line

//
// Private functions:
//
  void                  allocateObjects();
  void                  deleteObjects();                // Delete allocated objects
// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  GSMFading(int fadeSeed=0, int type=RA_TYPE, float velocity=DEFAULT_SPEED,
            float carrier=DEFAULT_CARRIER_FREQ, float sampling=DEFAULT_SAMP_FREQ);
  ~GSMFading();
  
/********************************
 * Resetting the channel.       *
 ********************************/
  void          resetChannel();

/**********************
 * Set parameters:    *
 **********************/
  void          setFadeSeed(int seed);
  void          setEnvironmentType(int type);
  void          setMobileVelocity(float velocity);
  void          setCarrierFrequency(float carrier);
  void          setSamplingFrequency(float sampling);

/**********************
 * Get parameters:    *
 **********************/
  int           fadeSeed()                      {return _fadeSeed;              }
  int           environmentType()               {return _environmentType;       }
  float         mobileVelocity()                {return _mobileVelocity;        }
  float         carrierFrequency()              {return _carrierFrequency;      }
  float         samplingFrequency()             {return _samplingFrequency;     }
  float         dopplerFrequency();

/************************
 * Processing input     *
 * through the channel: *
 ************************/
  Complex       processInput(const Complex &channelInput);              // process new sample through channel

};
#endif
