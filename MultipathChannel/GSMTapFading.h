#ifndef _GSM_TAP_FADING_H
#define _GSM_TAP_FADING_H    1
/************************************************************************
 *                                                                      *
 * This class implements the fading for one of the taps in multipath    *
 * channel model.  It uses a sum of sinusoids, or Jakes Generator.      *
 *                                                                      *
 * File:GSMTapFading.h                                                  *
 *                                                                      *
 ************************************************************************/
class   Complex;
class   JakesGenerator;

#define MIN_SAMPLING_FREQ       100.            // Minimum allowed sampling frequency

class GSMTapFading
{
private:
  int                   _dopplerFilterType;             // Type of doppler filter, classical, etc.
  
  float                 _dopplerFrequency;              // Doppler frequency in Hertz
  float                 _samplingFrequency;             // output sampling frequency in Hz

  JakesGenerator        *_jakesGen;                     // Jakes generator
//
// Private functions:
//
// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  GSMTapFading(int filterType, int fadeNumber, int numberSineWaves, int numberFaders, float doppler,
               float sampleFrequency);
  ~GSMTapFading();


/********************************
 * Resetting the tap.           *
 ********************************/
  void          resetTap();

/**********************
 * Set parameters:    *
 **********************/
  void          setFilterType(int filterType);
  void          setFadeNumber(int number);
  void          setNumberSineWaves(int number);
  void          setNumberFaders(int numberFaders);
  void          setDopplerFrequency(float doppler);
  void          setSamplingFrequency(float sampling);


/**********************
 * Get parameters:    *
 **********************/
  float         dopplerFrequency()              {return _dopplerFrequency;      }
  float         samplingFrequency()             {return _samplingFrequency;     }

/************************
 * Returns the current  *
 * fading for the tap.  *
 ************************/
  Complex       tapFading();

};
#endif
