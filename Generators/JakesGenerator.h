#ifndef _JAKESGENERATOR_H
#define _JAKESGENERATOR_H       1
/********************************************************************************
 * This class implements a fading simulator based on Jakes' model, resulting    *
 * in an envelope that is Rayleigh distributed. Provisions are made for the     *
 * user's _velocity.                                                            *
 * File: /User/frank/C++/Generators/JakesGenerator.h                            *
 *                                                                              *
 ********************************************************************************/
#if defined(WIN32)
#include <GNU/Complex.h>
#else
#include "Complex.h"
#endif


#ifndef LOGICAL
#define LOGICAL char
#endif

#define DEFAULT_NUMBER_SINE_WAVES       12                      // For Jakes model
#define DEFAULT_NUMBER_FADERS           12                      // number of independent faders
#define DEFAULT_DOPPLER                 200.                    // 200 Hz doppler
#define MIN_DOPPLER                     0.1
#define SAMPLE_STANDARD                 0.001
#define NUMBER_QUADRANTS                4

class JakesGenerator
{
private:
  int           _fadeNumber;                            // used for generating independent faders
  int           _numberSineWaves;                       // number of sine waves used in generator
  int           _oldNumberWaves;                        // old number for mallocs/deletes
  int           _numberFaders;                          // # of independent faders
  long          _randomSeed;                            // random number generator seed
  long          _numberSimulatePoints;                  // # of simulated points
  long          _oldNumberPoints;                       // Used to keep track of malloc'ing points

  float         _dopplerFreq;                           // Doppler frequency in Hz
  float         *_iData;                                // simulated magnitude data for I channel
  float         *_qData;                                // simulated magnitude data for Q channel
  float         *_magnitudeData;                        // sqrt(i**2 + q**2) given in dB

  double        _samplingTime;                          // Time between samples given in seconds
  double        _normalizingGain;                       // Gain to normalize output energy
  double        *_cosThetas;                            // Cosine wave arguments
  double        *_sinThetas;                            // Sine wave arguments
  double        *_cosPhiInits;                          // Initial cosine phases
  double        *_sinPhiInits;                          // Initial sine phases
  double        *_omegaMCosAlphaDeltaT;                 // omegaM*cos(alphaNK)*delta_time
  double        *_omegaMSinAlphaDeltaT;                 // omegaM*sin(alphaNK)*delta_time
  double        _tableScale;                            // Look up table scale factor

//
// private methods:
//
  double        urand();                                // Returns uniform random # for initial phases
  void          updateParameters();                     // Initialize the generator parameters
  

//
// public methods:
//
public:
  JakesGenerator(int fadeNumber, int numberSineWaves= DEFAULT_NUMBER_SINE_WAVES,
                 int numberFaders=DEFAULT_NUMBER_FADERS,
                 float dopplerFreq=DEFAULT_DOPPLER,
                 float sampleTime = SAMPLE_STANDARD);
                                                        // Constructor
  ~JakesGenerator();                                    // Destructor
//
// The following functions set parameters:
//
  void  initGenerator();
  void  setFadeNumber(int walsh);                       // Sets new Walsh number
  void  setNumberSineWaves(int number);                 // Sets new # of sine waves
  void  setNumberFaders(int number);                    // Sets new # of independent faders
  void  setSamplingTime(float sampleTime);              // Set new sampling time in seconds
  void  setDopplerFreq(float freq);                     // Set new doppler in Hertz

//
// The following functions generate simulation data.
// NOTE: generateIandQData() must be called before retrieving any data
//
  void          generateIAndQData(long numberSamples);  // Run the generator for I and Q data
  Complex       complexFading();                        // Get I and Q fading as complex
  const         float *calculateFadeMagnitude();        // Return the magnitude from previously calculated
                                                        // I and Q data
//
// The following functions get parameters:
//
  long   getNumberSimulatePoints()              {return _numberSimulatePoints;}
  float getSamplingTime()                       {return _samplingTime;}
  const float *getIData()                       {return _iData;}
  const float *getQData()                       {return _qData;}
  const float *getMagnitudeData()               {return _magnitudeData;}

};

#endif