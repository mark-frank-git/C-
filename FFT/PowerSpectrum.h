#ifndef _POWER_SPECTRUM_H
#define _POWER_SPECTRUM_H   1
/************************************************************************
 *                                                                      *
 * This subclass of object implements an object for computing PSDs,     *
 * autocorrelations, etc.                                               *
 *                                                                      *
 * File: PowerSpectrum.h                                                *
 *                                                                      *
 ************************************************************************/
#include    "FFT.h"

#define DB_THRESHOLD            15.             // Find mean index from this threshold

class   FFT;                                    // Class prototype
class   DataWindow;

    
//! Mean types enum
/*! The mean_type_enum categories the different types of mean calculations used in meanIndexOf */
enum mean_type_enum
{
  MEAN_FROM_SCM,                        //!< Calc mean index from spectral center mass
  MEAN_FROM_THRESHOLD                   //!< Calc from threshold dB down
};


//! Avg types enum
/*! The avg_type_enum categories the different types of averaging when calculating PSD over multiple blocks
    in the avgPowerSpectrum() functions */
enum avg_type_enum
{
  PSD_AVG,                              //!< average PSDs
  PSD_MIN,                              //!< take min of PSDs
  PSD_MAX                               //!< take max of PSDs
};

//! Calibrate types enum
/*! The calibrate_type_enum categories the different ways of scaling the PSD
    in the avgPowerSpectrum() function */
enum calibrate_type_enum
{
  CALIBRATE_FOR_CWS,                    //!< scale to get sine wave amplitudes correct
  CALIBRATE_FOR_NOISE                   //!< scale to get noise levels correct
};

//! PowerSpectrum Class
/*!
   The PowerSpectrum class offers the calculations of PSDs and autocorrelations.
*/
class PowerSpectrum
{
private:
  DataWindow    *_dataWindow;                   //!< Data windowing object
  FFT           *_fft;                          //!< FFT object
  int           _numberSpectrumPoints;          //!< Number of points in PSD
  int           _numberFloatPoints;             //!< Number of points in _outputFloat
  int           _numberAutoPoints;              //!< Number of points in autocorrelation
  int           _numberPhasePoints;             //!< Number of points in phase spectrum
  int           _oldSpectrumPoints;             //!< Old number of spectrum points
  int           _oldPhasePoints;                //!< Old number of phase points
  int           _oldAutoPoints;                 //!< Old number of autocorrelation points
  int           _meanIndexType;                 //!< Method to calc mean index
  int           _calibrateMode;                 //!< Calibrate PSDs for noise or CWs

  float         _twoSidedBandwidth;             //!< Calc'd in meanIndexOf
  float         *_outputSpectrum;               //!< results of FFT
  float         *_outputPhase;                  //!< Phase spectrum output
  float         *_autocorrelation;              //!< Autocorrelation function 
  float         *_outputFloat;                  //!< results of complexToFloat
//
// Private functions
//

public:

//
// Public functions
//
//! Class constructor, the sampling frequency is used for frequency axis scaling
  PowerSpectrum(float sampling=1000.);              // Class constructor
//!< Class destructor
 ~PowerSpectrum();                                  // Class destructor

/**********************
 * Set parameters:    *
 **********************/
//! setSamplingFrequency() allows the setting of a new sampling frequency
/**
  * a normal member taking one argument with no return value
  * @param frequency the new sampling frequency in Hertz
  * @see frequencyPoints()
  * @see twoSidedFrequency()
  * @return void
  */
  void setSamplingFrequency(float frequency)    {_fft->setSamplingFrequency(frequency);  return;}

//! setWindowType() allows the setting of a new data windowing type. The types are defined
//! in DataWindow.h
/**
  * a normal member taking one argument with no return value
  * @param type the new data windowing type
  * @return void
  */
  void setWindowType(int type);

//! setMeanIndexType() allows the setting of a new mean type.
/**
  * a normal member taking one argument with no return value
  * @param type the new mean index type
  * @see meanIndexOf()
  * @return void
  */
  void setMeanIndexType(int type)               {_meanIndexType  = type; return;}

//! setCalibrateMode() allows the setting of a PSD scaling calibration mode.
/**
  * a normal member taking one argument with no return value
  * @param mode the type of PSD scaling
  * @see avgPowerSpectrum()
  * @return void
  */
  void setCalibrateMode(int mode)               {_calibrateMode  = mode; return;}

//*********************************************************************************************
// The following functions return parameters
//
// windowType:           Returns the time series windowing function type
// frequencyPoints:      Returns an array of spectral frequency points
// twoSidedFrequency:    Similar to above, but for a two-sided spectrum
// twoSidedBandwidth:    Two sided power spectrum bandwidth. NOTE: must first call
//                          meanIndexOf().
// deltaFrequency:      Returns the spacing between frequency points
// isPowerOfTwo:        Returns whether or not the input is a power of two
//*****************************************************************************************
  int           windowType();
  int           calibrateMode()                 {return _calibrateMode;                 }
  float         *frequencyPoints()              {return _fft->frequencyPoints();        }
  float         *twoSidedFrequency()            {return _fft->twoSidedFrequency();      }
  float         twoSidedBandwidth()             {return _twoSidedBandwidth;             }
  float         deltaFrequency(int fftSize)     {return _fft->deltaFrequency(fftSize);  }
  LOGICAL       isPowerOfTwo(int n)             {return _fft->isPowerOfTwo(n);          }

//*********************************************************************************************
// The following functions calculate system parameters
//
// resolutionBandwidth:         Returns the resolution bandwidth in Hertz given the # of samples
// frequencyFromIndex:          Returns the frequency corresponding to input array index
// numberDataSamples:           Returns the number of data samples needed to achieve res bandwidth
// nextFFTSize:                 Returns the next largest FFT size given the # of samples
// indexFromFrequency:          Returns FFT array index corresponding to input frequency
//*****************************************************************************************
  float         resolutionBandwidth(int numberSamples);
  float         frequencyFromIndex(float index);
  int           numberDataSamples(float resolutionBandwidth);
  int           nextFFTSize(int numberSamples);
  int           indexFromFrequency(float frequency);

//*********************************************************************************************
// The following functions calculate power spectrums
//
// powerSpectrum:               Non-windowed power spectrum, no averaging
// psdFromAuto:                 Power spectrum from input autocorrelation
// avgPowerSpectrum:            Windowed and averaged power spectrum
// avgPowerSpectrum:            Windowed and averaged from Complex data
// avgPowerSpectrum:            Windowed and averaged power spectrum using Chirp Z
// avgPowerSpectrumFromReal:    Windowed and averaged power spectrum from real input data
// avgPowerSpectrumFromReal:    Windowed and averaged power spectrum from real input data using chirp z
//*****************************************************************************************
  float *powerSpectrum(float *input, int numberPoints);
  float *psdFromAuto(float *input, int numberPoints);
  float *avgPowerSpectrum(const float *realInput, const float *imagInput, int fftSize, int overlapSize, int numberAvgs,
                          int type=MAGNITUDE_IN_WATTS, int windowSize=0, LOGICAL initialize=YES,
                          int averageType=PSD_AVG);
  float *avgPowerSpectrum(const Complex *complexInput, int fftSize, int overlapSize, int numberAvgs,
                          int type=MAGNITUDE_IN_WATTS, int windowSize=0, LOGICAL initialize=YES,
                          int averageType=PSD_AVG);
  float *avgPowerSpectrum(const float *realInput, const float *imagInput, int inputSize, int overlapSize,       // Chirp Z
                          int numberAvgs, int numberOutputPts,
                          float startFreq, float endFreq, int type=MAGNITUDE_IN_WATTS,
                          LOGICAL initialize=YES, int averageType=PSD_AVG);
  float *avgPowerSpectrumFromReal(const float *realInput, int fftSize, int overlapSize, int numberAvgs,
                                  int type=MAGNITUDE_IN_WATTS, int windowSize=0, LOGICAL initialize=YES,
                                  int averageType=PSD_AVG);
  float *avgPowerSpectrumFromReal(const float *realInput, int inputSize, int overlapSize, int numberAvgs,       // Chirp Z
                          int numberOutputPts,
                          float startFreq, float endFreq, int type=MAGNITUDE_IN_WATTS,
                          LOGICAL initialize=YES, int averageType=PSD_AVG);
  
//*********************************************************************************************
// The following functions calculate phase spectrums
//
// phaseSpectrum:               Non-windowed phase spectrum, no averaging
//*****************************************************************************************
  float *phaseSpectrum(const float *realInput, const float *imagInput, int numberPoints);
  
//*********************************************************************************************
// The following functions calculate complex spectrums
//
// complexSpectrum:             Non-windowed complex spectrum, no averaging
//*****************************************************************************************
  void  complexSpectrum(float *realInput, float *imagInput, int numberPoints, LOGICAL windowData=0);

//*********************************************************************************************
// The following functions calculate autocorrelations
//
// autoFromPSD:                 Autocorrelation from input power spectrum
// autocorrelationOf:           Time domain autocorrelation
//*****************************************************************************************
  float *autoFromPSD(const float *inputSpectrum, int numberPoints, int type, float *scale);
  float *autocorrelationOf(const float *realInput, const float *imagInput, int numberPoints, int windowLength,
                           int lagLength, int type, float *scale);

//*********************************************************************************************
// The following functions perform operations on autocorrelations and power spectrums:
//
// meanIndexOf:                 Finds the middle of a power spectrum array
// maxBinFrequency:             Finds maximum bin frequency (see the PeakDetect class for better algorithms)
// smoothSpectrum:              Smooths the input spectrum by using an averaging window
// shiftRight:                  Circular shifts an array (PSD) right for 2 sided display
// shiftRight:                  Circular shift for complex array
// magnitudeToDB:               Converts watts to dB
// normalizeAuto:               Normalizes autocorrelation function
// complexToFloat               Converts a complex spectrum to floating point
// oneSidedToTwoSided:          Converts the previously calculated spectrum from one-sided to two-sided
//*****************************************************************************************
  float meanIndexOf(const float *inputSpectrum, int numberPoints);
  float maxBinFrequency(const float *inputSpectrum, int numberPoints, float flo, float fhi);
  void  smoothSpectrum(float *inputSpectrum, int numberPoints, int windowSize);
  void  shiftRight(float *inputArray, int shift, int arraySize);
  void  shiftRight(Complex *inputArray, int shift, int arraySize);
  void  magnitudeToDB(float *inputSpectrum, int numberPoints);
  void  normalizeAuto(float *autoc, int numberPoints, float scale);
  float *complexToFloat(const Complex *inputArray, int arraySize, int conversionType);
  float *oneSidedToTwoSided();
};

#endif
