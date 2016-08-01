#ifndef _QUANTIZER_SNR_H
#define _QUANTIZER_SNR_H 1
/********************************************************************************
 *                                                                              *
 * This class calculates the output SNR from an n level quantizer.              *
 *                                                                              *
 * File: /User/frank/C++/SNRCalc/QuantizerSNR.h                                 *
 *                                                                              *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif

class   NoiseCorrelation;                       // Class prototype
//
// The following are the correlator types:
//
#define DEFAULT_STEP            1.              // Quantizer step size
#define DEFAULT_SNR             -20.            // Default SNR in dB
#define DEFAULT_SAMPLING        20.46e6         // Default sampling freq in Hz

class QuantizerSNR
{
protected:
  NoiseCorrelation      *_noiseAuto;            // Calculates noise autocorrelation function
  int                   _quantizerLevels;       // Number of levels in the quantizer

  double                _quantizerStepSize;     // Quantizer spacing, delta
  double                _inputSNR;              // Input signal-to-noise ratio (not in dB)
//
// Private functions
//
  double        lfnh(double h, double k, double rho);   // cumulative binomial distribution function
  double        mean2Level();                   // Special case of 1 bit quantizer
  double        varianceIJ2Level(double tau);   // Special case of 1 bit quantizer
  double        varianceII2Level();             // Special case of 1 bit quantizer
public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  QuantizerSNR(int levels=4, double step=DEFAULT_STEP, double inputSNR=DEFAULT_SNR);    // Constructor
  virtual ~QuantizerSNR();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setNumberLevels(int levels)             {_quantizerLevels       = levels;       return;}
  void  setQuantizerStepSize(double step)       {_quantizerStepSize     = step;         return;}
  void  setInputSNR(double snr);
  void  setFilterType(int type);
  void  setFilterOrder(int order);
  void  setFilterCutoff(double cutoff);
  void  setNoiseDensity(double density);
  void  setCenterFrequency(double freq);

/*******************************
 * These methods get parameters:*
 *******************************/
  int           numberLevels()                  {return _quantizerLevels;}
  double        quantizerStepSize()             {return _quantizerStepSize;}
  double        inputSNR()                      {return _inputSNR;}
  int           filterType();
  int           filterOrder();
  double        filterCutoff();
  double        noiseDensity();
  double        noiseSigma();
  double        centerFrequency();

/********************************
 * These methods calculate and  *
 * return the output SNR, etc.  *
 * for the given settings.      *
 ********************************/
  double        meanValue();
  double        varianceIJ(double tau);
  double        varianceII();
};

#endif
