#ifndef _SPECTRAL_DISPLAY_H
#define _SPECTRAL_DISPLAY_H   1
/************************************************************************
 *                                                                      *
 * This subclass of object implements an object for calculating         *
 * spans, overlaps, etc. for spectral display using the Chirp Z         *
 * transform.                                                           *
 *                                                                      *
 * File: SpectralDisplay.h                                              *
 *                                                                      *
 ************************************************************************/

#include <Filters/DataWindow.h>

#define DEFAULT_SAMPLING        22.1184e6       // 66.3552/3
#define DEFAULT_WINDOW          BLACKMAN_HARRIS4
#define DEFAULT_ACQ_SPAN        10.e6
#define DEFAULT_SPAN            120.e6
#define DEFAULT_RBW             1.e6
#define DEFAULT_CENTER          900.e6
#define DEFAULT_STEP            614.4e3

#define MIN_SPAN                10.e3
#define MAX_SPAN                2.5e9
#define MIN_CENTER              25.e6
#define MAX_CENTER              1.2e9
#define MIN_RBW                 100.
#define MAX_RBW                 10.e6
#define MAX_SAMPLING            100.e6
#define MIN_SAMPLING            100.
#define DEFAULT_DISP_POINTS     501

class SpectralDisplay
{
private:
  DataWindow    *_dataWindow;                   // Data windowing object

  int           _displayPoints;                 // Desired # of display points

  float         _spanRange;                     // Total desired span in Hertz
  float         _centerFrequency;               // Center frequency of span in Hertz
  float         _acquisitionSpan;               // Max span in Hertz for a single acquisition
  float         _samplingFrequency;             // sampling frequency in Hertz for single acquisition
  float         _resolutionBW;                  // desired resolution BW in Hertz
  float         _synthesizerStepSize;           // Size of the RF LO step size in Hertz
//
// Private functions
//

public:
  SpectralDisplay(float span=DEFAULT_SPAN, float center=DEFAULT_CENTER,
                  float acqSpan=DEFAULT_ACQ_SPAN, float sampling=DEFAULT_SPAN,
                  float rbw=DEFAULT_RBW,
                  int windowType=DEFAULT_WINDOW);       // Class constructor
 ~SpectralDisplay();                                    // Class destructor

/**********************
 * Set parameters:    *
 **********************/
  void  setTotalSpan(float totalSpan);
  void  setCenterFrequency(float center);

  void  setAcquisitionSpan(float acqSpan);
  void  setSamplingFrequency(float samplingFreq);
  void  setResolutionBW(float rbw);
  void  setWindowType(int type);
  void  setSynthesizerStepSize(float step);

//*********************************************************************************************
// The following functions return previously set parameters
//
// windowType:          Returns the specified window function type
// spanRange:           Returns the specified span in Hertz
// centerFrequency      Returns the span center frequency in Hertz
// acquisitionSpan:     Returns the specified max span for a single acq in Hertz.
// samplingFrequency:   Returns the specified sampling frequency for a single acq in Hertz.
// resolutionBW:        Returns the specified resolution BW in Hertz.
// synthesizerStepSize  Returns the specified synthesizer step in Hertz.
//*****************************************************************************************
  int           windowType()                    {return _dataWindow->windowType();      }
  float         spanRange()                     {return _spanRange;                     }
  float         centerFrequency()               {return _centerFrequency;               }
  float         acquisitionSpan()               {return _acquisitionSpan;               }
  float         samplingFrequency()             {return _samplingFrequency;             }
  float         resolutionBW()                  {return _resolutionBW;                  }
  float         synthesizerStepSize()           {return _synthesizerStepSize;           }

//*********************************************************************************************
// The following functions return calculated parameters
//
// numberAcquisitionsFor: Returns # of acquisitions given an input overlap # spectral points
// dataSize:            Returns size of windowed data to achieve RBW
// synthesizerSteps     Returns the number of synthesizer steps between acquisitions
// frequencySpacing:    Returns the frequency spacing in Hz
// overlapRange:        Returns overlap range in frequency given synth step size, etc. 
//*****************************************************************************************
  int           numberAcquisitions();
  int           dataSize();
  int           synthesizerSteps();
  float         frequencySpacing();
  float         overlapRange();
  float         acqStartFrequencyFor(int acqNumber);
  float         acqCenterFrequencyFor(int acqNumber);
  float         acqEndFrequencyFor(int acqNumber);
  int           acqPointsFor(int acqNumber, float startFreq=-1.);
  int           acqOverlapPointsFor(int acqNumber=0, float lastEnd=-1.);

};


#endif
