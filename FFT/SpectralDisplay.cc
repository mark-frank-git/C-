/************************************************************************
 *                                                                      *
 * This subclass of object implements an object for calculating         *
 * spans, overlaps, etc. for spectral display using the Chirp Z         *
 * transform.                                                           *
 *                                                                      *
 * File: SpectralDisplay.cc                                             *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 09/07/99  - Started                                              *
 ************************************************************************/

/**************************
 * Include files:         *
 **************************/
#include        "SpectralDisplay.h"
#include        <math.h>
#include        <stdlib.h>


#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

// ############################# Class Constructor #################################
// SpectralDisplay -- Constructor for the SpectralDisplay class
//
// Input:       sampling:   Sampling Frequency in Hertz
//
// Output:                  None
//
// ############################# Class Constructor #################################
SpectralDisplay::SpectralDisplay(float span, float center, float acqSpan, float sampling,
                  float rbw, int windowType)

{
  _dataWindow            = new DataWindow();

  _displayPoints        = DEFAULT_DISP_POINTS;

  setSynthesizerStepSize(DEFAULT_STEP);
  setTotalSpan(span);
  setCenterFrequency(center);

  setAcquisitionSpan(acqSpan);
  setSamplingFrequency(sampling);
  setResolutionBW(rbw);
  setWindowType(windowType);
  return;
}

// ############################# Class Destructor #################################
// ~SpectralDisplay -- Destructor for the SpectralDisplay class
//
// Input:       sampling:   Sampling Frequency in Hertz
//
// Output:                  None
//
// ############################# Class Destructor #################################
SpectralDisplay::~SpectralDisplay()
{
  delete    _dataWindow;
  return;
}

// ############################# Public Method ###############################
// setTotalSpan -- Sets a new total span range in Hertz.
//
// Input:       span:           New span in Hertz
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void SpectralDisplay::setTotalSpan(float span)
{
  _spanRange    = MAX(MIN_SPAN, span);
  _spanRange    = MIN(MAX_SPAN, _spanRange);
  return;
}

// ############################# Public Method ###############################
// setTotalSpan -- Sets a new total span range in Hertz.
//
// Input:       span:           New span in Hertz
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void SpectralDisplay::setCenterFrequency(float center)
{
  _centerFrequency      = MAX(MIN_CENTER, center);
  _centerFrequency      = MIN(MAX_CENTER, _centerFrequency);
  return;
}

// ############################# Public Method ###############################
// setAcquisitionSpan -- Sets a new maximum span for a single acquisition.
//
// Input:       span:           New acquistion span in Hertz
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void SpectralDisplay::setAcquisitionSpan(float span)
{
  _acquisitionSpan      = MAX(MIN_SPAN, span);
  _acquisitionSpan      = MIN(MAX_SPAN, _acquisitionSpan);
  return;
}

// ############################# Public Method ###############################
// setSamplingFrequency -- Sets a new sampling frequency in Hertz.
//
// Input:       sampling:       New sampling frequency in Hertz
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void SpectralDisplay::setSamplingFrequency(float sampling)
{
  _samplingFrequency    = MAX(MIN_SAMPLING, sampling);
  _samplingFrequency    = MIN(MAX_SAMPLING, _samplingFrequency);
  return;
}

// ############################# Public Method ###############################
// setResolutionBW -- Sets a new RBW in Hertz.
//
// Input:       rbw:            New RBW in Hertz
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void SpectralDisplay::setResolutionBW(float rbw)
{
  _resolutionBW         = MAX(MIN_RBW, rbw);
  _resolutionBW         = MIN(MAX_RBW, _resolutionBW);
  return;
}

// ############################# Public Method ###############################
// setWindowType -- Set a new windowing type for the data window.
//
// Input:       type:           E.g., Hamming
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void SpectralDisplay::setWindowType(int type)
{
  _dataWindow->setWindowType(type);
  return;
}

// ############################# Public Method ###############################
// setSynthesizerStepSize -- Set a new synthesizer step size
//
// Input:       step:           Synthesizer step size in Hertz
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void SpectralDisplay::setSynthesizerStepSize(float step)
{
  _synthesizerStepSize  = MAX(0., step);
  return;
}

// ############################# Public Method ###############################
// numberAcquisitions -- Calculates and returns the # of acquisitions needed
//                       to achieve specified span.
//
// Input:                       None
//          
// Output:                      # of acquisitions
//
// Notes:
// ############################# Public Method ###############################
int SpectralDisplay::numberAcquisitions()
{
  int   acquisitions;
  float overlap_size, float_acqs;

  overlap_size  = overlapRange();
  if(_acquisitionSpan > overlap_size)
  {
    float_acqs          = (_spanRange-overlap_size)/(_acquisitionSpan -overlap_size);
    acquisitions        = ROUND(float_acqs+0.5);                                // Round up
  }
  else
    acquisitions        = 0;
  return acquisitions;
}

// ############################# Public Method ###############################
// dataSize -- Calculates and returns the data size needed to achieve specified RBW
//
// Input:                       None
//          
// Output:                      Size of windowed data
//
// Notes:
// 1. The # of data samples depends on the window type
// ############################# Public Method ###############################
int SpectralDisplay::dataSize()
{
  int   number_samples;
//
// Note: _resolutionBW has been constrained > 0.
//
  number_samples  = ROUND(_dataWindow->bandwidth3dB()*_samplingFrequency/_resolutionBW);
  return number_samples;
}

// ############################# Public Method ###############################
// synthesizerSteps -- Calculates and returns the number of synthesizer steps
//                      between acquisitions
//
// Input:                       None
//          
// Output:                      number of synthesizer steps, an integer in order
//                              to re-use TFR compensation
//
// Notes:
// ############################# Public Method ###############################
int SpectralDisplay::synthesizerSteps()
{
  int   number_steps;
  float float_steps;
//
// Round down required steps, so acquistions overlap:
//
  float_steps   = _acquisitionSpan/_synthesizerStepSize;
  number_steps  = ROUND((float_steps-0.5));
  return number_steps;
}

// ############################# Public Method ###############################
// frequencySpacing -- Calculates and returns the frequency spacing between
//                     Chirp Z output points.
//
// Input:                       None
//          
// Output:                      Chirp Z output spacing in Hertz
//
// Notes:
// ############################# Public Method ###############################
float SpectralDisplay::frequencySpacing()
{
  float freq_spacing;

  freq_spacing  = _spanRange/(_displayPoints-1);

  return freq_spacing;
}

// ############################# Public Method ###############################
// overlapRange -- Calculates and returns the frequency range of the overlap
//                 region
//
// Input:                       None
//          
// Output:                      overlap range in Hz
//
// Notes:
// ############################# Public Method ###############################
float SpectralDisplay::overlapRange()
{
  int   m;
  float overlap_range;

//
// get # of freq steps between acquisitions, and calculate overlap range from this
//
  m             = synthesizerSteps();
  overlap_range = _acquisitionSpan - m*_synthesizerStepSize;
  return overlap_range;
}

// ############################# Public Method ###############################
// acqStartFrequencyFor -- Calculates and returns the start frequency for the
//                         given acquisition
//
// Input:       acqNumber:      The acquisition # starting from 0
//              
//          
// Output:                      start frequency in Hz
//
// Notes:
// ############################# Public Method ###############################
float SpectralDisplay::acqStartFrequencyFor(int acqNumber)
{
  int   overlap;
  float last_end, spacing, start_freq;

  if(acqNumber < 1)
    return (_centerFrequency - _spanRange/2.);
//
// First get the overlap points with previous acq, and
// the end frequency of previous acq
//
  last_end      = acqEndFrequencyFor(acqNumber-1);
  overlap       = acqOverlapPointsFor(acqNumber, last_end);
  spacing       = frequencySpacing();
//
// Now, calculate start frequency from above:
//
  if(overlap < 2)
    start_freq  = last_end;
  else
    start_freq  = last_end - (overlap-1)*spacing;
  return start_freq;
}

// ############################# Public Method ###############################
// acqCenterFrequencyFor -- Returns the center frequency of given acq
//
// Input:       acqNumber:      The acquisition # starting from 0
//              
//          
// Output:                      center frequency in Hz
//
// Notes:
// ############################# Public Method ###############################
float SpectralDisplay::acqCenterFrequencyFor(int acqNumber)
{
  int   steps;
  float f_low, acq_center;

  f_low         = _centerFrequency - _spanRange/2.;
  steps         = synthesizerSteps();
  acq_center    = f_low + _acquisitionSpan/2. + acqNumber*steps*_synthesizerStepSize;
  return acq_center;
}

// ############################# Public Method ###############################
// acqEndFrequencyFor -- Calculates and returns the end frequency for the
//                         given acquisition
//
// Input:       acqNumber:      The acquisition # starting from 0
//              
//          
// Output:                      end frequency in Hz
//
// Notes:
// ############################# Public Method ###############################
float SpectralDisplay::acqEndFrequencyFor(int acqNumber)
{
  int   acq_points, number_acqs;
  float start, spacing;
//
// Check for end condition:
//
  number_acqs   = numberAcquisitions();
  if((acqNumber+1)==number_acqs)
  {
    return (_centerFrequency + _spanRange/2.);
  }
//
// End frequency is based on start frequency, # of acq points,
// and frequency spacing.
//
  start         = acqStartFrequencyFor(acqNumber);
  acq_points    = acqPointsFor(acqNumber, start);
  spacing       = frequencySpacing();
  return (start + acq_points*spacing);
}

// ############################# Public Method ###############################
// acqPointsFor -- Calculates and returns the number of spectral points for the
//                         given acquisition
//
// Input:       acqNumber:      The acquisition # starting from 0
//              startFreq:      Pre-calculated start frequency
//              
//          
// Output:                      # of acquisition points
//
// Notes:
// ############################# Public Method ###############################
int SpectralDisplay::acqPointsFor(int acqNumber, float startFreq)
{
  int   acq_points, number_acqs;
  float start, center, end_freq, spacing, float_points;
//
// Calculate spectral point spacing:
//
  spacing       = frequencySpacing();
//
// Check for end condition:
//
  number_acqs   = numberAcquisitions();
  if((acqNumber+1)==number_acqs)
  {
    end_freq    = acqEndFrequencyFor(acqNumber);
    start       = acqStartFrequencyFor(acqNumber);
    acq_points  = ROUND((end_freq-start)/spacing);
    return acq_points;
  }

  if(acqNumber < 1)
  {
    float_points        = _acquisitionSpan/spacing;
    acq_points          = ROUND(float_points - 0.5);    // round down
    return acq_points;
  }
//
// # of points is based on center frequency, start frequency, and acq span.
//
  center        = acqCenterFrequencyFor(acqNumber);
  if(startFreq < 0.)
    start       = acqStartFrequencyFor(acqNumber);
  else
    start       = startFreq;
  float_points  = (center - start + _acquisitionSpan/2.)/spacing;
  acq_points    = ROUND(float_points - 0.5);            // round down
  return acq_points;
}

// ############################# Public Method ###############################
// acqOverlapPointsFor -- Calculates and returns the number of spectral points
//                        in the overlap region (with last acq) for a given acquisition
//
// Input:       acqNumber:      The acquisition # starting from 0
//              lastEnd:        Previous acquisition end frequency
//              
//          
// Output:                      # overlap points
//
// Notes:
// ############################# Public Method ###############################
int SpectralDisplay::acqOverlapPointsFor(int acqNumber, float lastEnd)
{
  int   overlap_points;
  float last_center, last_end, spacing, float_points, overlap;

  if(acqNumber < 1)
    return 0;                                           // Not defined for 1st acq
//
// # of points is based on center frequency, start frequency, and acq span.
//
  spacing               = frequencySpacing();
  overlap               = overlapRange();
  if(lastEnd < 0.)
    last_end            = acqEndFrequencyFor(acqNumber-1);
  else
    last_end            = lastEnd;
  last_center           = acqCenterFrequencyFor(acqNumber-1);
  float_points          = (overlap - (last_center + _acquisitionSpan/2. - last_end))/spacing;
  overlap_points        = ROUND(float_points - 0.5);            // round down
  overlap_points++;
  return overlap_points;
}
