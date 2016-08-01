#ifndef _DATA_WINDOW_H
#define _DATA_WINDOW_H  1
/************************************************************************
 *                                                                      *
 * This class implements a data windowing function.  It is based on the *
 * paper, "On the use of windows for harmonic analysis with the dis-    *
 * crete Fourier transform," by F.J. Harris, proc. IEEE vol. 66, no.1.  *
 *                                                                      *
 * File:DataWindow.h                                                    *
 *                                                                      *
 ************************************************************************/
//
// The following are the types of windows that are available
//
#define RECTANGULAR             0               // w(n) = 1
#define TRIANGULAR              1               // Also known as Bartlett
#define HANNING                 2               // Cos()^alpha
#define HAMMING                 3               // Modified Hanning
#define BLACKMAN_HARRIS3        4               // 3 term Blackman-Harris
#define BLACKMAN_HARRIS4        5               // 4 term Blackman-Harris
#define GAUSSIAN_WINDOW         6               // Minimum time-bandwidth
#define KAISER                  7               // Kaiser Bessel (max main lobe energy)
#define TUKEY_WINDOW            8               // Tukey Cosine taper window (for transitioning)
#define RAISED_COSINE           9               // EDGE raised cosine window
#define MAX_WINDOW_TYPE         RAISED_COSINE

#define NUMBER_WINDOW_TYPES (RAISED_COSINE+1)

//
// The following are default values for window parameters
//
#define KAISER_ALPHA            2.5
#define GAUSSIAN_ALPHA          3.
#define HANNING_ALPHA           2
#define TUKEY_ALPHA             0.25
#define RAISED_SAMPLES          4

class   Complex;

#ifndef OBJECTIVE_C_CODE
class DataWindow
{
private:
  int           _windowType;                            // Windowing type, see above
  int           _hanningAlpha;

  float         _kaiserAlpha;
  float         _gaussianAlpha;
  float         _tukeyAlpha;
  float         _raisedSamples;                         // Samples/symbol for raised cosine

  double        *_window;                               // Window function
  double        _coherentGain;                          // Normalization factor for CW signals
  double        _noncoherentGain;                       // Normalization factor for noise like signals
//
// Private functions:
//
  Complex       dirichletKernel(double omega, int length); // Local for dftMagnitudeAt()
  double        dirichletMagnitude(double omega, int length); // Local for dftMagnitudeAt()
// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  DataWindow(int type=RECTANGULAR);
  ~DataWindow();
        
/**********************
 * Set parameters:    *
 **********************/
  void          setWindowType(int type)         {_windowType    = type;         return;}
  void          setHanningAlpha(int alpha)      {_hanningAlpha  = alpha;        return;}
  void          setGaussianAlpha(float alpha)   {_gaussianAlpha = alpha;        return;}
  void          setKaiserBeta(float alpha)      {_kaiserAlpha   = alpha;        return;}
  void          setTukeyAlpha(float alpha)      {_tukeyAlpha    = alpha;        return;}
  void          setRaisedCosineSamplesPerSymbol(float samples)
                                                {_raisedSamples = samples;      return;}

/**********************
 * Get parameters:    *
 **********************/
  int           windowType()                    {return _windowType;}
  double        coherentGain()                  {return _coherentGain;}
  double        noncoherentGain()               {return _noncoherentGain;}
  float         bandwidth6dB();
  float         bandwidth3dB();

/**********************
 * Getting Outputs:   *
 **********************/
  const         double *windowFunction(int length);         // Returns the window coefficients
  double        dftMagnitudeSquaredAt(double omega, int length);
private:
  void findTriangularWindow(
                            const int length,
                            double    *window) const;
  void findHanningWindow (
                            const int   length,
                                                                                                                const float hannningAlpha,
                            double       *window) const;
  void findHammingWindow (
                            const int   length,
                            double       *window) const;

  void findKaiserWindow (
                            const int   length,
                                                                                                                const float kaiserAlpha,
                            double       *window) const;

  void findGaussianWindow (
                            const int   length,
                                                                                                                const float gaussianAlpha,
                            double       *window) const;

  void findTukeyWindow (
                            const int   length,
                                                                                                                const float tukeyAlpha,
                            double       *window) const;

  void findRaisedCosineWindow (
                              const int         length,
                              const float samplesPerSymbol,
                              double     *window) const;
};

#endif
#endif
