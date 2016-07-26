#ifndef _OVSF_GENERATOR_H
#define _OVSF_GENERATOR_H       1
/********************************************************************************
 * This subclass of Object generates orthogonal variable spreading factor (VSF) *
 * codes for WCDMA and stores them in arrays for output.                        *
 *                                                                              *
 * File: OVSFGenerator.h                                                        *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 06/27/01 - Started                                                      *
 *                                                                              *
 * Notes:                                                                       *
 * 1. The VSF code functions have values of 0 or 1.                             *
 * 2. Reference is "Spreading codes for direct sequence CDMA and wideband       *
 *    CDMA cellular networks," Dinan and Jabbari, IEEE Comm. Mag. Sept. 1998.   *
 *                                                                              *
 ********************************************************************************/
#ifndef LOGICAL
#define LOGICAL char
#endif

#define DEFAULT_LENGTH  64

class OVSFGenerator
{
protected:
  int           _codeLength;                    // Length of code to be generated
  int           _codeRow;                       // Row of the code to be generated
  char          *_vsfCode;                      // The generated code

//
// Private functions:
//
  LOGICAL       isPowerOfTwo(int n);
  void          generateVSFCode(char *vsf, int length, int row);

//
// Public functions:
//
public:
    
/********************************
 * Constructors, destructors    *
 ********************************/
  OVSFGenerator(int length = DEFAULT_LENGTH, int row = 1);
  ~OVSFGenerator();

/********************************
 * Getting a code:              *
 ********************************/
  const char    *vsfCode()      {return _vsfCode;}

};
#endif

