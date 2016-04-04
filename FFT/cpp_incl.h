#ifndef _CPP_INCL_H
#define _CPP_INCL_H 1
/************************************************************************************
 *                                                                                  *
 * File: cpp_incl.h                                                                 *
 *                                                                                  *
 *  This include file is used for including all the C++ prototype files. It is used *
 *  so that DOS versions can replace this file with a set of DOS file names         *
 *  (8 chars)                                                                       *
 *                                                                                  *
 *  Revision history:                                                               *
 *  1. 12/20/94 - Started.                                                          *
 *  2. 12/30/94 - Add fade simulator objects                                        *
 *  3. 01/06/95 - Add C include files -> delete c_incl.h.                           *
 *  4. 04/04/95 - Add Modem.h, SVModem.h, ISUModem.h                                *
 *  5. 07/04/95 - Use local copies from GNU C++ library.                            *
 *  6. 07/05/95 - Add CallState.h, ISUCallState.h SVCallState.h.                    *
 *  7. 08/30/95 - Add BaseISU.h                                                     *
 *  8. 09/26/95 - Add EAT_CANDIDATE_LIST.                                           *
 *  9. 10/10/95 - Add LogHist.h                                                     *
 ************************************************************************************/
#include <C_Libraries/constants.h>                      // C type file
#include <Specfuns/specfuns.h>

#include <Filters/DataWindow.h>
#include "FFT.h"
#include "PowerSpectrum.h"


// GNU C++ objects
#include <GNU/Complex.h>


#endif
