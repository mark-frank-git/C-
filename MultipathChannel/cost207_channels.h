#ifndef _COST207_CHANNELS_H
#define _COST207_CHANNELS_H    1
/************************************************************************
 *                                                                      *
 * This header file defines the channels specified in COST 207.         *
 *                                                                      *
 *                                                                      *
 * File:cost207_channels.h                                              *
 *                                                                      *
 *                                                                      *
 * See: GSM 05.05 V8.2.0, Annex C.                                      *
 *                                                                      *
 ************************************************************************/
#include "environment_types.h"

//
// Number of paths is equal to the number of taps to use in the delay line:
//
int     number_paths[NUMBER_ENV_TYPES]  =
{
  4,                    // RA
  12,                   // HT
  12,                   // TU
  6,                    // EQ
  2,                    // TI
  1,                    // Static
  2,                    // Test 1
  2,                    // Test 2
  12,                   // Test 3
  4                     // Test 4
};
#define TOTAL_ENV_TAPS  57              // 4 + 12 + 12 + 6 + 2 + 1 + 2 + 2 + 12 + 4
//
// The tap times are given in microseconds
//
float tap_times[TOTAL_ENV_TAPS] =
{
  0., 0.2, 0.4, 0.6,                                                    // RA
  0., 0.2, 0.4, 0.6, 0.8, 2.0, 2.4, 15.0, 15.2, 15.8, 17.2, 20.,        // HT
  0., 0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2.4, 3.0, 3.2, 5.0,            // TU
  0., 3.2, 6.4, 9.6, 12.8, 16.0,                                        // EQ
  0., 0.4,                                                              // TI
  0.,                                                                   // Static
  0., 2.,                                                               // TEST1
  0., 2.,                                                               // TEST2
  0., 0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2.4, 3.0, 3.2, 5.0,            // TEST3
  0., 1., 2., 10.                                                       // TEST4
};
//
// The tap powers are given in relative dB
//
float tap_powers[TOTAL_ENV_TAPS] =
{
  0., -2.0, -10.0, -20.0,                                                       // RA
  -10., -8.0, -6.0, -4.0, 0.0, 0.0, -4.0, -8.0, -9.0, -10.0, -12.0, -14.0,      // HT
  -4.0, -3.0, 0.0, -2.0, -3.0, -5.0, -7.0, -5.0, -6.0, -9.0, -11.0, -10.0,      // TU
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                                                 // EQ
  0.0, 0.0,                                                                     // TI
  0.0,                                                                          // Static
  0.0, -80.,                                                                    // Test1
  -3., -80.,                                                                    // Test2
  -4.0, -3.0, 0.0, -2.0, -3.0, -5.0, -7.0, -5.0, -6.0, -9.0, -11.0, -10.0,      // Test3
  -80., -80., -80., 0.                                                          // Test4
};

int     filter_types[TOTAL_ENV_TAPS] =
{
  RICIAN_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE,                                          // RA
  CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE,
  CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE,       // HT
  CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE,
  CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE,       // TU
  CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE,       // EQ
  CLASSICAL_TYPE, CLASSICAL_TYPE,                                                                       // TI
  NO_FADING_TYPE,                                                                                       // Static
  CLASSICAL_TYPE, CLASSICAL_TYPE,                                                                       // TEST1
  CLASSICAL_TYPE, CLASSICAL_TYPE,                                                                       // TEST2
  CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE,
  CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE, CLASSICAL_TYPE,       // TEST3
  RICIAN_TYPE, RICIAN_TYPE, RICIAN_TYPE, RICIAN_TYPE                                                    // TEST4
};


#endif
