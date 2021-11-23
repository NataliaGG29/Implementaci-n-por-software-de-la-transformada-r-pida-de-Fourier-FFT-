/* Ruler 1         2         3         4         5         6         7        */

/**********************************  fft.h  ***********************************/
/*                                                                            */
/*   Purpose: Header to compute FFT                                           */
/*                                                                            */
/******************************************************************************/

#ifndef FFT_H
#define FFT_H



/********************************** Headers ***********************************/

/* ------------------------ Inclusion of Own Headers ------------------------ */

#include "complex.h"



/**************************** Symbolic  Constants *****************************/

/* ---------------------- Analysis in Frequency Domain ---------------------- */

#define FORWARD              1    /* Forward direction for FFT computation    */
#define REVERSE              0    /* Reverse direction for FFT computation    */



/************************** Prototypes of Functions ***************************/

/* ---------------------------- Public Functions ---------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

void Compute_FFT( int, int, COMPLEX_T * );
int Find_Power( long );
void Compute_FFT_fsm( int, int, COMPLEX_T * );
#ifdef __cplusplus
}
#endif
#endif /* FFT_H */
