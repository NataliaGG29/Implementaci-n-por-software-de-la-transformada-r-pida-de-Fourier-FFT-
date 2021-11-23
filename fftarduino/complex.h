/* Ruler 1         2         3         4         5         6         7        */

/********************************  complex.h  *********************************/
/*                                                                            */
/*   Purpose: Header to operate complex numbers                               */
/*                                                                            */
/******************************************************************************/

#ifndef COMPLEX_H
#define COMPLEX_H

/************************** Definition of Data Types **************************/

/* ---------------------------- Complex Numbers ----------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	double real; /* Real part of type complex number      */
	double imag; /* Imaginary part of type complex number */
} COMPLEX_T;


/************************** Prototypes of Functions ***************************/

/* ---------------------------- Public Functions ---------------------------- */

COMPLEX_T cplx_Set     ( double, double );
double    cplx_Get_Real( COMPLEX_T );
double    cplx_Get_Imag( COMPLEX_T );
int       cplx_Swap    ( COMPLEX_T *, COMPLEX_T * );

COMPLEX_T cplx_Null();
COMPLEX_T cplx_Multiply( COMPLEX_T, COMPLEX_T );
COMPLEX_T cplx_Divide  ( COMPLEX_T, COMPLEX_T );
COMPLEX_T cplx_Subtract( COMPLEX_T, COMPLEX_T );
COMPLEX_T cplx_Add     ( COMPLEX_T, COMPLEX_T );
COMPLEX_T cplx_Scale   ( COMPLEX_T, double );

double cplx_Magnitude  ( COMPLEX_T );
double cplx_Phase      ( COMPLEX_T );
double cplx_Ratio      ( COMPLEX_T, COMPLEX_T );

#ifdef __cplusplus
}
#endif

#endif /* COMPLEX_H */
