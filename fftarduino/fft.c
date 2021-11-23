/* Ruler 1         2         3         4         5         6         7        */

/********************************  complex.c  *********************************/
/*                                                                            */
/*   Purpose: Module to perform some basic operations with complex numbers    */
/*                                                                            */
/*   Origin:  Module tested and documented by jcgiraldo, since June 4th, 2013 */
/*                                                                            */
/*   Notes:   None                                                            */
/*                                                                            */
/*   DATE         RESPONSIBLE  COMMENT                                        */
/*   -----------------------------------------------------------------------  */
/*   Mar 01/2019  J.C.Giraldo  cplx_Divide was incorporated                   */
/*   Nov 07/2017  J.C.Giraldo  complex was changed to COMPLEX_T               */
/*   Nov 07/2017  J.C.Giraldo  cplx_Phase was incorporated                    */
/*   Nov 07/2017  J.C.Giraldo  cplx_Ratio was incorporated                    */
/*   Jun 05/2013  J.C.Giraldo  Initial implementation                         */
/*                                                                            */
/******************************************************************************/

/*********************** Directives of C Pre-processor ************************/

/************************** Conditional Compilation ***************************/

#define DEBUGGING                 /* Use directive when debugging some lines  */



/**************************** Symbolic  Constants *****************************/

/* ---------------------- Constants to Exit Functions ----------------------- */

#define SUCCEED              1    /* Succeeded ending in function execution   */
#define FAIL                 0    /* Failed    ending in function execution   */



/********************************** Headers ***********************************/

/* ------------------------ Inclusion of Std Headers ------------------------ */

#include <math.h>   /* Due to sqrt */

/* ----------------------- Inclusion of Other Headers ----------------------- */

#include "complex.h"
#include "fft.h"



/*****************************  Public Functions  *****************************/

/* ---------------------- Analysis in Frequency Domain ---------------------- */

/*FN****************************************************************************
*
*   int
*   Compute_FFT( int dir, int pow, COMPLEX_T *x );
*
*   Return:          Transformation in values passed by reference
*
*   Purpose:         Compute an in-place complex-to-complex FFT
*
*   Note:            This function computes an in-place complex-to-complex FFT
*                    x.real and y.imag are the real and imaginary arrays of
*                    2^pow points.
*                    dir = 1 or FORWARD gives forward transform
*                    dir = 0 or REVERSE gives reverse transform
*
*                    There is a modification by Peter Cusak to utilize the
*                    Microsoft complex type.
*
*   Plan:
*           Part 1: Calculate the number of points
*           Part 2: Do the bit reversal
*           Part 3: Compute the FFT
*           Part 4: Scale for forward transformation
*
*   Register of Revisions (Debugging Process):
*
*   DATE         RESPONSIBLE  COMMENT
*   -------------------------------------------------------------------------
*   Jun 04/2013  J.C.Giraldo  Incorporation of functions with complex numbers
*   May 03/2013  J.C.Giraldo  Readable identifiers
*   Jun --/1993  P. Bourkes   Initial implementation
*
*******************************************************************************/
void
Compute_FFT_fsm( int dir, int pow, COMPLEX_T *x )
{
	typedef enum { ZERO,ONE,TWO,TREE,FOUR,FIVE, SIX } STATE_T;
	STATE_T state = ZERO;
	int FIN =1 ;
	//INSTRCCION_A
	long    points, i, j, k, l, i1, i2, l1, l2;
	COMPLEX_T c, temp, u;
	i=0;
	points = 1; 	

	while( FIN == 1 ) {
		
    	switch( state ) {

    		case ZERO:
          
			if( i < pow) {
				//INSTRCCION_B
				points <<= 1;
				i++;
                		state = ZERO;
            		} 
			if (i >= pow){
				//INSTRCCION_C
                		i2 = points >> 1;
				j = 0;
				i= 0;
				state = ONE;
           		}
			break;
					
		case ONE:  				
			if( i< points-1 ) {
				//INSTRCCION_D
          			if(i<j) cplx_Swap( &x[i], &x[j] );
          			k=i2;
               			state = TWO;
           		}
           		if( i >= points-1 ) {
           			//INSTRCCION_H
          			c.real = -1.0;
          			c.imag = 0.0;
          			l2 = 1;
				l= 0;
               			state = TREE;
           		}
        		break;
			
		case TWO:  				
			if( k<= j ) {
				//INSTRCCION_F
          			j -= k;
				k >>= 1;
               			state = TWO;
           		}
           		else{
           			//INSTRCCION_G
				j +=k;
       				i++;
				state = ONE;
			}
        	break;
        	
		case TREE:  				
			if( l>=pow && dir != FORWARD ) {
               			FIN = 0;
           		}
           		if( l< pow ) {
           			//INSTRCCION_I
          			l1 = l2;
    				l2 <<= 1;
				u.real = 1.0;
   		 		u.imag = 0.0;
          			j=0;
               			state = FOUR;
           		}
           		if( l >= pow && dir == FORWARD ) {
          			i=0;
               			state = SIX;
           		}
           	break;
			   	
           	case FOUR:
           		if( j<l1 ) {
          			i=j;
               			state = FIVE;
			}
			if( j>=l1 ) {
				//INSTRCCION_L
        			c.imag = sqrt( ( 1.0 - c.real ) / 2.0 );
    				if( dir == FORWARD ) c.imag = -c.imag;
    				c.real = sqrt( ( 1.0 + c.real ) / 2.0 );
          			l++;
        	       		state = TREE;
           		}
        	break;
			case FIVE:
           		if( i<points ) {
    				//INSTRUCCION_J;
          			i1    = i + l1;
        	    		temp  = cplx_Multiply( u, x[i1] );
        	    		x[i1] = cplx_Subtract( x[i], temp );
        	    		x[i]  = cplx_Add( x[i], temp );
				i=i+l2;
        	       		state = FIVE;
			}
			if( i>=points ) {
    				//INSTRUCCION_K;
				u = cplx_Multiply( u, c );
          			j++;
               			state = FOUR;
           		}
        	break;
		case SIX:
           		if( i<points ) {
    				//INSTRUCCION_M;
    				x[i] = cplx_Scale( x[i], points );
				i++;
        	       		state = SIX;
			}
			if( i>=points ) {
          			FIN = 0;
           		}
        	break;
    	} /* switch */
	} /* while */
}


void
Compute_FFT( int dir, int pow, COMPLEX_T *x )
{
long    points, i, j, k, l, i1, i2, l1, l2;
COMPLEX_T c, temp, u;

/* Part 1: Calculate the number of points */

for( points = 1, i = 0; i < pow; i++ ) points <<= 1;

/* Part 2: Do the bit reversal */

i2 = points >> 1;
j = 0;
for( i = 0; i < points-1; i++ ) {
    if( i < j ) cplx_Swap( &x[i], &x[j] );
    k = i2;
    while( k <= j ) j -= k, k >>= 1;
    j += k;
}

/* Part 3: Compute the FFT */

c.real = -1.0;
c.imag =  0.0;

l2 =  1;
for( l = 0; l < pow; l++ ) {
    l1 = l2;
    l2 <<= 1;
    u.real = 1.0;
    u.imag = 0.0;
    for( j = 0; j < l1; j++) {
        for( i = j; i < points; i += l2 ) {
            i1    = i + l1;
            temp  = cplx_Multiply( u, x[i1] );
            x[i1] = cplx_Subtract( x[i], temp );
            x[i]  = cplx_Add( x[i], temp );
        }
        u = cplx_Multiply( u, c );
    }
    c.imag = sqrt( ( 1.0 - c.real ) / 2.0 );
    if( dir == FORWARD ) c.imag = -c.imag;
    c.real = sqrt( ( 1.0 + c.real ) / 2.0 );
}

/* Part 4: Scale for forward transformation */

if( dir == FORWARD )
    for( i = 0; i < points; i++ )
        x[i] = cplx_Scale( x[i], points );

} /* Compute_FFT */



/*FN****************************************************************************
*
*   int
*   Find_Power( int number );
*
*   Purpose:         Return the (integer) logarithm in base 2 of a number
*
*   Note:            Number must be an integer equal to zero or greater
*                    Another option is shifting and truncation
*
*   Register of Revisions (Debugging Process):
*
*   DATE       RESPONSIBLE  COMMENT
*   -----------------------------------------------------------------------
*   May 03/13  J.C.Giraldo  Readable identifiers
*   --- --/--  P. Bourkes   Initial implementation
*
*******************************************************************************/

int
Find_Power( long number )
{
int power = 0;

while( number != 1 ) {
    number = (long)( number/2 );
    power++;
}

return( power );

} /* Find_Power */
