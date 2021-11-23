/* ----------------------- Inclusion of Own Headers ------------------------ */
#include "complex.h"
#include "fft.h"
/* ---------------------- Constants to Exit Functions ----------------------- */
#define SUCCEED 1 /* Succeeded ending in function execution */
#define FAIL 0 /* Failed ending in function execution */
/* ---------------------- Analysis in Frequency Domain ---------------------- */
#define FORWARD 1 /* Forward direction for FFT computation */
#define REVERSE 0 /* Reverse direction for FFT computation */
/********** Symbolic Constants ***********/
/* -------- Default Values for Iteration in Time & Frequency Domain --------- */
#define SAMPLES 256


static COMPLEX_T Signal[SAMPLES];
static COMPLEX_T Signal_FFT[SAMPLES];

const int analoginput = A0;
int entrada = 0;
bool rd = false; // lectura realizada
int counter = 0;

int dir = FORWARD;
int power;


void transmitir_datos(COMPLEX_T* Signal){
  
  for ( int i = 0; i < SAMPLES / 2; i++ ){
    Serial.println(cplx_Magnitude( Signal[i])/100000);
  }
}

void setup() {
  Serial.begin(9600); // serial init
  // setup interrupts
  cli();//stop interrupts
  TCCR4A = 0;
  TCCR4B = 0;
  TCNT4 = 0;
  OCR4A = 50 / 1; 
  TCCR4B |= (1 << WGM12);
  TCCR4B |= (0 << CS42) | (1 << CS41) | (0 << CS40);
  TIMSK4 |= (1 << OCIE4A);  
  sei();//allow interrupts

  power = Find_Power( SAMPLES );
}

void loop() {
  if(rd){
    Signal[counter].real = entrada/1;
    counter++;
    rd = false;
  }
  if(counter > 255){
    cli();//stop interrupts
    for(int j = 0; j<SAMPLES;j++){
      Signal_FFT[j] = Signal[j];
    }
  Compute_FFT( dir, power, Signal_FFT );
  transmitir_datos(Signal_FFT);
  counter = 0;
  sei();//allow interrupts
  }
}

// Timer Interrupt
ISR(TIMER4_COMPA_vect) {
  entrada = analogRead(analoginput);
//  Serial.println(entrada);
  rd = true;
}
