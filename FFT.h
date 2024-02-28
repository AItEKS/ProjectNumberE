#ifndef FFT_h
#define FFT_h

typedef double Real;
extern double FFTSquareWorstError;
extern long AllocatedMemory;
void InitializeFFT(long MaxFFTLength);
void MulWithFFT(Real * ACoef, long ASize, Real * BCoef, long BSize, Real * CCoef);

#endif