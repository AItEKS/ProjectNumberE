#include "FFT.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926535897932384626

typedef struct ComplexTag {
    double real;
    double imaginary;
} Complex;

long FFTLengthMax;
Complex* OmegaFFT;
Complex* ArrayFFT0, * ArrayFFT1;
Complex* ComplexCoef;
double FFTSquareWorstError;
long AllocatedMemory;

void InitializeFFT(long MaxLength) {
    long i;
    double Step;
    FFTLengthMax = MaxLength;
    OmegaFFT = (Complex*)malloc(MaxLength / 2 * sizeof(Complex));
    ArrayFFT0 = (Complex*)malloc(MaxLength * sizeof(Complex));
    ArrayFFT1 = (Complex*)malloc(MaxLength * sizeof(Complex));
    ComplexCoef = (Complex*)malloc(MaxLength * sizeof(Complex));
    Step = 2.0 * PI / (double)MaxLength;
    for (i = 0; 2 * i < MaxLength; i++) {
        OmegaFFT[i].real = cos(Step * (double)i);
        OmegaFFT[i].imaginary = sin(Step * (double)i);
    }
    FFTSquareWorstError = 0.0;
    AllocatedMemory = 7 * MaxLength * sizeof(Complex) / 2;
}

void RecursiveFFT(Complex* Coef, Complex* FFT, long Length, long Step, long Sign) {
    long i, OmegaStep;
    Complex* FFT0, * FFT1, * Omega;
    double tmpReal, tmpImaginary;
    if (Length == 2) {
        FFT[0].real = Coef[0].real + Coef[Step].real;
        FFT[0].imaginary = Coef[0].imaginary + Coef[Step].imaginary;
        FFT[1].real = Coef[0].real - Coef[Step].real;
        FFT[1].imaginary = Coef[0].imaginary - Coef[Step].imaginary;
        return;
    }
    FFT0 = FFT;
    RecursiveFFT(Coef, FFT0, Length / 2, Step * 2, Sign);
    FFT1 = FFT + Length / 2;
    RecursiveFFT(Coef + Step, FFT1, Length / 2, Step * 2, Sign);
    Omega = OmegaFFT;
    OmegaStep = FFTLengthMax / Length;
    for (i = 0; 2 * i < Length; i++, Omega += OmegaStep) {
        tmpReal = Omega[0].real * FFT1[i].real - Sign * Omega[0].imaginary * FFT1[i].imaginary;
        tmpImaginary = Omega[0].real * FFT1[i].imaginary + Sign * Omega[0].imaginary * FFT1[i].real;
        FFT1[i].real = FFT0[i].real - tmpReal;
        FFT1[i].imaginary = FFT0[i].imaginary - tmpImaginary;
        FFT0[i].real = FFT0[i].real + tmpReal;
        FFT0[i].imaginary = FFT0[i].imaginary + tmpImaginary;
    }
}

void FFT(double* Coef, long Length, Complex* FFT, long NFFT) {
    long i;
    for (i = 0; i < Length; i++) {
        ComplexCoef[i].real = Coef[i];
        ComplexCoef[i].imaginary = 0.0;
    }
    for (; i < NFFT; i++) {
        ComplexCoef[i].real = ComplexCoef[i].imaginary = 0.0;
    }
    RecursiveFFT(ComplexCoef, FFT, NFFT, 1, 1);
}

void InverseFFT(Complex* FFT, long NFFT, double* Coef, long Length) {
    long i;
    double invNFFT = 1.0 / (double)NFFT;
    double tmp;
    RecursiveFFT(FFT, ComplexCoef, NFFT, 1, -1);
    for (i = 0; i < Length; i++) {
        tmp = invNFFT * ComplexCoef[i].real;
        Coef[i] = floor(0.5 + tmp);
        if ((tmp - Coef[i]) * (tmp - Coef[i]) > FFTSquareWorstError) {
            FFTSquareWorstError = (tmp - Coef[i]) * (tmp - Coef[i]);
        }
    }
}

void Convolution(Complex* A, Complex* B, long NFFT, Complex* C) {
    long i;
    double tmpReal, tmpImaginary;
    for (i = 0; i < NFFT; i++) {
        tmpReal = A[i].real * B[i].real - A[i].imaginary * B[i].imaginary;
        tmpImaginary = A[i].real * B[i].imaginary + A[i].imaginary * B[i].real;
        C[i].real = tmpReal;
        C[i].imaginary = tmpImaginary;
    }
}

void MulWithFFT(double* ACoef, long ASize, double* BCoef, long BSize, double* CCoef) {
    long NFFT = 2;
    while (NFFT < ASize + BSize) {
        NFFT *= 2;
    }
    FFT(ACoef, ASize, ArrayFFT0, NFFT);
    FFT(BCoef, BSize, ArrayFFT1, NFFT);
    Convolution(ArrayFFT0, ArrayFFT1, NFFT, ArrayFFT0);
    InverseFFT(ArrayFFT0, NFFT, CCoef, ASize + BSize - 1);
}