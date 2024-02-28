#include "FFT.c"
#include "BigInt.c"
#include <stdio.h>
#include <math.h>
#include <time.h>

BigInt tmpBigInt;

void BinarySplittingE(long N0, long N1, BigInt* P, BigInt* Q) {
    BigInt PP, QQ;
    long NMid;

    // Base case: when N1 - N0 equals 1
    if (N1 - N0 == 1) {
        P->Size = Q->Size = 1;
        P->Coef[0] = 1.;
        Q->Coef[0] = (double)N1;
        UpdateBigInt(P);
        UpdateBigInt(Q);
        return;
    }

    NMid = (N0 + N1) / 2;

    // Recursive call for the first half of the range
    BinarySplittingE(N0, NMid, P, Q);

    // Set the Coef and SizeMax for the second half of the range
    PP.Coef = P->Coef + P->Size;
    PP.SizeMax = P->SizeMax - P->Size;
    QQ.Coef = Q->Coef + Q->Size;
    QQ.SizeMax = Q->SizeMax - Q->Size;

    // Recursive call for the second half of the range
    BinarySplittingE(NMid, N1, &PP, &QQ);

    // Multiply P and QQ, and add the result to P
    MulBigInt(P, &QQ, &tmpBigInt);
    AddBigInt(&tmpBigInt, &PP, P);

    // Multiply Q and QQ
    MulBigInt(Q, &QQ, Q);
}

BigInt ECompute(long NbDec) {
    BigInt P, Q, tmp;
    long i, MaxSize, MaxFFTSize, SeriesSize;
    double logFactorial, logMax;

    // Calculate the maximum size for P and Q
    MaxSize = NbDec / NBDEC_BASE + 10 + (long)(2. * log((double)NbDec) / log(2.));

    // Initialize P and Q
    InitializeBigInt(&P, MaxSize);
    InitializeBigInt(&Q, MaxSize);

    // Calculate the maximum FFT size
    MaxFFTSize = 2;
    while (MaxFFTSize < MaxSize) MaxFFTSize *= 2;
    MaxFFTSize *= 2;

    // Initialize FFT and temporary BigInts
    InitializeFFT(MaxFFTSize);
    InitializeBigInt(&tmpBigInt, MaxFFTSize);
    InitializeBigInt(&tmp, MaxFFTSize);

    SeriesSize = 1;
    logFactorial = 0.;
    logMax = (double)NbDec * log(10.);

    // Calculate the series size based on the logFactorial and logMax
    while (logFactorial < logMax) {
        logFactorial += log((double)SeriesSize);
        SeriesSize++;
    }

    // Perform Binary Splitting algorithm
    BinarySplittingE(0, SeriesSize, &P, &Q);

    // Add P and Q
    AddBigInt(&P, &Q, &P);

    // Calculate the inverse of Q
    Inverse(&Q, &tmpBigInt, &tmp);

    // Multiply P and the inverse of Q
    MulBigInt(&P, &tmpBigInt, &tmpBigInt);

    // Set the size and coefficients of P
    P.Size = 1 + NbDec / NBDEC_BASE;
    for (i = 1; i <= P.Size; i++) P.Coef[P.Size - i] = tmpBigInt.Coef[tmpBigInt.Size - i];

    // Free memory for Q and temporary BigInts
    FreeBigInt(&Q);
    FreeBigInt(&tmpBigInt);

    return P;
}

int main() {
    double StartTime;
    BigInt E;
    long NbDec = 1000000;

    // Start the timer
    StartTime = (double)clock();

    // Compute E
    E = ECompute(NbDec);

    // Print the time taken and the value of E
    printf("Time: %.2lf seconds\n", ((double)clock() - StartTime) / CLOCKS_PER_SEC);
    PrintBigInt(&E);
}