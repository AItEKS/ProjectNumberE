#ifndef BigInt_h
#define BigInt_h
#include "FFT.h"

typedef struct BigIntTag {
	Real * Coef;
	long Size, SizeMax;
} BigInt;

#define BASE 10000.
#define invBASE 0.0001
#define NBDEC_BASE 4

void InitializeBigInt(BigInt * A, long MaxSize);
void FreeBigInt(BigInt * A);
void PrintBigInt(BigInt * A);
void UpdateBigInt(BigInt * A);
void AddBigInt(BigInt * A, BigInt * B, BigInt * C);
void MulBigInt(BigInt * A, BigInt * B, BigInt * C);
void Inverse(BigInt * A, BigInt * B, BigInt * tmpBigInt);

#endif