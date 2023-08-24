#ifndef MOTIF_VECTOR_H
#define MOTIF_VECTOR_H

#include "MotifHit.h"

// 动态数组结构体
typedef struct {
    MotifHit* hits;
    int size;
    int capacity;
} MotifHitVector;

size_t motifHitVectorSize(const MotifHitVector* vec);
void printMotifHitVector(const MotifHitVector* vec);
void initMotifHitVector(MotifHitVector* vec);
void pushMotifHitVector(MotifHitVector* vec, MotifHit hit);

int compareMotifHitsByPVal(const void* a, const void* b);
void sortMotifHitVectorByPVal(MotifHitVector* vec);

// Retains only the top k elements in the vector.
void retainTopKMotifHits(MotifHitVector* vec, size_t k);

void removeHitAtIndex(MotifHitVector* vec, size_t indx);

void clearMotifHitVector(MotifHitVector* vec);
void freeMotifVector(MotifHitVector* vec);

#endif /* MOTIF_VECTOR_H */
