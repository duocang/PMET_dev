#ifndef MOTIF_VECTOR_H
#define MOTIF_VECTOR_H

#include "pmet-index-MotifHit.h"

// 动态数组结构体
typedef struct {
    MotifHit* hits;
    int size;
    int capacity;
} MotifHitVector;


/**
 * Returns the current size (number of elements) of the MotifHitVector.
 * @param vec Pointer to the MotifHitVector.
 * @return size_t The size of the vector.
 */
size_t motifHitVectorSize(const MotifHitVector* vec);

/**
 * Prints the contents of the MotifHitVector.
 * @param vec Pointer to the MotifHitVector to be printed.
 */
void printMotifHitVector(const MotifHitVector* vec);

/**
 * Initializes the given MotifHitVector.
 * @param vec Pointer to the MotifHitVector to be initialized.
 */
void initMotifHitVector(MotifHitVector* vec);

/**
 * Adds a new MotifHit to the end of the MotifHitVector.
 * @param vec Pointer to the MotifHitVector.
 * @param hit The MotifHit to be added.
 */
void pushMotifHitVector(MotifHitVector* vec, const MotifHit* hit);

/**
 * Compares two MotifHits based on their p-values.
 * @param a Pointer to the first MotifHit.
 * @param b Pointer to the second MotifHit.
 * @return int Negative if 'a' is less than 'b', 0 if they're equal, positive if 'a' is greater than 'b'.
 */
int compareMotifHitsByPVal(const void* a, const void* b);

/**
 * Sorts the MotifHitVector based on the p-values of the contained MotifHits.
 * @param vec Pointer to the MotifHitVector to be sorted.
 */
void sortMotifHitVectorByPVal(MotifHitVector* vec);

/**
 * Retains only the top 'k' MotifHits in the MotifHitVector based on their order.
 * @param vec Pointer to the MotifHitVector.
 * @param k The number of MotifHits to retain.
 */
void retainTopKMotifHits(MotifHitVector* vec, size_t k);

/**
 * Removes the MotifHit at the specified index from the MotifHitVector.
 * @param vec Pointer to the MotifHitVector.
 * @param indx Index of the MotifHit to be removed.
 */
void removeHitAtIndex(MotifHitVector* vec, size_t indx);

/**
 * Clears all the MotifHits from the MotifHitVector without freeing the vector itself.
 * @param vec Pointer to the MotifHitVector to be cleared.
 */
void freeMotifHitVector(MotifHitVector* vec);

#endif /* MOTIF_VECTOR_H */
