#ifndef FIMOFILE_H
#define FIMOFILE_H

// Standard library headers
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

// Custom headers specific to this project
#include "utils.h"
#include "ScoreLabelPairVector.h"
#include "PromoterLength.h"
#include "FileRead.h"
#include "MotifHit.h"
#include "MotifHitVector.h"
// #include "Node.h"
#include "HashTable.h"
#include "MemCheck.h"

/**
 * Struct representing a simple pair with an index and a score.
 */
typedef struct
{
  long idx;
  double score;
} Pair;

/**
 * Struct to hold motif names and associated bin scores.
 */
typedef struct
{
  char *motifName;
  double binScore;
} btVals;

/**
 * Struct representing the size details of a promoter.
 */
typedef struct
{
  char *promoter;
  long len;
} promSizes;

/**
 * Dynamic string array that can grow in size.
 */
typedef struct
{
  char **items;
  int size;
  int capacity;
} DynamicArray;

/**
 * Represents details and metadata about a Fimo file.
 */
typedef struct
{
  int numLines;
  char *motifName;
  int motifLength;
  char *fileName;
  char *outDir;
  bool hasMotifAlt;  // Indicates if there's an alternative motif present.
  bool binScore;     // Indicates if bin scores are used.
  HashTable *ht;
} FimoFile;

/**
 * Create an FimoFile with default values.
 */
FimoFile* createFimoFile();

/**
 * Initializes a FimoFile structure with the provided values.
 *
 * @param file        A pointer to the FimoFile structure to be initialized.
 * @param numLines    The number of lines in the FIMO file.
 * @param motifName   The name of the motif.
 * @param motifLength The length of the motif.
 * @param fileName    The name of the FIMO file.
 * @param outDir      The output directory for FIMO results.
 * @param hasMotifAlt Boolean flag indicating if there's an alternate motif.
 * @param binScore    Boolean flag indicating if binomial scoring is to be used.
 */
void initFimoFile(FimoFile *file,
                  int numLines,
                  char *motifName,
                  int motifLength,
                  char *fileName,
                  char *outDir,
                  bool hasMotifAlt,
                  bool binScore);

// /**
//  * Copies the contents of a source FimoFile into a destination FimoFile.
//  * @param source The source FimoFile.
//  * @param dest The destination FimoFile.
//  * @return Boolean indicating success of the operation.
//  */
// bool copyFimoFile(const FimoFile *source, FimoFile *dest);

/**
 * Read the contents of a FimoFile.
 * @param file The FimoFile to read from.
 * @return Boolean indicating success of the reading operation.
 */
bool readFimoFile(FimoFile *file);

/**
 * Determine if two motifs overlap.
 * @param m1 First Motif.
 * @param m2 Second Motif.
 * @return Boolean indicating if the motifs overlap.
 */
bool motifsOverlap(MotifHit *m1, MotifHit *m2);

/**
 * Run geometric bionomial test on a vector of motif hits.
 * @param hitsVec Vector containing motif hits.
 * @param ... Additional parameters used in the test.
 * @return Pair result of the test.
 */
Pair geometricBinTest(MotifHitVector *hitsVec, size_t promoterLength, size_t motifLength);

/**
 * Calculate the cumulative distribution function for a binomial distribution.
 * probability of X < numPVals
 * @param numPVals Number of p-values.
 * @param numLocations Number of locations.
 * @param gm Geometric mean.
 * @return Double value of the cumulative distribution function.
 * @example
 * binomialCDF(1, 22, 0.01) = 0.801631
 * Probability of X = 1 events	0.17814013100868 (17.81%)
 * Probability of X ≤ 1 events	0.97977072054772 (97.98%)
 * Probability of X > 1 events	0.02022927945228 (2.02%)
 * Probability of X < 1 events	0.80163058953905 (80.16%)
 * Probability of X ≥ 1 events	0.19836941046095 (19.84%)
 */
double binomialCDF(size_t numPVals, size_t numLocations, double gm);

/**
 * Frees all memory allocated for a FimoFile' content only.
 * @param file The FimoFile to free.
 */
void deleteFimoFileContents(FimoFile *file);



/**
 * Frees all memory allocated for a FimoFile.
 * @param file The FimoFile to free.
 */
void deleteFimoFile(FimoFile *file);

/**
 * Process the contents of a FimoFile and potentially output some results.
 * @param fimoFile The FimoFile to process.
 * @param k Parameter for processing.
 * @param N Parameter for processing.
 * @param promSizes List containing promoter sizes.
 */
void processFimoFile(FimoFile *fimoFile, int k, int N, PromoterList *promSizes);

/**
 * Create a mock fimo file.
 * @param fileName The FimoFile to save.
 */
void createMockFimoFile(const char *fileName);

#endif /* FIMOFILE_H */
