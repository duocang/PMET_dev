#ifndef FIMOFILE_H
#define FIMOFILE_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "utils.h"
#include "ScoreLabelPairVector.h"
#include "PromoterLength.h"
#include "FileRead.h"
#include "MotifHit.h"
#include "MotifHitVector.h"
#include "Node.h"

typedef struct {
    long idx;
    double score;
} Pair;

typedef struct {
  char* motifName;
  double binScore;
} btVals;

typedef struct {
  char* promoter;
  long len;
} promSizes;

typedef struct {
    char** items;
    int size;
    int capacity;
} DynamicArray;


// FimoFile结构体
typedef struct
{
  int numLines;
  char *motifName;
  int motifLength;
  char *fileName;
  char *outDir;
  bool hasMotifAlt;
  bool binScore;
  NodeStore*  nodeStore;
} FimoFile;

// // Function prototypes
// FimoFile* FimoFile_create();
void initFimoFile_(FimoFile *file);
void initFimoFile(FimoFile *file,
                  int numLines,
                  char *motifName,
                  int motifLength,
                  char *fileName,
                  char *outDir,
                  bool hasMotifAlt,
                  bool binScore);
bool copyFimoFile(const FimoFile *source, FimoFile *dest);

bool readFimoFile(FimoFile *file);
bool motifsOverlap(MotifHit *m1, MotifHit *m2);
Pair geometricBinTest( MotifHitVector *hitsVec,  long promoterLength,  long motifLength);
double binomialCDF(long numPVals, long numLocations, double gm);
// void FimoFile_process(FimoFile *file, long k, long N, const char *promSizesKey, long promSizesValue, char **outString, double *outDouble);
// void FimoFile_destroy(FimoFile *file);
void freeFimoFile(FimoFile *file);

// btVals processFimoFile(FimoFile *fimoFile, int k, int N, PromoterList* promSizes);
void processFimoFile(FimoFile *fimoFile, int k, int N, PromoterList* promSizes);

#endif /* FIMOFILE_H */
