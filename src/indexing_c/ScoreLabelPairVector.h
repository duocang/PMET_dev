#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    double score;
    char* label;
} ScoreLabelPair;

typedef struct {
    ScoreLabelPair* items;
    size_t size;
    size_t capacity;
} ScoreLabelPairVector;


ScoreLabelPairVector* createScoreLabelPairVector();
bool pushBack(ScoreLabelPairVector* vec, double score, char* label);
void retainTopN(ScoreLabelPairVector* vec, size_t N);
int comparePairs(const void* a, const void* b);
int labelExists(const ScoreLabelPairVector* vec, const char* searchLabel);
double findScoreByLabel(ScoreLabelPairVector* vec, const char* searchLabel);
void sortVector(ScoreLabelPairVector* vec);
void printVector(ScoreLabelPairVector* vec);
void freeScoreLabelPairVector(ScoreLabelPairVector* vec);
void writeScoreLabelPairVectorToTxt(ScoreLabelPairVector* vector, const char* filename);
