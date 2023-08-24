#ifndef PROMOTER_LENGTH_H
#define PROMOTER_LENGTH_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PROMOTER_NAME_LENGTH 100  // Define a macro for maximum promoter name length


typedef struct Promoter {
    char* promoterName;
    int length;
    struct Promoter* next;
} Promoter;

typedef struct {
    Promoter* head;
} PromoterList;


void initPromoterList(PromoterList* list);
long findPromoterLength(PromoterList* list, const char* promoterName);
void freePromoterList(PromoterList* list);
void readPromoterLengthFile(PromoterList* list, const char* filename);






// #include "FileRead.h"

// #define TABLE_SIZE 30007 // Use a prime number for better distribution

// typedef struct
// {
//   char *promoterName;
//   int length;
// } HashEntry;

// typedef struct
// {
//   HashEntry *table[TABLE_SIZE];
// } HashTable;

// void initHashTable(HashTable *ht);
// unsigned int hash(const char *str);
// void readDataFromFile(HashTable* ht, const char* filename);
// void insert(HashTable *ht, const char *promoterName, int length);
// int find(HashTable *ht, const char *promoterName);
// void freeHashTable(HashTable *ht);

#endif /* PROMOTER_LENGTH_H */

