#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "PromoterLength.h"


void initPromoterList(PromoterList* list) {
    list->head = NULL;
}

void insertPromoter(PromoterList* list, const char* promoterName, int length) {
    Promoter* newPromoter = (Promoter*)malloc(sizeof(Promoter));
    if (!newPromoter) {
        perror("Failed to allocate memory for new promoter");
        return;
    }
    newPromoter->promoterName = strdup(promoterName);
    if (!newPromoter->promoterName) {
        perror("Failed to allocate memory for promoter name");
        free(newPromoter);
        return;
    }
    newPromoter->length = length;
    newPromoter->next = list->head;
    list->head = newPromoter;
}


long findPromoterLength(PromoterList* list, const char* promoterName) {
    Promoter* current = list->head;
    while (current) {
        if (strcmp(current->promoterName, promoterName) == 0) {
            return current->length;
        }
        current = current->next;
    }
    return -1; // Not found
}


void freePromoterList(PromoterList* list) {
    Promoter* current = list->head;
    while (current) {
        Promoter* toDelete = current;
        free(toDelete->promoterName);
        current = current->next;
        free(toDelete);
    }
    list->head = NULL;
}

void readPromoterLengthFile(PromoterList* list, const char* filename) {
    // 初始化PromoterList
    initPromoterList(list);
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        return;
    }

    char promoterName[MAX_PROMOTER_NAME_LENGTH];
    long length;
    while (fscanf(file, "%99s %ld", promoterName, &length) == 2) {  // Read max 99 chars to prevent overflow
        if (strlen(promoterName) >= MAX_PROMOTER_NAME_LENGTH - 1) {
            fprintf(stderr, "Promoter name too long: %s\n", promoterName);
            continue;  // Skip this line and move to next
        }
        insertPromoter(list, promoterName, length);
    }

    fclose(file);
}




// void initHashTable(HashTable *ht)
// {
//  int i = 0
//   for (i = 0; i < TABLE_SIZE; i++)
//   {
//     ht->table[i] = NULL;
//   }
// }

// unsigned int hash(const char *str)
// {
//   unsigned int value = 0;
//   for (; *str != '\0'; str++)
//   {
//     value = (value << 5) + *str;
//   }
//   return value % TABLE_SIZE;
// }

// void readDataFromFile(HashTable *ht, const char *filename)
// {
//   char *fileContent;
//   long numLines = readFileAndCountLines(filename, &fileContent);

//   if (numLines <= 0)
//   {
//     free(fileContent);
//     printf("Error: Cannot open file %s\n", filename);
//     exit(1);
//   }

//   // Now tokenize the content line by line
//   char *line = strtok(fileContent, "\n");
//   int itemsParsed = 0;

//   while (line)
//   {
//     char promoterName[255];
//     long length;
//     itemsParsed = sscanf(line, "%255s %ld", promoterName, &length);
//     insert(ht, promoterName, length);
//     // printf("%s %ld\n", promoterName, length);

//     line = strtok(NULL, "\n"); // Get next line
//   }

//   free(fileContent); // Free the content after processing
// }

// void insert(HashTable *ht, const char *promoterName, int length)
// {
//   unsigned int slot = hash(promoterName);
//   while (ht->table[slot] != NULL && strcmp(ht->table[slot]->promoterName, promoterName) != 0)
//   {
//     slot = (slot + 1) % TABLE_SIZE; // Linear probing
//   }
//   if (ht->table[slot] != NULL)
//   {
//     free(ht->table[slot]->promoterName); // Release the old key
//   }
//   else
//   {
//     ht->table[slot] = malloc(sizeof(HashEntry));
//   }
//   ht->table[slot]->promoterName = strdup(promoterName);
//   ht->table[slot]->length = length;
// }

// int find(HashTable *ht, const char *promoterName)
// {
//   unsigned int slot = hash(promoterName);
//   while (ht->table[slot] != NULL && strcmp(ht->table[slot]->promoterName, promoterName) != 0)
//   {
//     slot = (slot + 1) % TABLE_SIZE;
//   }
//   if (ht->table[slot] == NULL)
//     return -1;
//   return ht->table[slot]->length;
// }

// void freeHashTable(HashTable *ht)
// {
//   for (int i = 0; i < TABLE_SIZE; i++)
//   {
//     if (ht->table[i] != NULL)
//     {
//       free(ht->table[i]->promoterName);
//       free(ht->table[i]);
//     }
//   }
// }
