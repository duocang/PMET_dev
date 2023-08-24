#ifndef MOTIF_HIT_VECTOR_NODE_H
#define MOTIF_HIT_VECTOR_NODE_H

#include "MotifHitVector.h"

// Node structure for the linked list based store.
typedef struct Node
{
  char *key;
  MotifHitVector *value;
  struct Node *next;
} Node;

typedef struct
{
  Node *head;
} NodeStore;

// Initializes a KeyValueStore.
void initNodeStore(NodeStore *store);
// Finds a node with the given key in the store.
Node *findNodeInStore(NodeStore *store, const char *key);
// Inserts a MotifHit into the store.
void insertIntoNodeStore(NodeStore *store, MotifHit hit);
// Frees memory allocated for the KeyValueStore.
void freeNodeStore(NodeStore *store);
void printNodeStore(NodeStore *store);
long countNodesInStore(NodeStore *store);
size_t countAllMotifHitsInStore(NodeStore *store);
void deleteNodeByKeyStore(NodeStore *store, const char *key);
void writeMotifHitsToFile(const NodeStore* store, const char* filename);

#endif /* MOTIF_HIT_VECTOR_NODE_H */
