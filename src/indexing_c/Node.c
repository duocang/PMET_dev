#include "Node.h"
#include <stdlib.h>
#include <string.h>

void initNodeStore(NodeStore *store)
{
  store->head = NULL;
}

Node *findNodeInStore(NodeStore *store, const char *key)
{
  Node *current = store->head;
  while (current != NULL)
  {
    if (strcmp(current->key, key) == 0)
    {
      return current;
    }
    current = current->next;
  }
  return NULL;
}

void insertIntoNodeStore(NodeStore *store, MotifHit hit)
{
  Node *node = findNodeInStore(store, hit.sequence_name);

  // If node with the key is not found, create a new one.
  if (!node)
  {
    node = malloc(sizeof(Node));
    node->key = strdup(hit.sequence_name);
    node->value = malloc(sizeof(MotifHitVector)); // Allocate memory for MotifHitVector*
    initMotifHitVector(node->value);
    node->next = store->head;
    store->head = node;
  }

  // Insert the hit into the node's vector.
  pushMotifHitVector(node->value, hit);
}

void freeNodeStore(NodeStore *store)
{
  Node *current = store->head;
  while (current)
  {
    Node *temp = current;
    free(current->key);              // Free the dynamic string
    freeMotifVector(current->value); // Free the vector's dynamic memory
    free(current->value);            // Free the MotifHitVector*
    current = current->next;         // Move to the next node
    free(temp);                      // Free the node itself
  }
  store->head = NULL; // Ensure the head of the store is NULL after all nodes are deleted
}

void printNodeStore(NodeStore *store)
{
  Node *node = store->head;
  while (node)
  {
    printf("Key: %s\n", node->key);
    for (size_t i = 0; i < node->value->size; ++i)
    {
      printMotifHit(stdout, &node->value->hits[i]);
      // printf("\n");
    }
    printf("------------------\n\n");
    node = node->next;
  }
}

long countNodesInStore(NodeStore *store)
{
  long count = 0;
  Node *current = store->head;
  while (current)
  {
    count++;
    current = current->next;
  }
  return count;
}

size_t countAllMotifHitsInStore(NodeStore *store)
{
  Node *current = store->head;
  if (current == NULL)
  {
    fprintf(stderr, "Error: The provided Node pointer is NULL.\n");
    return 0; // 返回0或其他适当的错误代码
  }

  size_t totalCount = 0; // 用于累计所有MotifHit的计数

  // 遍历链表
  while (current != NULL)
  {
    if (current->value != NULL)
    { // 检查MotifHitVector是否为NULL
      totalCount += current->value->size;
    }
    else
    {
      fprintf(stderr, "Warning: Encountered a Node with a NULL MotifHitVector. Skipping.\n");
    }
    current = current->next;
  }

  return totalCount;
}

void deleteNodeByKeyStore(NodeStore *store, const char *key)
{
  if (store == NULL || store->head == NULL)
    return;

  Node *current = store->head;
  Node *prev = NULL;

  while (current != NULL)
  {
    if (strcmp(current->key, key) == 0)
    { // 当前节点的key与给定的key匹配
      if (prev == NULL)
      { // 我们正在删除头结点
        store->head = current->next;
      }
      else
      {
        prev->next = current->next;
      }

      // 清除资源
      clearMotifHitVector(current->value); // 假设MotifHitVector有一个free函数
      free(current->key);
      free(current);
      return; // 结束删除操作
    }

    prev = current;
    current = current->next;
  }
}

void writeMotifHitsToFile(const NodeStore *store, const char *filename)
{
  FILE *file = fopen(filename, "w");
  if (file == NULL)
  {
    fprintf(stderr, "Failed to open the file for writing.\n");
    return;
  }

  Node *currentNode = store->head;
  while (currentNode)
  {
    MotifHitVector *vec = currentNode->value;
    for (size_t i = 0; i < vec->size; i++)
    {
      MotifHit hit = vec->hits[i];

      // Assuming you want to write each field of MotifHit to the file.
      // Adjust the format as needed.
      // fprintf(file, "%s\t%s\t%s\t%ld\t%ld\t%c\t%f\t%f\t%s\t%f\n",
      //         hit.motif_id,
      //         hit.motif_alt_id,
      //         hit.sequence_name,
      //         hit.startPos,
      //         hit.stopPos,
      //         hit.strand,
      //         hit.score,
      //         hit.pVal,
      //         hit.sequence,
      //         hit.binScore);
      fprintf(file, "%s\t%s\t%ld\t%ld\t%c\t%f\t%.3e\t%s\n",
              hit.motif_id,
              hit.sequence_name,
              hit.startPos,
              hit.stopPos,
              hit.strand,
              hit.score,
              hit.pVal,
              hit.sequence);
    }
    currentNode = currentNode->next;
  }

  fclose(file);
}
