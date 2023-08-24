#include "ScoreLabelPairVector.h"

ScoreLabelPairVector *createScoreLabelPairVector()
{
  ScoreLabelPairVector *vec = (ScoreLabelPairVector *)malloc(sizeof(ScoreLabelPairVector));
  vec->items = NULL;
  vec->size = 0;
  vec->capacity = 0;
  return vec;
}

int comparePairs(const void *a, const void *b)
{
  double diff = ((ScoreLabelPair *)a)->score - ((ScoreLabelPair *)b)->score;
  return (diff > 0.0) - (diff < 0.0);
}

void sortVector(ScoreLabelPairVector *vec)
{
  qsort(vec->items, vec->size, sizeof(ScoreLabelPair), comparePairs);
}

bool pushBack(ScoreLabelPairVector *vec, double score, char *label)
{
  if (vec->size >= vec->capacity)
  {
    size_t newCapacity = vec->capacity == 0 ? 4 : vec->capacity * 2;
    ScoreLabelPair *newItems = (ScoreLabelPair *)realloc(vec->items, newCapacity * sizeof(ScoreLabelPair));
    if (!newItems)
      return false;
    vec->items = newItems;
    vec->capacity = newCapacity;
  }
  vec->items[vec->size].score = score;
  // 分配内存并复制标签
  /*
      内存泄漏：你在pushBack函数里直接将传入的label指针赋值给vec->items[vec->size].label。
      如果你传入的是一个栈上的局部字符串，那么当这个字符串超出其作用域时，你的动态数组里的指针将
      会指向无效的内存。为了解决这个问题，你应该在pushBack函数里为label分配堆内存，并复制字符
      串到这片内存。
  */
  vec->items[vec->size].label = strdup(label);
  if (!vec->items[vec->size].label)
    return false;

  vec->size++;
  return true;
}

void retainTopN(ScoreLabelPairVector *vec, size_t N)
{
  // If N is greater than or equal to the current size, do nothing
  if (N >= vec->size)
  {
    return;
  }

  // Free memory for labels that are beyond the Nth element
  for (size_t i = N; i < vec->size; ++i)
  {
    free(vec->items[i].label);
  }

  // Update the size of the vector to N
  vec->size = N;

  // Resize the items array to match the new size
  vec->items = realloc(vec->items, N * sizeof(ScoreLabelPair));
  if (vec->items == NULL)
  {
    // Handle reallocation failure (for example, by exiting or returning an error)
    exit(1);
  }
  vec->capacity = N; // Update the capacity to N
}

int labelExists(const ScoreLabelPairVector* vec, const char* searchLabel) {
    for (size_t i = 0; i < vec->size; i++) {
        if (strcmp(vec->items[i].label, searchLabel) == 0) {
            return 1;  // Label found
        }
    }
    return 0;  // Label not found
}

double findScoreByLabel(ScoreLabelPairVector *vec, const char *searchLabel)
{
  for (size_t i = 0; i < vec->size; i++)
  {
    if (strcmp(vec->items[i].label, searchLabel) == 0)
    {
      return vec->items[i].score;
    }
  }
  return -1.0; // Return a special value to indicate the label was not found
}

void printVector(ScoreLabelPairVector *vec)
{
  for (size_t i = 0; i < vec->size; i++)
  {
    printf("Score: %.15f, Label: %s\n", vec->items[i].score, vec->items[i].label);
  }
}

void freeScoreLabelPairVector(ScoreLabelPairVector *vec)
{
  if (vec)
  {
    for (size_t i = 0; i < vec->size; i++)
    {
      free(vec->items[i].label);
    }
    free(vec->items);
    free(vec);
  }
}


void writeScoreLabelPairVectorToTxt(ScoreLabelPairVector* vector, const char* filename) {
    if (!vector || !filename) {
        fprintf(stderr, "Null pointer provided to saveScoreLabelPairVectorToTxt.\n");
        return;
    }

    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Unable to open file for writing");
        return;
    }

    // Assuming the first line to be headers
    fprintf(file, "Score\tLabel\n");

    for (size_t i = 0; i < vector->size; i++) {
        fprintf(file, "%.15f\t%s\n", vector->items[i].score, vector->items[i].label);
    }

    fclose(file);
}
