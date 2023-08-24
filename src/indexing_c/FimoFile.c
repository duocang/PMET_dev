#include "FimoFile.h"

#include <string.h>  // For strtok
#include <stdbool.h> // For bool
#include <math.h>
#include <float.h> // for DBL_MAX

void initFimoFile_(FimoFile *file)
{
  file->numLines = 0;
  file->motifName = NULL;
  file->motifLength = 0;
  file->fileName = NULL;
  file->outDir = NULL;
  file->binScore = false;
  file->hasMotifAlt = false;
  initNodeStore(file->nodeStore);
}

void initFimoFile(FimoFile *file,
                  int numLines,
                  char *motifName,
                  int motifLength,
                  char *fileName,
                  char *outDir,
                  bool hasMotifAlt,
                  bool binScore)
{
  file->numLines = numLines;
  file->motifLength = motifLength;
  file->hasMotifAlt = hasMotifAlt;
  file->binScore = binScore;

  // 使用strdup为字符串创建新的内存副本
  // 如果传递的指针为NULL，strdup将返回NULL，所以这些操作是安全的
  file->motifName = strdup(motifName);
  file->fileName = strdup(fileName);
  file->outDir = strdup(outDir);

  // 初始化NodeStore
  // 注意：在此我假定`NodeStore`结构体有一个对应的初始化函数
  // 如果没有，请确保正确初始化file->nodeStore
  file->nodeStore = malloc(sizeof(NodeStore));
  if (file->nodeStore)
    initNodeStore(file->nodeStore); // 假设你有这样的函数
}

bool copyFimoFile(const FimoFile *source, FimoFile *dest)
{
  if (!source || !dest)
  {
    return false;
  }

  initFimoFile_(dest);

  dest->numLines = source->numLines;
  dest->motifLength = source->motifLength;
  dest->hasMotifAlt = source->hasMotifAlt;

  if (source->motifName)
  {
    dest->motifName = strdup(source->motifName);
    if (!dest->motifName)
    {
      return false; // Memory allocation failed.
    }
  }

  if (source->fileName)
  {
    dest->fileName = strdup(source->fileName);
    if (!dest->fileName)
    {
      return false; // Memory allocation failed.
    }
  }

  if (source->outDir)
  {
    dest->outDir = strdup(source->outDir);
    if (!dest->outDir)
    {
      return false; // Memory allocation failed.
    }
  }

  // Here, we should also copy the content of nodes.
  // Depending on your specific needs, you might want to deep copy the list or just copy the pointers.

  // Assuming a deep copy:
  Node *currentSrcNode = source->nodeStore->head;
  while (currentSrcNode)
  {
    for (size_t i = 0; i < currentSrcNode->value->size; i++)
    {
      // 将每个MotifHit插入到目标的NodeStore中
      insertIntoNodeStore(dest->nodeStore, currentSrcNode->value->hits[i]);
    }
    currentSrcNode = currentSrcNode->next;
  }

  return true;
}

bool readFimoFile(FimoFile *fimoFile)
{
  char *fileContent;
  long numLines = readFileAndCountLines(fimoFile->fileName, &fileContent);

  if (numLines <= 0)
  {
    free(fileContent);
    return false;
  }

  fimoFile->numLines = numLines;

  // Now tokenize the content line by line
  char *line = strtok(fileContent, "\n");
  line = strtok(NULL, "\n"); // skip header
  int lineNum = 1;
  while (line)
  {
    char motif[256], motifAlt[256], geneID[256], sequence[256];
    int start, stop, binScore;
    char strand;
    double score, pval;

    int itemsParsed;
    MotifHit hit;
    if (fimoFile->binScore)
    {
      itemsParsed = sscanf(line, "%255s %255s %255s %d %d %c %lf %lf %255s %d",
                           motif, motifAlt, geneID, &start, &stop, &strand, &score, &pval, sequence, &binScore);
      // Use the parsed values to initialize the hit.
      initMotifHit(&hit, motif, motifAlt, geneID, start, stop, strand, score, pval, sequence, binScore);
    }
    else
    {
      itemsParsed = sscanf(line, "%255s %255s %255s %d %d %c %lf %lf %255s",
                           motif, motifAlt, geneID, &start, &stop, &strand, &score, &pval, sequence);
      // Use the parsed values to initialize the hit.
      initMotifHit(&hit, motif, motifAlt, geneID, start, stop, strand, score, pval, sequence, -1);
    }

    fimoFile->motifLength = (stop - start) + 1;

    // Insert the hit into the NodeStore of the FimoFile
    insertIntoNodeStore(fimoFile->nodeStore, hit);

    line = strtok(NULL, "\n"); // Get next line

    // printf("Reading line %d out of %ld lines\n", lineNum, numLines);
    lineNum++;
  }

  long genesNum = countNodesInStore(fimoFile->nodeStore);
  printf("%ld genes and %ld hits found related to %s\n\n", genesNum, numLines, fimoFile->motifName);

  free(fileContent); // Free the content after processing
  return true;
}

void processFimoFile(FimoFile *fimoFile, int k, int N, PromoterList *promSizes)
{
  // printf("%ld\n", countNodesInStore(fimoFile->nodeStore));
  ScoreLabelPairVector *binThresholds = createScoreLabelPairVector();

  Node *currentNode = fimoFile->nodeStore->head;
  long numDone = 0;

  while (currentNode)
  {
    char *geneID = currentNode->key;
    long geneNum = countNodesInStore(fimoFile->nodeStore);
    MotifHitVector* vec = currentNode->value;

    // Sort hits within the gene based on their p-values.
    sortMotifHitVectorByPVal(vec);

    // Top k motifHit in vector with overlap with any one in the vector
    size_t currentIndex = 0;
    while (currentIndex < vec->size && currentIndex < k)
    {
      size_t nextIndex = currentIndex + 1;

      while (nextIndex < vec->size)
      {
        if (motifsOverlap(&vec->hits[currentIndex], &vec->hits[nextIndex]))
        {
          removeHitAtIndex(vec, nextIndex);
          // Do not increment nextIndex here because after removing
          // an element, the next element shifts to the current nextIndex
        }
        else
        {
          nextIndex++; // No overlap, move to next hit
        }
      }
      currentIndex++;
    }
    if (vec->size > k)
    {
      retainTopKMotifHits(currentNode->value, k);
    }

    // if (strcmp(geneID, "AT1G20680") == 0) {
    //   printMotifHitVector(currentNode->value);
    // }

    // Find the promoter size for the current gene in promSizes map
    long promterLength = findPromoterLength(promSizes, geneID);
    if (promterLength == -1)
    {
      printf("Error: Sequence ID: %s not found in promoter lengths file!\n", geneID);
      exit(1);
    }

    // Calculate the binomial p-value and the corresponding bin value for this gene
    Pair binom_p = geometricBinTest(currentNode->value, promterLength, fimoFile->motifLength);

    // if (strcmp(geneID, "AT1G20680") == 0) {
    //   printf(" binom_p.score：%f\n",  binom_p.score);

    //   printf(" binom_p.idx：%ld\n",  binom_p.idx);

    //   printf(" currentNode->value->size：%d\n",  currentNode->value->size);
    // }
    // Save the best bin value for this gene in the binThresholds vector
    pushBack(binThresholds, binom_p.score, geneID);

    // Resize the hits for this gene if necessary based on the bin value
    if (currentNode->value->size > (binom_p.idx + 1))
    {
      retainTopKMotifHits(currentNode->value, binom_p.idx + 1);
    }

    currentNode = currentNode->next; // go for next gene
  } // while (currentNode)

  // Sort the binThresholds vector by ascending score
  sortVector(binThresholds);

  // Save the Nth best bin value and gene ID to the thresholds file
  if (binThresholds->size > N)
    retainTopN(binThresholds, N);
  // printVector(binThresholds);
  writeScoreLabelPairVectorToTxt(binThresholds, paste(3, "", fimoFile->outDir, "/",  "binomial_thresholds.txt"));

  // printf("countNodesInStore(fimoFile->nodeStore): %ld\n", countNodesInStore(fimoFile->nodeStore));

  DynamicArray genesDeleted;
  genesDeleted.items = malloc(sizeof(char*) * 10); // initial capacity, say 10
  genesDeleted.size = 0;
  genesDeleted.capacity = 10;

  currentNode = fimoFile->nodeStore->head;
  while (currentNode)
  {
    char *geneID = currentNode->key;

    // Check if the gene is in binThresholds.
    if (labelExists(binThresholds, geneID) != 1)
    {
      // Resize if necessary
      if (genesDeleted.size == genesDeleted.capacity) {
          genesDeleted.capacity *= 2; // Double the capacity
          genesDeleted.items = realloc(genesDeleted.items, sizeof(char*) * genesDeleted.capacity);
      }
      // Append geneID
      genesDeleted.items[genesDeleted.size] = strdup(geneID);
      genesDeleted.size++;
    }
    currentNode = currentNode->next; // go for next gene
  } // while (currentNode)

  // 遍历之前存储在动态数组中的geneID
  int i = 0;
  for (i = 0; i < genesDeleted.size; i++) {
      char *geneIDToDelete = genesDeleted.items[i];
      deleteNodeByKeyStore(fimoFile->nodeStore, geneIDToDelete);
  }
  // Release the memory of a dynamic array.
  for (i = 0; i < genesDeleted.size; i++) {
    free(genesDeleted.items[i]); // Free each strdup'ed string
  }
  free(genesDeleted.items);


  // write the bin thresholds to file
  char* outputDir = removeTrailingSlashAndReturn(fimoFile->outDir);
  printf("Write all processed fimo results to %s\n", paste(4, "", outputDir, "/", fimoFile->motifName, ".txt"));
  writeMotifHitsToFile(fimoFile->nodeStore, paste(4, "", outputDir, "/", fimoFile->motifName, ".txt"));

  free(outputDir);
  freeScoreLabelPairVector(binThresholds); // Release the memory of a dynamic array.
}

Pair geometricBinTest(MotifHitVector *hitsVec, long promoterLength, long motifLength)
{
  long possibleLocations = 2 * (promoterLength - motifLength + 1);

  double *pVals = (double *)malloc(hitsVec->size * sizeof(double));

  for (size_t i = 0; i < hitsVec->size; i++)
  {
    pVals[i] = hitsVec->hits[i].pVal; // Assuming MotifHit has a member named pVal
  }

  double lowestScore = DBL_MAX;
  long lowestIdx = hitsVec->size - 1;
  double product = 1.0;

  for (long k = 0; k < hitsVec->size; k++)
  {
    product *= pVals[k];
    double geom = pow(product, 1.0 / (k + 1.0));
    double binomP = 1 - binomialCDF(k + 1, possibleLocations, geom);

    if (lowestScore > binomP)
    {
      lowestScore = binomP;
      lowestIdx = k;
    }
  }

  free(pVals);

  Pair result;
  result.idx = lowestIdx;
  result.score = lowestScore;

  return result;
}

bool motifsOverlap(MotifHit* m1, MotifHit* m2)
{
  return !(m2->startPos > m1->stopPos || m2->stopPos < m1->startPos);
}

double binomialCDF(long numPVals, long numLocations, double gm)
{
  // gm is geometric mean
  double cdf = 0.0;
  double b = 0.0;

  double logP = log(gm);
  double logOneMinusP = log(1 - gm);

  for (long k = 0; k < numPVals; k++)
  {
    if (k > 0)
      b += log((double)(numLocations - k + 1)) - log((double)k);
    cdf += exp(b + k * logP + (numLocations - k) * logOneMinusP);
  }
  return cdf;
}

void freeFimoFile(FimoFile *file)
{
  // Free the nodes store
  freeNodeStore(file->nodeStore);

  // Free motifName if allocated
  if (file->motifName)
  {
    free(file->motifName);
    file->motifName = NULL;
  }

  // Free fileName if allocated
  if (file->fileName)
  {
    free(file->fileName);
    file->fileName = NULL;
  }

  // Free outDir if allocated
  if (file->outDir)
  {
    free(file->outDir);
    file->outDir = NULL;
  }
}
