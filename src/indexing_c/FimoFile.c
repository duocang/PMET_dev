#include "FimoFile.h"

#include <string.h>  // For strtok
#include <stdbool.h> // For bool
#include <math.h>
#include <float.h> // for DBL_MAX

void initFimoFile_(FimoFile *file)
{
  if (!file)
  {
    fprintf(stderr, "Error: Invalid FimoFile provided.\n");
    return;
  }
  // Setting the entire structure to zero ensures that all members are `NULL` or `0`:
  memset(file, 0, sizeof(FimoFile));

  file->numLines = 0;
  file->motifName = NULL;
  file->motifLength = 0;
  file->fileName = NULL;
  file->outDir = NULL;
  file->binScore = false;
  file->hasMotifAlt = false;

  file->ht = NULL;
  file->ht = createHashTable();
  if (!file->ht)
  {
    fprintf(stderr, "Error: Memory allocation failed for Hash table *ht in initFimoFile_ function.\n");
    exit(1); // Or handle the error in another way if you prefer
  }
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
  if (file == NULL)
  {
    fprintf(stderr, "Error: Provided FimoFile pointer is NULL in initFimoFile.\n");
    return;
  }
  // Setting the entire structure to zero ensures that all members are `NULL` or `0`:
  memset(file, 0, sizeof(FimoFile));

  file->numLines = numLines;
  file->motifLength = motifLength;
  file->hasMotifAlt = hasMotifAlt;
  file->binScore = binScore;

  // 使用strdup为字符串创建新的内存副本
  // If the passed pointer is NULL, strdup will return NULL, so these operations are safe.
  // this can be different in other platforms, so we do as following
  file->motifName = motifName ? strdup(motifName) : NULL;
  file->fileName = fileName ? strdup(fileName) : NULL;
  file->outDir = outDir ? strdup(outDir) : NULL;
  // check ram allocating
  if ((motifName && !file->motifName) || (fileName && !file->fileName) || (outDir && !file->outDir))
  {
    fprintf(stderr, "Error: Memory allocation failed for strings in initFimoFile.\n");
    // Clean up any allocated memory to avoid memory leaks
    free(file->motifName);
    free(file->fileName);
    free(file->outDir);
    return;
  }

  file->ht = createHashTable();
  if (file->ht == NULL)
  {
    fprintf(stderr, "Error: Memory allocation failed for hash table in initFimoFile.\n");

    // Clean up the allocated memory to avoid memory leaks
    free(file->motifName);
    free(file->fileName);
    free(file->outDir);
    return;
  }
}


bool readFimoFile(FimoFile *fimoFile)
{
  // Ensure fimoFile and its hash table is initialized.
  if (!fimoFile || !fimoFile->ht || !fimoFile->fileName || !fimoFile->motifName || !fimoFile->outDir)
  {
    fprintf(stderr, "Error: Invalid FimoFile provided.\n");
    return false;
  }

  char *fileContent;
  long numLines = readFileAndCountLines(fimoFile->fileName, &fileContent) - 1;

  if (numLines <= 0)
  {
    free(fileContent);
    return false;
  }

  fimoFile->numLines = numLines;

  char *saveptr; // Used for strtok_r
  char *line = strtok_r(fileContent, "\n", &saveptr);
  line = strtok_r(NULL, "\n", &saveptr); // skip header
  int lineNum = 0;

  MotifHitVector *currentVec = createMotifHitVector();
  char prevGeneID[256] = "NO_GENE_YET";

  while (line)
  {
    char motif[256], motifAlt[256], geneID[256], sequence[256];
    int start, stop, binScore;
    char strand;
    double score, pval;
    MotifHit hit;
    // parse each line of file
    if (fimoFile->binScore)
    {
      if (sscanf(line, "%255s %255s %255s %d %d %c %lf %lf %255s %d", motif, motifAlt, geneID, &start, &stop, &strand, &score, &pval, sequence, &binScore) != 10)
      {
        free(fileContent);
        fprintf(stderr, "Error: Failed to parse line %d.\n", lineNum);
        return false;
      }
      initMotifHit(&hit, motif, motifAlt, geneID, start, stop, strand, score, pval, sequence, binScore);
    }
    else
    {
      if (sscanf(line, "%255s %255s %255s %d %d %c %lf %lf %255s", motif, motifAlt, geneID, &start, &stop, &strand, &score, &pval, sequence) != 9)
      {
        free(fileContent);
        fprintf(stderr, "Error: Failed to parse line %d.\n", lineNum);
        return false;
      }
      initMotifHit(&hit, motif, motifAlt, geneID, start, stop, strand, score, pval, sequence, -1);
    }
    fimoFile->motifLength = (stop - start) + 1;

    // 如果是第一行或者 geneID 发生了变化
    if (strcmp(prevGeneID, geneID) != 0)
    {
      // 如果 currentVec 不为空，则保存到哈希表
      if (currentVec && currentVec->size > 0)
      {
        putHashTable(fimoFile->ht, prevGeneID, currentVec);
      }
      // 创建新的 MotifHitVector
      currentVec = createMotifHitVector();
      strcpy(prevGeneID, geneID); // 更新prevGeneID
    }
    // 将 hit 加入到 currentVec
    pushMotifHitVector(currentVec, &hit);

    line = strtok_r(NULL, "\n", &saveptr); // Get next line using strtok_r
    lineNum++;

  } // while (line)
  // Add last MotifHitVector (gene)
  if (currentVec)
  {
    putHashTable(fimoFile->ht, prevGeneID, currentVec);
  }

  free(fileContent); // Free the content after processing
  return true;
}

void processFimoFile(FimoFile *fimoFile, int k, int N, PromoterList *promSizes)
{
  if (!fimoFile || !fimoFile->ht || !promSizes)
  {
    fprintf(stderr, "Error: Null pointer passed to processFimoFile.\n");
    exit(1);
  }

  ScoreLabelPairVector *binThresholds = createScoreLabelPairVector();

  if (!binThresholds)
  {
    fprintf(stderr, "Error: Failed to create binThresholds vector.\n");
    exit(1);
  }

  // 遍历每一个桶（bucket）
  for (int i = 0; i < TABLE_SIZE; ++i)
  {
    struct kv *current = fimoFile->ht->table[i];
    // 如果当前桶非空，遍历其链表
    while (current != NULL)
    {
      // // 这里，current->key 和 current->value 就是链表中当前节点的键和值。
      // printf("Key: %s, Value:\n", current->key);
      // printMotifHitVector(current->value);

      char *geneID = current->key;
      if (!geneID)
      {
        fprintf(stderr, "Error: Encountered a hash table with a null key.\n");
        deletePromoterLenListContent(binThresholds);
        exit(1);
      }

      MotifHitVector *vec = current->value;
      if (!vec)
      {
        fprintf(stderr, "Error: Encountered a hash table with a null value.\n");
        deletePromoterLenListContent(binThresholds);
        exit(1);
      }
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
        retainTopKMotifHits(current->value, k);
      }

      // Find the promoter size for the current gene in promSizes map
      size_t promterLength = findPromoterLength(promSizes, geneID);
      if (promterLength == -1)
      {
        printf("Error: Sequence ID: %s not found in promoter lengths file!\n", geneID);
        exit(1);
      }

      // Calculate the binomial p-value and the corresponding bin value for this gene
      Pair binom_p = geometricBinTest(current->value, promterLength, fimoFile->motifLength);

      // Save the best bin value for this gene in the binThresholds vector
      pushBack(binThresholds, binom_p.score, geneID);

      // Resize the hits for this gene if necessary based on the bin value
      // void` pointers are general-purpose pointers that can point to any
      // data type, but they don't carry type information, so you can't
      // directly dereference them or access their members
      if (((MotifHitVector *)(current->value))->size > (binom_p.idx + 1))
      {
        retainTopKMotifHits(current->value, binom_p.idx + 1);
      }

      current = current->next;
    }
  }

  // Sort the binThresholds vector by ascending score
  sortVector(binThresholds);

  // Save the Nth best bin value and gene ID to the thresholds file
  if (binThresholds->size > N)
    retainTopN(binThresholds, N);

  /****************************************************************************
   * Write PMET index result
   ****************************************************************************/
  int i;
  for (i = 0; i < binThresholds->size; i++)
  {
    char *binThresholdName = binThresholds->items[i].label;
    // printf("%s\t%f\n", binThresholdName, binThresholds->items[i].score);

    MotifHitVector *vec = getHashTable(fimoFile->ht, binThresholdName);
    writeVectorToFile(vec,
                      paste(4, "", removeTrailingSlashAndReturn(fimoFile->outDir), "/", fimoFile->motifName, ".txt"));
  }

  /****************************************************************************
   * Write "binomial_thresholds.txt"
   ****************************************************************************/
  // 打开文件以供写入binomial_thresholds.txt。如果文件不存在，将创建一个新文件。
  FILE *file = fopen(paste(3, "", removeTrailingSlashAndReturn(fimoFile->outDir), "/", "binomial_thresholds.txt"),
                     "w");
  // 检查文件是否成功打开。
  if (file == NULL)
  {
    fprintf(stderr, "Failed to open the file for writing.\n");
    exit(EXIT_FAILURE);
  }
  // return Nth best value to save in thresholds file
  double thresholdScore = binThresholds->items[binThresholds->size - 1].score;
  fprintf(file, "%s\t%f\n", fimoFile->motifName, thresholdScore);

  deleteScoreLabelVector(binThresholds);
}

Pair geometricBinTest(MotifHitVector *hitsVec, size_t promoterLength, size_t motifLength)
{
  // Null check for hitsVec
  if (!hitsVec)
  {
    fprintf(stderr, "Error: Null hitsVec provided to geometricBinTest.\n");
    exit(1); // Consider if you want to exit or handle the error differently
  }

  // Data integrity checks
  if (promoterLength <= 0 || motifLength <= 0 || motifLength > promoterLength || !hitsVec->hits)
  {
    fprintf(stderr, "Error: Invalid data provided to geometricBinTest.\n");
    exit(1); // Consider if you want to exit or handle the error differently
  }

  size_t possibleLocations = 2 * (promoterLength - motifLength + 1);

  double *pVals = (double *)malloc(hitsVec->size * sizeof(double));

  // Check memory allocation
  if (!pVals)
  {
    fprintf(stderr, "Error: Memory allocation failed in geometricBinTest.\n");
    exit(1);
  }

  for (size_t i = 0; i < hitsVec->size; i++)
  {
    pVals[i] = hitsVec->hits[i].pVal;
  }

  double lowestScore = DBL_MAX;
  size_t lowestIdx = hitsVec->size - 1;
  double product = 1.0;

  for (size_t k = 0; k < hitsVec->size; k++)
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

bool motifsOverlap(MotifHit *m1, MotifHit *m2)
{
  // Check if motifs are NULL
  if (!m1 || !m2)
  {
    return false; // Consider if you want to exit or return false in this scenario
  }

  // Check data integrity for both motifs
  if (m1->startPos > m1->stopPos || m2->startPos > m2->stopPos)
  {
    return false; // Consider if you want to exit or return false in this scenario
  }

  // Overlapping condition
  return !(m1->stopPos < m2->startPos || m1->startPos > m2->stopPos);
}

double binomialCDF(size_t numPVals, size_t numLocations, double gm)
{
  // gm is geometric mean
  if (gm <= 0.0 || gm >= 1.0 || numPVals < 0 || numPVals > numLocations)
  {
    fprintf(stderr, "Error: Invalid parameters in binomialCDF.\n");
    return -1.0; // Error value
  }

  // gm is geometric mean
  double cdf = 0.0;
  double b = 0.0;

  double logP = log(gm);
  double logOneMinusP = log(1 - gm);

  for (size_t k = 0; k < numPVals; k++)
  {
    if (k > 0)
      b += log((double)(numLocations - k + 1)) - log((double)k);
    cdf += exp(b + k * logP + (numLocations - k) * logOneMinusP);
  }
  return cdf;
}

void deleteFimoFileContents(FimoFile *file)
{
  if (!file)
    return; // Check if the provided pointer is not NULL

  // Free the hash table
  if (file->ht)
  {
    deleteHashTable(file->ht);
    file->ht = NULL;
  }

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

  printf("FimoFile (contents only) has been released!\n");
}

void deleteFimoFile(FimoFile *file)
{
  deleteFimoFileContents(file);

  // Finally, free the FimoFile struct itself
  free(file);

  printf("FimoFile (including itself) has been released!\n");
}

// Create a fake Fimo file
void createMockFimoFile(const char *fileName)
{
  FILE *file = fopen(fileName, "w");
  if (!file)
  {
    fprintf(stderr, "Error: Unable to create mock Fimo file.\n");
    return;
  }

  // a simple mockup of the motif hits. Please modify it to fit your real format.
  //  fprintf(file, "HEADER LINE\n");
  fprintf(file, "motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence\n");
  fprintf(file, "MOTIF11 MOTIF1-ALT GENE1 1 3 + 0.5 0.011 AAAC 1\n");
  fprintf(file, "MOTIF11 MOTIF2-ALT GENE1 1 4 + 0.6 0.202 TTTT 0\n");
  fprintf(file, "MOTIF11 MOTIF2-ALT GENE2 1 4 + 0.6 0.302 TTTT 0\n");
  fprintf(file, "MOTIF11 MOTIF2-ALT GENE2 1 4 + 0.6 0.602 TTTT 0\n");
  fprintf(file, "MOTIF11 MOTIF2-ALT GENE2 1 4 + 0.6 0.102 TTTT 0\n");
  fprintf(file, "MOTIF11 MOTIF1-ALT GENE2 2 3 + 0.5 0.001 AAAT 1\n");
  fprintf(file, "MOTIF11 MOTIF2-ALT GENE2 2 4 + 0.6 0.888 TTTT 0\n");
  fclose(file);
}