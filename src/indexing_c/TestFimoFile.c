#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "FileRead.h"
#include "MotifHit.h"
#include "MotifHitVector.h"
#include "Node.h"
#include "PromoterLength.h"

#include "FimoFile.h"

void testInitFimoFile_()
{
  // Scenario 1: Normal initialization of a FimoFile
  FimoFile *file = malloc(sizeof(FimoFile));
  if (!file)
  {
    fprintf(stderr, "Error: Memory allocation failed for FimoFile in testInitFimoFile_.\n");
    exit(1);
  }
  initFimoFile_(file);

  assert(file->numLines == 0);
  assert(file->motifName == NULL);
  assert(file->motifLength == 0);
  assert(file->fileName == NULL);
  assert(file->outDir == NULL);
  assert(file->binScore == false);
  assert(file->hasMotifAlt == false);
  assert(file->nodeStore != NULL);
  assert(file->nodeStore->head == NULL);

  freeFimoFileContents(file); // Assuming you have a function that can free the resources held by FimoFile

  // Scenario 2: Pass a NULL FimoFile to initFimoFile_
  // This should print an error message but not crash.
  initFimoFile_(NULL);

  printf("All tests passed for initFimoFile_!\n\n\n");
}

void testInitFimoFile()
{
  // 创建一个FimoFile对象的指针
  FimoFile *testFile = malloc(sizeof(FimoFile));
  if (!testFile)
  {
    fprintf(stderr, "Failed to allocate memory for testFile.\n");
    exit(1);
  }

  // 调用初始化函数
  initFimoFile(testFile, 10, "motif_test", 5, "file_test", "dir_test", true, true);

  // 进行断言测试
  assert(testFile->numLines == 10);
  assert(testFile->motifLength == 5);
  assert(testFile->hasMotifAlt == true);
  assert(testFile->binScore == true);

  assert(strcmp(testFile->motifName, "motif_test") == 0);
  assert(strcmp(testFile->fileName, "file_test") == 0);
  assert(strcmp(testFile->outDir, "dir_test") == 0);

  freeFimoFile(testFile);

  printf("All tests for initFimoFile passed!\n\n\n");
}

void testReadFimoFile()
{
  const char *mockFileName = "test_data/mock_fimo.txt";
  createMockFimoFile(mockFileName);

  FimoFile testFile;
  printf("Value of file->nodeStore: %p\n", testFile.nodeStore);

  initFimoFile_(&testFile);

  testFile.fileName = strdup(mockFileName);
  testFile.motifName = strdup(mockFileName);
  testFile.outDir = strdup(mockFileName);
  testFile.binScore = true; // 请根据模拟文件内容进行设置

  assert(readFimoFile(&testFile) == true);
  assert(testFile.numLines == 2);
  assert(testFile.motifLength == 3);

  freeFimoFileContents(&testFile);

  printf("All tests in testReadFimoFile passed!\n");
}

// void testCopyFimoFile() {
//   const char *mockFileName = "test_data/mock_fimo.txt";
//   createMockFimoFile(mockFileName);

//   // Initialize source FimoFile
//   FimoFile sourceFile;
//   initFimoFile_(&sourceFile);
//   sourceFile.fileName = strdup(mockFileName);
//   sourceFile.binScore = true;  // 根据模拟文件内容进行设置
//   assert(readFimoFile(&sourceFile) == true);

//   // Copy to destination FimoFile
//   FimoFile destFile;
//   initFimoFile_(&destFile);
//   assert(copyFimoFile(&sourceFile, &destFile) == true);

//   // Check if the files are equal
//   assert(areFimoFilesEqual(&sourceFile, &destFile));

//   // Cleanup
//   freeFimoFile(&sourceFile);
//   freeFimoFile(&destFile);
// }


int main()
{
  testInitFimoFile_();
  testInitFimoFile();
  testReadFimoFile();

  FimoFile *myFimoFile = malloc(sizeof(FimoFile));

  // 使用参数初始化FimoFile结构体
  initFimoFile(myFimoFile,
               10,                    // numLines
               "CCA1",               // motifName
               0,                     // motifLength
               "test_data/CCA1.txt", // fileName
               "./test_result",       // outDir
               false,                 // hasMotifAlt
               false                  // binScore
  );

  if (!readFimoFile(myFimoFile))
  {
    fprintf(stderr, "Error reading fimo.txt!\n");
    freeFimoFile(myFimoFile);
    return 1;
  }

  Node *currentNode = myFimoFile->nodeStore->head;

  PromoterList *list = malloc(sizeof(PromoterList));
  // initPromoterList(list);
  readPromoterLengthFile(list, "test_data/promoter_lengths.txt");
  printf("Length of AT1G01010: %ld\n", findPromoterLength(list, "AT1G01010"));

  // printNodeStore(myFimoFile->nodeStore);

  processFimoFile(myFimoFile, 5, 5000, list);

  long genesNum = countNodesInStore(myFimoFile->nodeStore);
  size_t hitsNum = countAllMotifHitsInStore(myFimoFile->nodeStore);
  printf("\nAftering filtering..\n");
  printf("%ld genes and %ld hits found related to %s\n\n", genesNum, hitsNum, myFimoFile->motifName);

  // printNodeStore(myFimoFile->nodeStore);

  freePromoterList(list);
  freeFimoFile(myFimoFile);

  return 0;
}

// clang -o test TestFimoFile.c FileRead.c Node.c MotifHit.c MotifHitVector.c FimoFile.c PromoterLength.c ScoreLabelPairVector.c utils.c