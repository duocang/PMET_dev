#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "FileRead.h"
#include "MotifHit.h"
#include "MotifHitVector.h"
#include "PromoterLength.h"

#include "FimoFile.h"
#include "HashTable.h"
#include "MemCheck.h"

void testInitFimoFile_()
{
  // Scenario 1: Normal initialization of a FimoFile
  printf("Testing ininatization:\n");
  printf("\n    Scenario 1: Normal initialization of a FimoFile:");
  // FimoFile *file = malloc(sizeof(FimoFile));
  FimoFile *file;
  file = (FimoFile *)new_malloc(sizeof(FimoFile));
  // file = (FimoFile *)malloc(sizeof(FimoFile));
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
  assert(file->ht != NULL);

  deleteFimoFileContents(file); // Assuming you have a function that can free the resources held by FimoFile
  new_free(file);

  // Scenario 2: Pass a NULL FimoFile to initFimoFile_
  // This should print an error message but not crash.
  printf("\n    Scenario 2: Pass a NULL FimoFile to initFimoFile_:\n");
  initFimoFile_(NULL);

  // Scenario 3: No free of the resources held by FimoFile
  FimoFile *file2;
  file2 = (FimoFile *)new_malloc(sizeof(FimoFile));


  printf("All tests passed for initFimoFile_!\n\n\n\n");
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

  assert(testFile->ht != NULL);

  deleteFimoFile(testFile);

  printf("All tests for initFimoFile passed!\n\n\n");
}

void testReadFimoFile()
{
  // const char *mockFileName = "test_data/mock_fimo.txt";
  // createMockFimoFile(mockFileName);
  const char *mockFileName = "test_data/MYB46_2_short.txt";

  FimoFile testFile;
  printf("Value of file->ht: %p\n", testFile.ht);

  initFimoFile_(&testFile);

  testFile.fileName = strdup(mockFileName);
  testFile.motifName = strdup(mockFileName);
  testFile.outDir = strdup(mockFileName);
  testFile.binScore = false; // 请根据模拟文件内容进行设置

  assert(readFimoFile(&testFile) == true);
  assert(testFile.numLines == 8);
  assert(testFile.motifLength == 8);

  MotifHitVector *vec = getHashTable(testFile.ht, "AT1G01020");
  printMotifHitVector(vec);

  vec = getHashTable(testFile.ht, "AT1G01010");
  printMotifHitVector(vec);

  printf("All tests in testReadFimoFile passed!\n");
}

void testProcess()
{
  printf("Testing testProcess...\n");
  FimoFile *myFimoFile = malloc(sizeof(FimoFile));

  // 使用参数初始化FimoFile结构体
  initFimoFile(myFimoFile,
               0,                     // numLines
               "MYB52",               // motifName
               0,                     // motifLength
               "test_data/MYB52.txt", // fileName
               "./test_result",       // outDir
               false,                 // hasMotifAlt
               false                  // binScore
  );

  if (!readFimoFile(myFimoFile))
  {
    fprintf(stderr, "Error reading fimo.txt!\n");
    deleteFimoFile(myFimoFile);
  }

  MotifHitVector *vec = getHashTable(myFimoFile->ht, "AT1G01010");
  printMotifHitVector(vec);

  HashTable *ht = myFimoFile->ht;
  // for (size_t i = 0; i < TABLE_SIZE; i++)
  // {
  //   struct kv *current = myFimoFile->ht->table[i];
  //   while (current != NULL)
  //   {
  //     printf("Key: %s, Value:\n", current->key);
  //     printMotifHitVector(current->value);
  //     current = current->next;
  //   }
  // }

  PromoterList *list = malloc(sizeof(PromoterList));
  // initPromoterList(list);
  readPromoterLengthFile(list, "test_data/promoter_lengths.txt");
  printf("Length of AT1G01010: %ld\n", findPromoterLength(list, "AT1G01010"));

  processFimoFile(myFimoFile, 5, 5000, list);

  deletePromoterLenListContent(list);
  deleteFimoFile(myFimoFile);

  printf("All tests in testReadFimoFile passed!\n");
}

int main()
{
#ifdef DEBUG
  atexit(show_block); // 在程序结束后显示内存泄漏报告
#endif
  testInitFimoFile_();
  // testInitFimoFile();
  // testReadFimoFile();

  // testProcess();

  return 0;
}

/*
clang -DDEBUG \
      -o test \
      TestFimoFile.c  \
      FileRead.c  \
      HashTable.c  \
      MemCheck.c  \
      Node.c  \
      MotifHit.c  \
      MotifHitVector.c  \
      FimoFile.c  \
      PromoterLength.c  \
      ScoreLabelPairVector.c  \
      utils.c

*/
