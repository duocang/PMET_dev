#include <stdio.h>
#include "MotifHit.h"
#include "MotifHitVector.h"
#include "Node.h"

int main()
{
  NodeStore store;
  initNodeStore(&store);

  MotifHit hit1, hit2, hit3, hit4;
  initMotifHit(&hit1, "AHL12", "孙悟空", "AT1G01020", 614, 621, '+', 7.85401, 0.000559, "AAATAATT", 0);
  initMotifHit(&hit2, "AHL15", "唐三藏", "AT1G01021", 700, 708, '-', 8.5001, 0.000600, "AAGGTTAA", 1);
  initMotifHit(&hit3, "AHL18", "猪八戒", "AT1G01022", 800, 808, '+', 6.7000, 0.000700, "TTAACCAA", 2);
  initMotifHit(&hit4, "AHL20", "沙僧", "AT1G01020", 650, 658, '-', 7.1000, 0.000800, "GGGTTTCC", 3);

  // Add the test data to the store
  insertIntoNodeStore(&store, hit1);
  insertIntoNodeStore(&store, hit2);
  insertIntoNodeStore(&store, hit3);
  insertIntoNodeStore(&store, hit4);

  // Print the store content to check if the data has been added correctly.
  printNodeStore(&store);




  long genesNum = countNodesInStore(&store);
  size_t hitsNum  = countAllMotifHitsInStore(&store);
  printf("%ld genes and %ld hits found related to %s\n\n", genesNum, hitsNum, "ALH15");





  printf("Testing deleteNodeByKeyStore\n\n");
  deleteNodeByKeyStore(&store, "AT1G01022");
  printNodeStore(&store);


  printf("Testing writeMotifHitsToFile\n\n");
  writeMotifHitsToFile(&store, "test_result/TestNode_resut.txt");


  genesNum = countNodesInStore(&store);
  hitsNum  = countAllMotifHitsInStore(&store);
  printf("Aftering filtering..\n");
  printf("%ld genes and %ld hits found related to %s\n\n", genesNum, hitsNum, "ALH15");




  // Free allocated memory
  freeMotifHit(&hit1);
  freeMotifHit(&hit2);
  freeMotifHit(&hit3);
  freeNodeStore(&store);

  return 0;
}

// clang -o test TestNode.c Node.c MotifHit.c MotifHitVector.c