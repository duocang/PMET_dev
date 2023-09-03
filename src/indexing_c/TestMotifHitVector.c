#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "MotifHit.h"
#include "MotifHitVector.h"

int main()
{
  // 创建三个MotifHit数据
  MotifHit hit1, hit2, hit3;

  initMotifHit(&hit1, "AHL12", "孙悟空", "AT1G01020", 614, 621, '+', 7.85401, 0.000559, "AAATAATT", 0);
  initMotifHit(&hit2, "AHL13", "猪八戒", "AT1G01020", 715, 722, '-', 6.85401, 0.009959, "TTAATAAT", 1);
  initMotifHit(&hit3, "AHL14", "沙和尚", "AT1G01020", 816, 823, '+', 5.85401, 0.000759, "GGGTAAGG", 2);
  initMotifHit(&hit3, "AHL14", "唐三藏", "AT1G01020", 816, 823, '+', 5.85401, 0.000009, "GGGTAAGG", 2);

  // 初始化MotifVector
  MotifHitVector vec;
  initMotifHitVector(&vec);

  // 添加到vector
  pushMotifHitVector(&vec, &hit1);
  pushMotifHitVector(&vec, &hit2);
  pushMotifHitVector(&vec, &hit3);

  // 确保数据正确地添加到了vector
  assert(strcmp(vec.hits[0].motif_id, "AHL12") == 0);
  assert(strcmp(vec.hits[1].motif_id, "AHL13") == 0);
  assert(strcmp(vec.hits[2].motif_id, "AHL14") == 0);

  // Test the new functions
  printf("Vector size: %zu\n", motifHitVectorSize(&vec)); // Expected output: 3

  printf("Printing MotifHits in the vector:\n");
  printMotifHitVector(&vec);

  printf("Testing ordering function:\n");
  // 对vector中的MotifHit按pVal排序
  sortMotifHitVectorByPVal(&vec);
  printMotifHitVector(&vec);

  printf("Testing TopK function:\n");
  retainTopKMotifHits(&vec, 2);
  printMotifHitVector(&vec);


  printf("Testing delete function:\n");
  removeHitAtIndex(&vec, 1);
  printMotifHitVector(&vec);


  printf("All tests passed for MotifVector!\n");

  // 清理
  freeMotifHit(&hit1);
  freeMotifHit(&hit2);
  freeMotifHit(&hit3);
  deleteMotifHitVectorContent(&vec);

  return 0;
}

// clang MotifHit.c MotifHitVector.c TestMotifHitVector.c -o test