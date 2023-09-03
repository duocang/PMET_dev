#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "HashTable.h"
#include "MotifHit.h"
#include "MotifHitVector.h"

// 要放入哈希表中的结构体
struct Student
{
  int age;
  float score;
  char name[32];
  char data[1024 * 1024 * 10];
};

// 结构体内存释放函数
static void free_student(void *stu)
{
  free(stu);
}

// 显示学生信息的函数
static void show_student(struct Student *p)
{
  printf("姓名:%s, 年龄:%d, 学分:%.2f\n", p->name, p->age, p->score);
}

void testStudent()
{
  // 新建一个HashTable实例
  HashTable *ht = createHashTable();
  if (NULL == ht)
  {
    printf("Fail to create a hash table!\n");
  }

  // 向哈希表中加入多个学生结构体
  for (int i = 0; i < 100; i++)
  {
    struct Student *stu = (struct Student *)malloc(sizeof(struct Student));
    stu->age = 18 + rand() % 5;
    stu->score = 50.0f + rand() % 100;
    sprintf(stu->name, "同学%d", i);
    putHashTable2(ht, stu->name, stu, free_student);
  }

  // 根据学生姓名查找学生结构
  for (int i = 0; i < 10; i++)
  {
    char name[32];
    sprintf(name, "同学%d", i);
    struct Student *stu = (struct Student *)getHashTable(ht, name);
    show_student(stu);
  }

  // 销毁哈希表实例
  deleteHashTable(ht);
}


void testMotif()
{
  MotifHit hit1, hit2, hit3, hit4;
  initMotifHit(&hit1, "AHL12", "孙悟空", "西游记", 614, 621, '+', 7.8501, 0.000559, "AAATAATT", 0);
  initMotifHit(&hit2, "AHL12", "唐三藏", "西游记", 700, 708, '-', 8.5001, 0.000600, "AAGGTTAA", 1);
  initMotifHit(&hit3, "AHL12", "诸葛亮", "三国志", 800, 808, '+', 6.7000, 0.000700, "TTAACCAA", 2);
  initMotifHit(&hit4, "AHL12", "刘皇叔", "三国志", 650, 658, '-', 7.1000, 0.000800, "GGGTTTCC", 3);

  HashTable *ht = createHashTable();
  if (NULL == ht)
  {
    printf("Fail to create a hash table!\n");
  }

  putHashTable(ht, "西游记", &hit1);
  putHashTable(ht, "西游记", &hit2);
  putHashTable(ht, "三国志", &hit3);
  putHashTable(ht, "三国志", &hit4);


  MotifHit *hit = getHashTable(ht, "西游记");

  printMotifHit(stdout, hit);

  deleteHashTable(ht);
}


void freeVec(void *ptr) {
    MotifHitVector *vec = (MotifHitVector *)ptr;
    if(vec) {
        deleteMotifHitVectorContent(vec);  // 释放 hits 数组
        free(vec);        // 释放 MotifHitVector 结构
    }
}

void testMotifVector()
{
  MotifHit hit1, hit2, hit3, hit4;
  initMotifHit(&hit1, "AHL12", "孙悟空", "西游记", 614, 621, '+', 7.8501, 0.000559, "AAATAATT", 0);
  initMotifHit(&hit2, "AHL12", "唐三藏", "西游记", 700, 708, '-', 8.5001, 0.000600, "AAGGTTAA", 1);
  initMotifHit(&hit3, "AHL12", "诸葛亮", "三国志", 800, 808, '+', 6.7000, 0.000700, "TTAACCAA", 2);
  initMotifHit(&hit4, "AHL12", "刘皇叔", "三国志", 650, 658, '-', 7.1000, 0.000800, "GGGTTTCC", 3);

  // MotifHitVector vec1, vec2;
  MotifHitVector vec1;
  MotifHitVector *vec2 = createMotifHitVector();
  initMotifHitVector(&vec1);
  // initMotifHitVector(&vec2);

  pushMotifHitVector(&vec1, &hit1);
  pushMotifHitVector(&vec1, &hit2);
  pushMotifHitVector(vec2, &hit3);
  pushMotifHitVector(vec2, &hit4);

  HashTable *ht = createHashTable();
  if (NULL == ht)
  {
    printf("Fail to create a hash table!\n");
  }

  putHashTable2(ht, "西游记", &vec1, NULL);
  putHashTable2(ht, "三国志", vec2, freeVec);

  printf("提取三国志：\n");
  MotifHitVector *vec = getHashTable(ht, "三国志");
  printMotifHitVector(vec);

  deleteHashTable(ht);
}





int main()
{
  // testStudent();
  // testMotif();
  testMotifVector();
  return 0;
}

// clang -o test TestHashTable.c HashTable.c MotifHit.c MotifHitVector.c