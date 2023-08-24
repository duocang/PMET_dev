#include <stdio.h>

#include "FileRead.h"
#include "MotifHit.h"
#include "MotifHitVector.h"
#include "Node.h"
#include "PromoterLength.h"

#include "FimoFile.h"
int main() {
    FimoFile* myFimoFile = malloc(sizeof(FimoFile));

    // 使用参数初始化FimoFile结构体
    initFimoFile(myFimoFile,
                 10,                         // numLines
                 "MYB46_2",                 // motifName
                 0,                         // motifLength
                 "test_data/MYB46_2_short.txt",   // fileName
                 "./test_result",                   // outDir
                 false,                     // hasMotifAlt
                 false                      // binScore
    );

    // 读取fimo.txt内容
    if(!readFimoFile(myFimoFile)) {
        fprintf(stderr, "Error reading fimo.txt!\n");
        freeFimoFile(myFimoFile);
        return 1;
    }

    printNodeStore(myFimoFile->nodeStore);

    // 打印出读取到的MotifHit信息
    Node *currentNode = myFimoFile->nodeStore->head;
    while(currentNode) {
        printf("Sequence ID: %s\n", currentNode->key);
        printMotifHitVector(currentNode->value);
        currentNode = currentNode->next;
    }



    PromoterList* list = malloc(sizeof(PromoterList));
    readPromoterLengthFile(list, "test_data/promoter_lengths.txt");
    // printf("Length of AT1G01010: %ld\n", findPromoterLength(list, "AT1G01010"));


    processFimoFile(myFimoFile, 5, 1000, list);

    currentNode = myFimoFile->nodeStore->head;
    while(currentNode) {
        printf("Sequence ID: %s\n", currentNode->key);
        printMotifHitVector(currentNode->value);
        currentNode = currentNode->next;
    }

    freePromoterList(list);
    // 清理FimoFile中的内存
    freeFimoFile(myFimoFile);

    return 0;
}

// clang -o test TestFimoFile_small.c FileRead.c Node.c MotifHit.c MotifHitVector.c FimoFile.c PromoterLength.c ScoreLabelPairVector.c utils.c