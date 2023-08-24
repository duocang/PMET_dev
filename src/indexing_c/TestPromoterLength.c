#include "PromoterLength.h"

int main() {

    PromoterList* list = malloc(sizeof(PromoterList));
    initPromoterList(list);

    readPromoterLengthFile(list, "test_data/promoter_lengths.txt");

    printf("Length of AT1G01010: %d\n", findPromoterLength(list, "AT1G01010"));
    // ... you can test with other gene names

    freePromoterList(list);
    return 0;

    // PromoterList list;
    // initPromoterList(&list);

    // readPromoterLengthFile(&list, "test_data/promoter_lengths.txt");

    // printf("Length of AT1G01010: %d\n", findPromoterLength(&list, "AT1G01010"));
    // // ... you can test with other gene names

    // freePromoterList(&list);
    // return 0;
}
// clang TestPromoterLength.c PromoterLength.c -o test



// int main()
// {
//   HashTable ht;
//   initHashTable(&ht);

//   readDataFromFile(&ht, "test_data/promoter_lengths.txt");

//   printf("Length of AT1G01010: %d\n", find(&ht, "AT1G01010"));
//   // ... you can test with other gene names

//   freeHashTable(&ht);
//   return 0;
// }

