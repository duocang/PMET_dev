#include "PromoterLength.h"

int main() {

    PromoterList* list = malloc(sizeof(PromoterList));
    initPromoterList(list);

    readPromoterLengthFile(list, "test_data/promoter_lengths.txt");

    printf("Length of AT1G01010: %zu\n", findPromoterLength(list, "AT1G01010"));
    // ... you can test with other gene names

    deletePromoterLenListContent(list);
    return 0;
}
// clang TestPromoterLength.c PromoterLength.c -o test
