#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "ScoreLabelPairVector.h"

int main() {
    ScoreLabelPairVector* vec = createScoreLabelPairVector();

    // Add items to the vector
    pushBack(vec, 3.4, "C");
    pushBack(vec, 1.2, "A");
    pushBack(vec, 2.7, "B");
    pushBack(vec, 4.8, "D");

    printf("Before sorting:\n");
    printVector(vec);

    // Sort and print
    sortVector(vec);

    printf("\nAfter sorting:\n");
    printVector(vec);

    // Find score by label
    char* searchLabel = "B";
    double foundScore = findScoreByLabel(vec, searchLabel);
    if (foundScore != -1.0) {
        printf("\nScore for label %s: %f\n", searchLabel, foundScore);
    } else {
        printf("\nLabel %s not found.\n", searchLabel);
    }

    // Clean up
    deletePromoterLenListContent(vec);
    // for (size_t i = 0; i < vec->size; i++) {
    //     free(vec->items[i].label); // If dynamically allocated in a real scenario
    // }
    // free(vec->items);
    // free(vec);
    return 0;
}

// clang TestScoreLabelPairVector.c ScoreLabelPairVector.c -o test
