#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

int main() {
    char* result = paste(5, "-", "Hello", "world", "this", "is", "C", NULL);
    printf("%s\n", result); // Outputs: Hello-world-this-is-C
    free(result);

    result = paste(5, NULL, "Hello", "world", "this", "is", "C");
    printf("%s\n", result); // Outputs: Hello-world-this-is-C
    free(result);

    result = paste(5, "", "Hello", "world", "this", "is", "C");
    printf("%s\n", result); // Outputs: Hello-world-this-is-C
    free(result);



    result = paste2("Hello", "world", "-");
    printf("%s\n", result); // Outputs: Hello-world-this-is-C
    free(result);

    result = paste2("", "Hello", "world");
    printf("%s\n", result); // Outputs: Hello-world-this-is-C
    free(result);

    result = paste2(NULL, "Hello", "world");
    printf("%s\n", result); // Outputs: Hello-world-this-is-C
    free(result);
    return 0;
}


// clang -o test TestUtils.c utils.c
