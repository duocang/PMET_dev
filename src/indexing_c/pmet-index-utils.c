#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

char* paste2(const char* spe, const char* string1, const char* string2) {
    if (!string1 || !string2) return NULL;

    int len1 = strlen(string1);
    int len2 = strlen(string2);
    spe = spe ? spe : "";  // If spe is NULL, treat it as an empty string
    int lenSep = strlen(spe);

    char* result = malloc(len1 + len2 + lenSep + 1);  // +1 for the null-terminator
    if (!result) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    strcpy(result, string1);

    if (spe) {
        strcat(result, spe);
    }

    strcat(result, string2);

    return result;
}


char* paste(int numStrings, const char* sep, ...) {
    va_list args;
    size_t length = 0;
    sep = sep ? sep : "";  // If sep is NULL, treat it as an empty string
    int lenSep = strlen(sep);

    // First pass: calculate the total length
    va_start(args, sep);
    int i = 0;
    for (i = 0; i < numStrings; i++) {
        const char* str = va_arg(args, const char*);
        length += strlen(str);
        if (i < numStrings - 1) {
            length += lenSep;
        }
    }
    va_end(args);

    // Allocate memory for the result
    char* result = (char*)malloc(length + 1);
    if (!result) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
    result[0] = '\0';

    // Second pass: concatenate strings
    va_start(args, sep);
    int i = 0;
    for (i = 0; i < numStrings; i++) {
        const char* str = va_arg(args, const char*);
        strcat(result, str);
        if (i < numStrings - 1 && lenSep > 0) {
            strcat(result, sep);
        }
    }
    va_end(args);

    return result;
}
