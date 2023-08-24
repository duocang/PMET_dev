#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

char* paste2(const char* sep, const char* string1, const char* string2);
char* paste(int numStrings, const char* sep, ...);