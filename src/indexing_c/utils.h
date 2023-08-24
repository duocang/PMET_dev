#ifndef PMET_INDEX_UTILS_H
#define PMET_INDEX_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

char* paste2(const char* sep, const char* string1, const char* string2);
char* paste(int numStrings, const char* sep, ...);
char* getFilenameNoExt(const char* path);
void removeTrailingSlash(char *path);
char* removeTrailingSlashAndReturn(const char *path);

#endif /* PMET_INDEX_UTILS_H */
