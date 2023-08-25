#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

char *paste2(const char *spe, const char *string1, const char *string2)
{
  if (!string1 || !string2)
    return NULL;

  int len1 = strlen(string1);
  int len2 = strlen(string2);
  spe = spe ? spe : ""; // If spe is NULL, treat it as an empty string
  int lenSep = strlen(spe);

  char *result = malloc(len1 + len2 + lenSep + 1); // +1 for the null-terminator
  if (!result)
  {
    perror("Memory allocation failed");
    exit(EXIT_FAILURE);
  }

  // Copy the first string
  strcpy(result, string1);
  // Copy the separator
  strcpy(result + len1, spe);
  // Copy the second string
  strcpy(result + len1 + lenSep, string2);

  return result;
}

char *paste(int numStrings, const char *sep, ...)
{
  va_list args;
  size_t length = 0;
  sep = sep ? sep : ""; // If sep is NULL, treat it as an empty string
  int lenSep = strlen(sep);

  // First pass: calculate the total length
  va_start(args, sep);
  for (int i = 0; i < numStrings; i++)
  {
    const char *str = va_arg(args, const char *);
    length += strlen(str);
    if (i < numStrings - 1)
    {
      length += lenSep;
    }
  }
  va_end(args);

  // Allocate memory for the result
  char *result = (char *)malloc(length + 1);
  if (!result)
  {
    perror("Memory allocation failed");
    exit(EXIT_FAILURE);
  }
  char *currentPos = result; // Pointer to track our position in the result string

  // Second pass: concatenate strings
  va_start(args, sep);
  for (int i = 0; i < numStrings; i++)
  {
    const char *str = va_arg(args, const char *);
    size_t lenStr = strlen(str);

    memcpy(currentPos, str, lenStr);
    currentPos += lenStr;

    if (i < numStrings - 1 && lenSep > 0)
    {
      memcpy(currentPos, sep, lenSep);
      currentPos += lenSep;
    }
  }
  va_end(args);

  *currentPos = '\0'; // Null-terminate the result string

  return result;
}

char* removeTrailingSlashAndReturn(const char *path) {
    if (!path) return NULL;

    size_t len = strlen(path);
    if (len == 0) return strdup(""); // 返回一个空字符串的复制品

    char *newPath = strdup(path); // 复制原始字符串
    if (!newPath) return NULL; // 如果内存分配失败

    if (newPath[len - 1] == '/') {
        newPath[len - 1] = '\0';  // 删除最后的 '/'
    }

    return newPath;
}
