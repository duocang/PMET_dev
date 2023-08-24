#include "utils.h"

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

  strcpy(result, string1);

  if (spe)
  {
    strcat(result, spe);
  }

  strcat(result, string2);

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
  result[0] = '\0';

  // Second pass: concatenate strings
  va_start(args, sep);
  for (int i = 0; i < numStrings; i++)
  {
    const char *str = va_arg(args, const char *);
    strcat(result, str);
    if (i < numStrings - 1 && lenSep > 0)
    {
      strcat(result, sep);
    }
  }
  va_end(args);

  return result;
}

char *getFilenameNoExt(const char *path)
{
  const char *base = strrchr(path, '/');
  if (!base)
    base = path; // 如果没有'/'，则整个字符串是文件名
  else
    base++;

  char *dot = strrchr(base, '.');
  size_t len = dot ? (size_t)(dot - base) : strlen(base);

  char *filename = (char *)malloc(len + 1);
  if (!filename)
  {
    // 如果需要，处理内存分配失败的情况
    fprintf(stderr, "Failed to allocate memory for filename.\n");
    exit(EXIT_FAILURE);
  }

  strncpy(filename, base, len);
  filename[len] = '\0';

  return filename;
}

void removeTrailingSlash(char *path) {
    if (!path) return;  // 检查path是否为NULL

    size_t len = strlen(path);
    if (len == 0) return;  // 检查字符串是否为空

    if (path[len - 1] == '/') {
        path[len - 1] = '\0';  // 将最后的'/'替换为字符串终止符
    }
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
