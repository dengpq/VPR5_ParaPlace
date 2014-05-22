#ifndef READLINE_H
#define READLINE_H

#include "util.h"

char** ReadLineTokens(INOUT FILE* InFile,
                      INOUT int* LineNum);
int CountTokens(IN char** Tokens);
void FreeTokens(INOUT char** *TokensPtr);

#endif
