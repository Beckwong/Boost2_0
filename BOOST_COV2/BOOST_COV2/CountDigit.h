#pragma once
#ifndef _COUNTDIGIT_H
#define _COUNTDIGIT_H

static unsigned int wordbit[65536];
int bitCount(uint64 i);
int popcount(uint64 i);
void CountDigit();
#endif