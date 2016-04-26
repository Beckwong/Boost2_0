#include"CountDigit.h"
#include"main.h"


void CountDigit()
{
	for(uint64 i =0;i<65535;i++)
	{
		wordbit[i]=bitCount(i);
	}
}

int bitCount(uint64 i)
{
	i = i - ((i >> 1) & 0x5555555555555555);
	i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
	i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
	i = i + (i >> 8);
	i = i + (i >> 16);
	i = i + (i >> 32);
	return (int)i & 0x7f;
}

int popcount(uint64 i)
{
	return (wordbit[i&0xFFFF]+wordbit[(i>>16)&0xFFFF]+wordbit[(i>>32)&0xFFFF]+wordbit[i>>48]);
}