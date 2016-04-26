#pragma once
#ifndef _MAIN_H
#define _MAIN_H
#include<iostream>
#include<fstream>
#include <vector>
#include<math.h>
#include<string>

typedef long long int64;
typedef unsigned long long uint64;

class genotype{
public:
	vector<vector<uint64>>geno;

	genotype():geno(4,vector<uint64>()){};
};
#endif