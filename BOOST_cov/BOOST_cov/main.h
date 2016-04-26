#ifndef _MAIN_H
#define _MAIN_H
#include<iostream>
#include<fstream>
#include <vector>
#include<math.h>
#include<string>
#include<stdlib.h>

typedef long long int64;
typedef unsigned long long uint64;


using namespace std;

int popcount(uint64 i);
double Abs(double a);
int bitCount(uint64 i);

class genotype{
public:
	vector<vector<uint64>>geno;

	genotype():geno(4,vector<uint64>()){};
};

class MarginalDistr{
public:
	vector<int>MarginalDistrSNP;
	vector<vector<int>>MarginalDistrSNP_Y_G;
	MarginalDistr():MarginalDistrSNP(3),MarginalDistrSNP_Y_G(4,vector<int>(3)){};
};
#endif