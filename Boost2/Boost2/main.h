#pragma once
#ifndef _MAIN_H
#define _MAIN_H

#include<vector>
#include<math.h>
#include<time.h>
#include<vector>
#include<time.h>
#include<fstream>
#include<iostream>
#include<algorithm>



using namespace std;
typedef long long  int64;
typedef unsigned long long uint64;

static int LengthLongType = 64;

static int popcount(uint64 i);
double Abs(double a);
int bitCount(int i);

class genotype{
public:
	vector<uint64>genocase;
	vector<uint64>genoctrl;

	genotype(){};
};

class MarginalDistr{
public:
	vector<int> MarginalDistrSNP;
	vector<vector<int>>MarginalDistrSNP_Y; 
	MarginalDistr()£ºMarginalDistrSNP(3),MarginalDistrSNP_Y(3,vector<int>(2)){};
};
#endif