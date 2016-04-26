#ifndef _CALCULATEGENOJOINTDISTR_H
#define _CALCULATEGENOJOINTDISTR_H
#include"main.h"

using namespace std;

void CalculateGenoJointDistr(vector<vector<genotype*>>&pgeno,vector<vector<int>>&GenoJointDistr,int snp1,int snp2,vector<MarginalDistr*>&pMarginalDistr);

#endif