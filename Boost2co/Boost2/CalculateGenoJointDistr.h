#pragma once
#ifndef _CALCULATEGENOJOINTDISTR_H
#define _CALCULATEGENOJOINTDISTR_H

#include "main.h"
#include<vector>

void CalculateGenoJointDistr(vector<vector<genotype*>>&pgeno,int nSNPs,int nlongintcase, int nlongintctrl,vector<vector<int>> &GenoJointDistr,int snp1,int snp2,vector<MarginalDistr*>&pMarginalDistr);






#endif