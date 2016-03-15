#pragma once
#ifndef _CALCULATEMARGINALDISTR_H
#define _CALCULATEMARGINALDISTR_H

#include "main.h"

#include<vector>

void CalculateMarginalDistr(vector<vector<genotype*>>&pgeno,int nSNPs, int nsamples, int nlongintcase, int nlongintctrl, vector<MarginalDistr*>&pMarginalDistr);

#endif