#pragma once
#ifndef _CALCULATEMARGINALENTROPY_H
#define _CALCULATEMARGINALENTROPY_H

#include<iostream>
#include "main.h"
#include<vector>

void CalculateMarginalEntropy(vector<vector<genotype*>>&pgeno,vector<vector<int> >&DataSize,vector<int>&nCase_Ctrl,vector<double>&MarginalEntropySNP,vector<double>MarginalEntropySNP_Y);

#endif