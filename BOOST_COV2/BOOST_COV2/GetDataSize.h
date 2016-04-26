#pragma once
#ifndef _GETDATASIZE_H
#define _GETDATASIZE_H
#include<vector>
#include "main.h"
void GetDataSize(string filename,vector<int>&DataSize,vector<vector<int>>&nCase_Gender);
void count_sample(ifstream fp_i,vector<int>&DataSize,int index,vector<vector<int>>&nCase_Gender);
#endif