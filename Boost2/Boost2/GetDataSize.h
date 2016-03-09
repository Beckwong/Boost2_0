#pragma once
#ifndef _GETDATASIZE_H
#define _GETDATASIZE_H

#include<vector>
#include<string>
#include<time.h>
#include<fstream>
#include<iostream>

using namespace std;

void GetDataSize(string filename,vector<vector<int>>&DataSize,vector<int>&nCase_Ctrl);
void count_sample(ifstream& fp_i,vector<vector<int> >&DataSize,int& index,vector<int>&nCase_Ctrl);
#endif