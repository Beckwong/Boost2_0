// Boost2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "main.h"
#include "GetDataSize.h"
typedef long long  int64;
typedef unsigned long long uint64;
using namespace std;


class genotype{
	vector<uint64>genocase();
	vector<uint64>genoctrl();
};

class MarginalDistr{
	vector<int> MarginalDistrSNP();
	vector<vector<int> > MarginalDistrSNP_Y(3,vector<int>(2)); 
};


double Abs(double a){
	return (a<0)?-a:a;
};

int bitCount(int i){
	i = i - ((i >> 1) & 0x5555555555555555);
	i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
	i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
	i = i + (i >> 8);
	i = i + (i >> 16);
	i = i + (i >> 32);
	return (int)i & 0x7f;
}

int _tmain(int argc, _TCHAR* argv[])
{
	time_t start,end;
	vector<struct genotype*>pgeno;
	vector<struct MarginalDistr*>pMarginalDistr;

	/**** Data size in data files, [0] denotes number of samples [1]denotes number of SNPs     ***/
	vector<vector<int> >DataSize;
	vector<int>nCase_Ctrl;
	int nDataSet=0;

	int LengthLongType = sizeof(unsigned long long);
	int nsamples=0,nSNPs=0;
	int ncase, nctrl, nlongintcase, nlongintctrl;
	int icase, ictrl;

	vector<int>GenoJointDistr;
	vector<int>AllelleJointDistr;

	vector<double>MarginalEntropySNP;
	vector<double>MarginalEntropySNP_Y;
	vector<double>MarginalAssociation;

	int wordbit[65536];

	double JointEntropyTwoSNP, JointEntropyTwoSNP_Y, MarginalEntropyY;





	//start loading data

	ifstream fp,fp_i;
	string filename_i;
	string filename = "filenamelist.txt";
	fp.open(filename);

	for(int i=0;i<65536;i++){
		wordbit[i] = bitCount(i);
	}
	
	GetDataSize(filename,DataSize,nCase_Ctrl);

	/*** check output ***/
	cout<<"Number of samples:"<<DataSize[0][0]<<endl;
	cout<<"Number of SNPs in txt1:"<<DataSize[0][1]<<endl;
	cout<<"Number of samples:"<<DataSize[1][0]<<endl;
	cout<<"Number of SNPs in txt2:"<<DataSize[1][1]<<endl;
	/****           ****/

	nDataSet = DataSize.size();
	nsamples = DataSize[0][0];
	for(int i=0;i<DataSize.size();i++){
		nSNPs += DataSize[i][1];
	
	}

	nlongintcase = ceil((double)ncase/LengthLongType);
	nlongintctrl = ceil((double)nctrl/LengthLongType);



	while(1);
	return 0;
}

