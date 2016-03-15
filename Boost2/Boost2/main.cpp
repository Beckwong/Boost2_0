// Boost2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "main.h"
#include "LoadData.h"
#include "GetDataSize.h"
#include "CalculateMarginalEntropy.h"
#include "CalculateMarginalDistr.h"
#include "CalculateGenoJointDistr.h"
typedef long long  int64;
typedef unsigned long long uint64;
using namespace std;

static unsigned int wordbit[65536];

static int popCount(uint64 i)
{
	return (wordbit[i&0xFFFF]+wordbit[(i>>16)&0xFFFF]+wordbit[(i>>32)&0xFFFF]+wordbit[i>>48]);
}

double Abs(double a){
	return (a<0)?-a:a;
}

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
	vector<vector<genotype*>>pgeno;
	vector<MarginalDistr*>pMarginalDistr;

	/**** Data size in data files, [0] denotes number of samples [1]denotes number of SNPs     ***/
	vector<vector<int> >DataSize;
	vector<int>nCase_Ctrl;
	int nDataSet=0;

	
	int nsamples=0,nSNPs=0;
	int nlongintcase, nlongintctrl;


	vector<double>MarginalEntropySNP;
	vector<double>MarginalEntropySNP_Y;
	vector<double>MarginalAssociation;

	vector<vector<int>> GenoJointDistr(2,vector<int>(9));
	vector<vector<double>>Pab(3,vector<double>(3,0));
	vector<vector<double>>Pbc(2,vector<double>(3,0));
	vector<vector<double>>Pca(2,vector<double>(3,0));
	vector<pair<int,int>>InteractionPairs;
	vector<double>InteractionMeasurePairs;
	vector<int>DistrCollection(1001);
	int InteractionCount = 0;
	double maxInteraction = - 9999999;
	double minInteraction =  9999999;
	double thresholdRecord = 30;

	

	double JointEntropyTwoSNP, JointEntropyTwoSNP_Y, MarginalEntropyY,ptmp1,ptmp2;





	//start loading data

	ifstream fp;
	string filename = "filenamelist.txt";

	for(int i=0;i<65536;i++)
	{
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

	//memory allocation
	nlongintcase = ceil((double)nCase_Ctrl[0]/LengthLongType);
	nlongintctrl = ceil((double)nCase_Ctrl[1]/LengthLongType);

	/*pgeno.resize(nSNPs);
	for(int i=0;i<nSNPs;i++){
		pgeno[i].resize(3);
	}*/
	
	for(int i=0;i<nSNPs;i++)
	{
		vector<genotype*>tmp;
		for(int j=0;j<3;j++)
		{
			genotype *tmp_genotype = new genotype;
			tmp.push_back(tmp_genotype);
		}
		pgeno.push_back(tmp);
	}
	

	for(int i=0;i<nSNPs;i++)
	{
		for(int j=0;j<3;j++)
		{
			pgeno[i][j]->genocase.resize(nlongintcase);
			pgeno[i][j]->genoctrl.resize(nlongintctrl);
		}
	}	


	//load data to bit representation
	LoadData(filename,pgeno,DataSize);
	cout<<endl;
	cout<<hex<<pgeno[0][0]->genocase[0]<<endl;
	cout<<hex<<pgeno[0][1]->genocase[0]<<endl;
	cout<<hex<<pgeno[0][2]->genocase[0]<<endl;


	//calculate marginal distribution
	ptmp1=(double)nCase_Ctrl[0]/(nCase_Ctrl[0]+nCase_Ctrl[1]);
	MarginalEntropyY = -ptmp1 * log(ptmp1)-(1-ptmp1) * log(1-ptmp1);
	
	MarginalEntropySNP.resize(nSNPs);
	MarginalEntropySNP_Y.resize(nSNPs);
	MarginalAssociation.resize(nSNPs);

	CalculateMarginalEntropy(pgeno,DataSize,nCase_Ctrl,MarginalEntropySNP,MarginalEntropySNP_Y);
	for( int i=0;i<nSNPs;i++)
	{
		MarginalAssociation[i] = (-MarginalEntropySNP_Y[i] + MarginalEntropySNP[i] + MarginalEntropyY)*nsamples*2;
	}


	for(int i=0;i<nSNPs;i++)
	{
		MarginalDistr* tmp=new MarginalDistr();
		pMarginalDistr.push_back(tmp);
	} 
	  
	CalculateMarginalDistr(pgeno,nSNPs,nsamples,nlongintcase,nlongintctrl,pMarginalDistr);

	DistrCollection.resize(nSNPs);

	//calculating interaction 
	for(int snp1=0;snp1<nSNPs-1;snp1++)
	{
		for(int snp2=snp1+1;snp2<nSNPs;snp2++)
		{
			//calculate joint distribution
			CalculateGenoJointDistr(pgeno,nSNPs,nlongintcase,nlongintctrl,GenoJointDistr,snp1,snp2,pMarginalDistr);

			//P(A|B) transpose
			/*index
				B=0  B=1  B=2
			A=0  0    1    2
			A=1  3    4    5
			A=2  6    7    8
			*/

			//i for B(snp2), j for A(snp1)
			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					Pab[i][j] = (double)(GenoJointDistr[0][3*j+i]+GenoJointDistr[1][3*j+i])/pMarginalDistr[snp2]->MarginalDistrSNP[i];
				}
			}
			//P(B|C)
			/*
				index 
				C=0    C=1
			B=0  0       3
			B=1  1       4
			B=2  2       5
			*/
			//i for C(case or control), j for B(snp2)

			for(int i=0;i<2;i++)
			{
				for(int j=0;j<3;j++)
				{
					Pbc[i][j] = (double)(pMarginalDistr[snp2]->MarginalDistrSNP_Y[j][i])/nCase_Ctrl[i];
				}
			}

			//P(C|A)
			/*index
					A=0		A=1		A=2
			C = 0    0		2		4
			C = 1	 1		3		5
			*/
			//i for A, j for C
			for(int i=0;i<3;i++)
			{
				for(int j=0;j<2;j++)
				{
					Pca[i][j] = (double)(pMarginalDistr[snp1]->MarginalDistrSNP_Y[i][j])/pMarginalDistr[snp1]->MarginalDistrSNP[i];
				}
			}

			//Papprx = Pab*Pbc*Pca
			// tao = sum(Papprx)
			// sum(p.* log(p) - p.*log(Pappr) + p.* log(tao))
			double tao = 0.0;
			double InteractionMeasure;
			    //i for A,j for B, k for C
			for(int i = 0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					for(int k=0;k<2;k++)
					{
						ptmp1 = (double)GenoJointDistr[k][i*3+j]/nsamples;
						if(ptmp1>0)
						{
							InteractionMeasure += ptmp1 * log(ptmp1);  
						}
						ptmp2 = Pab[j][i]*Pbc[k][j]*Pca[i][k];
						if(ptmp2>0)
						{
							InteractionMeasure += -ptmp1*log(ptmp2);
							tao+=ptmp2;
						}
					}
				}
				
			}

			InteractionMeasure = (InteractionMeasure + log(tao))*nsamples*2;
			
			maxInteraction = max(InteractionMeasure,maxInteraction);
			minInteraction = min(InteractionMeasure,minInteraction);

			if(InteractionMeasure > thresholdRecord)
			{
				InteractionPairs.push_back(make_pair(snp1,snp2));
				InteractionMeasurePairs.push_back(InteractionMeasure);
			}

			if(InteractionMeasure>100)
			{
				DistrCollection[1000]++;
			}
			else if(InteractionMeasure>0)
			{
				DistrCollection[(int)(InteractionMeasure/0.1)]++;
			}
		}
		if((snp1+1)%100 == 0)
		{
			cout<<"Iteration "<<snp1+1<<endl;
		}
	}
    
	//post correction
	cout<<"Start post correction: Excat Gtest..."<<endl;

	for(int i=0;i<InteractionMeasurePairs.size();i++)
	{
		int snp1 = InteractionPairs[i].first;

	}

	 
	while(1);
	return 0;
}

