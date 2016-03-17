#include "CalculateMarginalDistr.h"
#include "stdafx.h"
#include<vector>
#include"main.h"
using namespace std;

void CalculateMarginalDistr(vector<vector<genotype*>>&pgeno,int nSNPs, int nsamples, int nlongintcase, int nlongintctrl, vector<MarginalDistr*>&pMarginalDistr)
{
	int count=0;

	for(int i=0;i<nSNPs;i++)
	{
		for(int j=0;j<3;j++)
		{
			count = 0;
			for(int k=0;k<nlongintcase;k++)
			{
				count +=bitCount(pgeno[i][j]->genocase[k]);
			}
			pMarginalDistr[i]->MarginalDistrSNP_Y[j][0] = count;
			
			count = 0;
			for(int k=0;k<nlongintctrl;k++)
			{
				count+= bitCount(pgeno[i][j]->genoctrl[k]);
			}
			pMarginalDistr[i]->MarginalDistrSNP_Y[j][1] = count;

			pMarginalDistr[i]->MarginalDistrSNP[j] = pMarginalDistr[i]->MarginalDistrSNP_Y[j][0] + pMarginalDistr[i]->MarginalDistrSNP_Y[j][1];
		}
	}
}