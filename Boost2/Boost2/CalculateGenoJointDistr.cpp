#include "CalculateGenoJointDistr.h"
#include "stdafx.h"
#include<vector>
#include "main.h"
using namespace std;


void CalculateGenoJointDistr(vector<vector<genotype*>>&pgeno,int nSNPs,int nlongintcase, int nlongintctrl,vector<vector<int>> &GenoJointDistr,int snp1,int snp2,vector<MarginalDistr*>&pMarginalDistr)
{
	register int count;

	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			count = 0;

			for(int k=0;k<nlongintcase;k++)
			{
				count += popcount(pgeno[snp1][i]->genocase[k] & pgeno[snp2][j]->genocase[k]);
			}

			GenoJointDistr[0][i*3+j] = count;
			count = 0;

			for(int k=0;k<nlongintctrl;k++)
			{
				count+= popcount(pgeno[snp1][i]->genoctrl[k] & pgeno[snp2][j]->genoctrl[k]);
			}

			GenoJointDistr[1][i*3+j] = count;

		}
	}

	//for case;
	GenoJointDistr[0][2] = pMarginalDistr[snp1]->MarginalDistrSNP_Y[0][0] - GenoJointDistr[0][0] - GenoJointDistr[0][1];
	GenoJointDistr[0][5] = pMarginalDistr[snp1]->MarginalDistrSNP_Y[1][0] - GenoJointDistr[0][3] - GenoJointDistr[0][4];
	
	GenoJointDistr[0][6] = pMarginalDistr[snp2]->MarginalDistrSNP_Y[0][0] - GenoJointDistr[0][0] - GenoJointDistr[0][3];
	GenoJointDistr[0][7] = pMarginalDistr[snp2]->MarginalDistrSNP_Y[1][0] - GenoJointDistr[0][1] - GenoJointDistr[0][4];
	
	GenoJointDistr[0][8] = pMarginalDistr[snp2]->MarginalDistrSNP_Y[2][0] - GenoJointDistr[0][2] - GenoJointDistr[0][5];

	// for control:
	GenoJointDistr[1][2] = pMarginalDistr[snp1]->MarginalDistrSNP_Y[0][1] - GenoJointDistr[1][0] - GenoJointDistr[1][1];
	GenoJointDistr[1][5] = pMarginalDistr[snp1]->MarginalDistrSNP_Y[1][1] - GenoJointDistr[1][3] - GenoJointDistr[1][4];
	
	GenoJointDistr[1][6] = pMarginalDistr[snp2]->MarginalDistrSNP_Y[0][1] - GenoJointDistr[1][0] - GenoJointDistr[1][3];
	GenoJointDistr[1][7] = pMarginalDistr[snp2]->MarginalDistrSNP_Y[1][1] - GenoJointDistr[1][1] - GenoJointDistr[1][4];
	
	GenoJointDistr[1][8] = pMarginalDistr[snp1]->MarginalDistrSNP_Y[2][1] - GenoJointDistr[1][6] - GenoJointDistr[1][7];

}