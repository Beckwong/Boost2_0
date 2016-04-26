#include "CalculateGenoJointDistr.h"
#include"main.h"
using namespace std;

void CalculateGenoJointDistr(vector<vector<genotype*>>&pgeno,vector<vector<int>>&GenoJointDistr,int snp1,int snp2,vector<MarginalDistr*>&pMarginalDistr)
{
	register int count;

	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			for(int m=0;m<2;m++)
			{
				for(int n=0;n<2;n++)
				{
					count=0;
					for(int index=0;index<pgeno[snp1][i]->geno[m*2+n].size();index++)
					{
						count+=popcount(pgeno[snp1][i]->geno[m*2+n][index] & pgeno[snp2][j]->geno[m*2+n][index]);
					}
					GenoJointDistr[m*2+n][i*3+j] = count;
				}
			}
		}
	}
	//control female
	GenoJointDistr[0][2] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[0][0]-GenoJointDistr[0][0]-GenoJointDistr[0][1];
	GenoJointDistr[0][5] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[0][1]-GenoJointDistr[0][3]-GenoJointDistr[0][4];
	GenoJointDistr[0][6] = pMarginalDistr[snp2]->MarginalDistrSNP_Y_G[0][0]-GenoJointDistr[0][0]-GenoJointDistr[0][3];
	GenoJointDistr[0][7] = pMarginalDistr[snp2]->MarginalDistrSNP_Y_G[0][1]-GenoJointDistr[0][1]-GenoJointDistr[0][4];
	GenoJointDistr[0][8] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[0][2]-GenoJointDistr[0][6]-GenoJointDistr[0][7];

	//control male
	GenoJointDistr[1][2] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[1][0]-GenoJointDistr[1][0]-GenoJointDistr[1][1];
	GenoJointDistr[1][5] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[1][1]-GenoJointDistr[1][3]-GenoJointDistr[1][4];
	GenoJointDistr[1][6] = pMarginalDistr[snp2]->MarginalDistrSNP_Y_G[1][0]-GenoJointDistr[1][0]-GenoJointDistr[1][3];
	GenoJointDistr[1][7] = pMarginalDistr[snp2]->MarginalDistrSNP_Y_G[1][1]-GenoJointDistr[1][1]-GenoJointDistr[1][4];
	GenoJointDistr[1][8] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[1][2]-GenoJointDistr[1][6]-GenoJointDistr[1][7];

	//case female
	GenoJointDistr[2][2] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[2][0]-GenoJointDistr[2][0]-GenoJointDistr[2][1];
	GenoJointDistr[2][5] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[2][1]-GenoJointDistr[2][3]-GenoJointDistr[2][4];
	GenoJointDistr[2][6] = pMarginalDistr[snp2]->MarginalDistrSNP_Y_G[2][0]-GenoJointDistr[2][0]-GenoJointDistr[2][3];
	GenoJointDistr[2][7] = pMarginalDistr[snp2]->MarginalDistrSNP_Y_G[2][1]-GenoJointDistr[2][1]-GenoJointDistr[2][4];
	GenoJointDistr[2][8] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[2][2]-GenoJointDistr[2][6]-GenoJointDistr[2][7];

	//case male
	GenoJointDistr[3][2] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[3][0]-GenoJointDistr[3][0]-GenoJointDistr[3][1];
	GenoJointDistr[3][5] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[3][1]-GenoJointDistr[3][3]-GenoJointDistr[3][4];
	GenoJointDistr[3][6] = pMarginalDistr[snp2]->MarginalDistrSNP_Y_G[3][0]-GenoJointDistr[3][0]-GenoJointDistr[3][3];
	GenoJointDistr[3][7] = pMarginalDistr[snp2]->MarginalDistrSNP_Y_G[3][1]-GenoJointDistr[3][1]-GenoJointDistr[3][4];
	GenoJointDistr[3][8] = pMarginalDistr[snp1]->MarginalDistrSNP_Y_G[3][2]-GenoJointDistr[3][6]-GenoJointDistr[3][7];
}