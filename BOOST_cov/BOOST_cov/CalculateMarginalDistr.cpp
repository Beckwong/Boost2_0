#include"CalculateMarginalDistr.h"

using namespace std;

void CalculateMarginalDistr(vector<vector<genotype*>>&pgeno,int nsnps,vector<MarginalDistr*>&pMarginalDistr)
{
	int count=0;

	for(int i=0;i<nsnps;i++)
	{
		for(int j=0;j<3;j++)
		{
			
			for(int m=0;m<2;m++)
			{
				for(int n=0;n<2;n++)
				{
					count=0;
					for(int index=0;index<pgeno[i][j]->geno[m*2+n].size();index++)
						count+=bitCount(pgeno[i][j]->geno[m*2+n][index]);
					pMarginalDistr[i]->MarginalDistrSNP_Y_G[m*2+n][j]=count;
				}
			}
			for(int m=0;m<2;m++)
			{
				for(int n=0;n<2;n++)
				{
					pMarginalDistr[i]->MarginalDistrSNP[j] += pMarginalDistr[i]->MarginalDistrSNP_Y_G[m*2+n][j];
				}
			}
			
		}
	}
}