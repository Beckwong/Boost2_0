#include"CalculateMargnialEntropy.h"
#include"main.h"
using namespace std;

void CalculateMarginalEntropy(vector<vector<genotype*>>&pgeno,vector<int>&DataSize,vector<vector<int>>&nCase_Gender,vector<double>&MarginalEntropySNP,vector<double>&MarginalEntropySNP_Y)
{
	int nsnps=0;
	int nsamples=nCase_Gender[0][0]+nCase_Gender[0][1]+nCase_Gender[1][0]+nCase_Gender[1][1];
	int count=0;
	vector<vector<int>>GenoMarginalDistr(3,vector<int>(2));
	for(int i=0;i<DataSize.size();i++)
	{
		nsnps+=DataSize[i];
	}

	for(int i=0;i<nsnps;i++)
	{
		for(int j=0;j<3;j++)
		{
			count=0;
			for(int k=0;k<pgeno[i][j]->geno[0].size();k++)
			{
				count+=bitCount(pgeno[i][j]->geno[0][k]);
			}
			for(int k=0;k<pgeno[i][j]->geno[1].size();k++)
			{
				count+=bitCount(pgeno[i][j]->geno[1][k]);
			}
			GenoMarginalDistr[j][0]=count;
			
			count=0;
			for(int k=0;k<pgeno[i][j]->geno[2].size();k++)
			{
				count+=bitCount(pgeno[i][j]->geno[2][k]);
			}
			for(int k=0;k<pgeno[i][j]->geno[3].size();k++)
			{
				count+=pgeno[i][j]->geno[3][k];
			}
			GenoMarginalDistr[j][1]=count;
		}

		for(int j=0;j<3;j++)
		{
			double tmp = GenoMarginalDistr[j][0]+GenoMarginalDistr[j][1];
			double ptmp;
			if(tmp>0)
			{
				ptmp = tmp/nsamples;
				MarginalEntropySNP[i]+= -(ptmp)*log(ptmp);
			}
			if(GenoMarginalDistr[j][0]>0)
			{
				ptmp = (double)GenoMarginalDistr[j][0]/nsamples;
				MarginalEntropySNP_Y[i]+= -(ptmp)*log(ptmp);
			}
			if(GenoMarginalDistr[j][1]>0)
			{
				ptmp = (double)GenoMarginalDistr[j][1]/nsamples;
				MarginalEntropySNP_Y[i]+= -(ptmp)*log(ptmp);
			}
		}
	}

}