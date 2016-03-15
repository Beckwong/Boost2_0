#include "CalculateMarginalEntropy.h"

using namespace std;

void CalculateMarginalEntropy(vector<vector<genotype*>>&pgeno,vector<vector<int> >&DataSize,vector<int>&nCase_Ctrl,vector<double>&MarginalEntropySNP,vector<double>&MarginalEntropySNP_Y)
{
	int nsnp=0;
	int nsample = DataSize[0][0];
	int count = 0;
	vector<vector<int>>GenoMarginalDistr(3,vector<int>(2));

	for(int i=0;i<DataSize.size();i++)
	{
		nsnp+=DataSize[i][1];
	}

	for(int i=0;i<nsnp;i++)
	{
		for(int j=0;j<3;j++)
		{
			count = 0;
			for(int k=0;k<pgeno[i][j]->genocase.size();k++)
			{
				count+=bitCount(pgeno[i][j]->genocase[k]);
			}
			GenoMarginalDistr[j][0] = count;
			count=0;
			
			for(int k=0;k<pgeno[i][j]->genoctrl.size();k++)
			{
				count+=bitCount(pgeno[i][j]->genoctrl[k]);
			}
			GenoMarginalDistr[j][1] = count;
		}

		for(int j=0;j<3;j++)
		{
			double tmp = GenoMarginalDistr[j][0]+GenoMarginalDistr[j][1];
			double ptmp;
			if(tmp>0)
			{
				ptmp = ptmp/nsample;
				MarginalEntropySNP[i]+= -(ptmp)*log(ptmp);
			}
			if(GenoMarginalDistr[j][0]>0)
			{
				ptmp = (double)GenoMarginalDistr[j][0]/nsample;
				MarginalEntropySNP_Y[i]+= -(ptmp)*log(ptmp);
			}
			if(GenoMarginalDistr[j][1]>0)
			{
				ptmp = (double)GenoMarginalDistr[j][1]/nsample;
				MarginalEntropySNP_Y[i]+= -(ptmp)*log(ptmp);
			}
		}

	}




}