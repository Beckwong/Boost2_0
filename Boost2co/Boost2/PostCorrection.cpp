#include "PostCorrection.h"
#include "stdafx.h"
#include<vector>
#include"main.h"
using namespace std;

double PostCorrection(vector<vector<int>>&GenoJointDistr,int nsamples)
{
	double mu[3][3][2] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	double mu0[3][3][2] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double mu_ij[3][3];
	double mu_jk[3][2];
	double mu_ik[3][2];
	double muError=0;
	double tao=0.0;
	double InteractionMeasure=0; 
	double ptmp1,ptmp2; 

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<2;k++)
			{
				muError += Abs(mu[i][j][k]-mu0[i][j][k]);
			}
		}
	}
	//iteration method to get mu[i][j][k]
	while(muError > 0.001)
	{
		for(int i=0;i<3;i++)
	   {
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					mu0[i][j][k]=mu[i][j][k];
				}
			}
	   }
  
	  
		//mu_ij
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				mu_ij[i][j] = mu[i][j][0]+mu[i][j][1]; 
			}
		}

		//mu_ijk = mu_ijk*n_ij/mu_ij

		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					if(mu_ij[i][j]>0)
						mu[i][j][k] = mu[i][j][k]*(GenoJointDistr[0][3*i+j]+GenoJointDistr[1][3*i+j])/mu_ij[i][j];
					else
						mu[i][j][k]=0;
				}
			}
		}

		//mu_ik
		for(int i=0;i<3;i++)
		{
			for(int k=0;k<2;k++)
			{
				mu_ik[i][k] = mu[i][0][k]+mu[i][1][k]+mu[i][2][k];

			}		
		}
		//mu_ijk = mu_ijk*n_ik/mu_ik

		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					if(mu_ik[i][k]>0)
						mu[i][j][k] = mu[i][j][k]*(GenoJointDistr[k][i*3]+GenoJointDistr[k][i*3+1]+GenoJointDistr[k][i*3+2])/mu_ik[i][k];
					else
						mu[i][j][k]=0;
				}
			}
		}

		//mu_jk
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<2;k++)
			{
				mu_jk[j][k] = mu[0][j][k]+mu[1][j][k]+mu[2][j][k];
			}
		}

		//mu_ijk = mu_ijk*n_jk/mu_jk
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					if(mu_jk[j][k]>0)
						mu[i][j][k] = mu[i][j][k]*(GenoJointDistr[k][0*3+j]+GenoJointDistr[k][1*3+j]+GenoJointDistr[k][2*3+j])/mu_jk[j][k];
				}
			}
		}

		muError = 0.0;
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					muError += Abs(mu[i][j][k]-mu0[i][j][k]);
				}
			}
		}


	}
	

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<2;k++)
			{
				ptmp1 = (double)GenoJointDistr[k][i*3+j]/nsamples;
				if(ptmp1>0)
				{
					InteractionMeasure += ptmp1*log(ptmp1);
				}
				ptmp2 = mu[i][j][k]/nsamples;
				if(ptmp2>0)
				{
					InteractionMeasure += -ptmp1*log(ptmp2);
					tao += ptmp2;
				}
			}
		}
	}

	InteractionMeasure = (InteractionMeasure + log(tao))*nsamples*2;

	return InteractionMeasure;

}