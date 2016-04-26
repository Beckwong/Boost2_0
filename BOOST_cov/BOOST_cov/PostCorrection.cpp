#include "PostCorrection.h"

using namespace std;

double PostCorrection(vector<vector<int>>&GenoJointDistr,int nsamples,bool flag)
{
	vector<double>mu(36,1);
	vector<double>mu0(36,0);
	vector<double>mu_ij(9,0);
	vector<double>mu_ik(6,0);
	vector<double>mu_jk(6,0);
	vector<double>mu_kt(4,0);
	vector<double>mu_it(6,0);
	vector<double>mu_jt(6,0);
	vector<double>mu_ijk(18,0);


	vector<double>n_ij(9,0);
	vector<double>n_ik(6,0);
	vector<double>n_jk(6,0);
	vector<double>n_it(6,0);
	vector<double>n_jt(6,0);
	vector<double>n_kt(4,0);
	vector<double>n_ijk(18,0);

	double muError=0.0;
	double tao=0;
	double InteractionMeasure = 0.0;
	double ptmp1=0,ptmp2=0;
	
	for(int index=0;index<36;index++)
	{
		muError += abs(mu[index]-mu0[index]);
	}

	while(muError>0.001)
	{
		for(int i=0;i<36;i++)
		{
			mu0[i]=mu[i];
		}

		//step1: mu_ij and n_ij
		for(int i=0;i<9;i++)
		{
			mu_ij[i]=0;
			n_ij[i]=0;
		}
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					for(int t=0;t<2;t++)
					{
						n_ij[i*3+j]+=GenoJointDistr[k*2+t][i*3+j];
						mu_ij[i*3+j]+=mu[i*12+j*4+k*2+t];
					}
				}
			}
		}
		//mu_ijkt = mu_ijkt*n_ij/mu_ij
	    for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					for(int t=0;t<2;t++)
					{
						if(mu_ij[i*3+j]>0)
						{
							mu[i*12+j*4+k*2+t] = mu[i*12+j*4+k*2+t]*n_ij[i*3+j]/mu_ij[i*3+j];
						}
						else
						{
							mu[i*12+j*4+k*2+t]=0;
						}
					}
				}
			}
		}

		//step2: mu_ik and n_ik
		for(int i=0;i<6;i++)
		{
			mu_ik[i]=0;
			n_ik[i]=0;
		}
		for(int i=0;i<3;i++)
		{
			for(int k=0;k<2;k++)
			{
				for(int j=0;j<3;j++)
				{
					for(int t=0;t<2;t++)
					{
						mu_ik[i*2+k]+=mu[i*12+j*4+k*2+t];
						n_ik[i*2+k] +=GenoJointDistr[k*2+t][i*3+j];
					}
				}
			}
		}
		
		//mu_ijkt=mu_ijkt*n_ik/mu_ik

		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					for(int t=0;t<2;t++)
					{
						if(mu_ik[i*2+k]>0)
							mu[i*12+j*4+k*2+t] = mu[i*12+j*4+k*2+t]*n_ik[i*2+k]/mu_ik[i*2+k];
						else
							mu[i*12+j*4+k*2+t] =0;
					}
				}
			}
		}

		//step3:mu_jk and n_jk
		for(int i=0;i<6;i++)
		{
			mu_jk[i]=0;
			n_jk[i]=0;
		}
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<2;k++)
			{
				for(int i=0;i<3;i++)
				{
					for(int t=0;t<2;t++)
					{
						mu_jk[j*2+k]+= mu[i*12+j*4+k*2+t];
						n_jk[j*2+k] += GenoJointDistr[k*2+t][i*3+j];
					}
				}
			}
		}

		//mu_ijkt = mu_ijkt*n_jk/mu_jk
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					for(int t=0;t<2;t++)
					{
						if(mu_jk[j*2+k]>0)
							mu[i*12+j*4+k*2+t] = mu[i*12+j*4+k*2+t]*n_jk[j*2+k]/mu_jk[j*2+k];
						else
							mu[i*12+j*4+k*2+t]=0;
					}
				}
			}
		}

		//step4: mu_kt and n_kt
		for(int i=0;i<4;i++)
		{
			mu_kt[i] = 0;
			n_kt[i] = 0;
		}
		for(int k=0;k<2;k++)
		{
			for(int t=0;t<2;t++)
			{
				for(int i=0;i<3;i++)
				{
					for(int j=0;j<3;j++)
					{
						mu_kt[k*2+t] += mu[i*12+j*4+k*2+t];
						n_kt[k*2+t] += GenoJointDistr[k*2+t][i*3+j];
					}
				}
			}
		}

		//mu_ijkt = mu_ijkt*n_kt/mu_kt
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					for(int t=0;t<2;t++)
					{
						if(mu_kt[k*2+t]>0)
							mu[i*12+j*4+k*2+t] = mu[i*12+j*4+k*2+t]*n_kt[k*2+t]/mu_kt[k*2+t];
						else
							mu[i*12+j*4+k*2+t]=0;
					}
				}
			}
		}

		//step 5 u_it and n_it
		for(int index=0;index<6;index++)
		{
			mu_it[index]=0;
			n_it[index]=0;
		}
		for(int i=0;i<3;i++)
		{
			for(int t=0;t<2;t++)
			{
				for(int j=0;j<3;j++)
				{
					for(int k=0;k<2;k++)
					{
						mu_it[i*2+t] += mu[i*12+j*4+k*2+t];
						n_it[i*2+t] +=  GenoJointDistr[k*2+t][i*3+j];
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
					for(int t=0;t<2;t++)
					{
						if(mu_it[i*2+t]>0)
							mu[i*12+j*4+k*2+t] = mu[i*12+j*4+k*2+t]*n_it[i*2+t]/mu_it[i*2+t];
						else
							mu[i*12+j*4+k*2+t] = 0;
					}
				}
			}
		}

		//step 6 u_jt and n_jt
		for(int index=0;index<6;index++)
		{
			mu_jt[index]=0;
			n_jt[index]=0;
		}
		for(int j=0;j<3;j++)
		{
			for(int t=0;t<2;t++)
			{
				for(int i=0;i<3;i++)
				{
					for(int k=0;k<2;k++)
					{
						mu_jt[j*2+t] += mu[i*12+j*4+k*2+t];
						n_jt[j*2+t] += GenoJointDistr[k*2+t][i*3+j];
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
					for(int t=0;t<2;t++)
					{
						if(mu_jt[j*2+t]>0)
							mu[i*12+j*4+k*2+t] = mu[i*12+j*4+k*2+t]*n_jt[j*2+t]/mu_jt[j*2+t];
						else
							mu[i*12+j*4+k*2+t] = 0;
					}
				}
			}
		}


		//step7: for association detection
		//mu_ijkt = mu_ijkt*n_ijk/mu_ijk
		if(flag)
		{
			for(int i=0;i<18;i++)
			{
				mu_ijk[i]=0;
				n_ijk[i] =0;
			}

			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					for(int k=0;k<2;k++)
					{
						for(int t=0;t<2;t++)
						{
							mu_ijk[i*6+j*2+k]+=mu[i*12+j*4+k*2+t];
							n_ijk[i*6+j*2+k]+=GenoJointDistr[k*2+t][i*3+j];
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
						for(int t=0;t<2;t++)
						{
							if(mu_ijk[i*6+j*2+k]>0)
								mu[i*12+j*4+k*2+t] = mu[i*12+j*4+k*2+t]*n_ijk[i*6+j*2+k]/mu_ijk[i*6+j*2+k];
							else
								mu[i*12+j*4+k*2+t]=0;
						}
					}
				}
			}
		}

		muError=0.0;
		for(int index=0;index<36;index++)
		{
			muError += abs(mu[index]-mu0[index]);
		}



	}

	double Likelihood=0;
	double tmp=0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<2;k++)
			{
				for(int t=0;t<2;t++)
				{
					tmp=GenoJointDistr[k*2+t][i*3+j];
					if(mu[i*12+j*4+k*2+t]>0)
					{
						Likelihood+=tmp*log(mu[i*12+j*4+k*2+t]);
					}
				}
			}
		}
	}

	return Likelihood;


}