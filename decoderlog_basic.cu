#include <stdio.h>
#include <stdlib.h>
#include "hashdefined.h"

void just_decoder_basic(float snr, int iter_num, float *outbits_x, float *outbits_y1, float *outbits_y2, int *decisionbits, int trellis_term)
{
void result(float LC, float *outbits_x, int *decisionbits, float *Lext12, float *Lext21);
void basic_decoder(float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y,  float LC, int which_decoder, int trellis_term);
float Lext12[DATASIZE], Lext21[DATASIZE];
extern int inv_permutation_bits[DATASIZE];
float LC;
int iter;
int i;
LC = (4.0f/3.0f)*pow(10,snr/10.0);
for(i=0; i<DATASIZE; i++)
{
	Lext21[i] = 0;
	Lext12[i] = 0;
}
for(iter =1; iter<= iter_num; iter++)
{
	basic_decoder(Lext21, Lext12, outbits_x, outbits_y1, LC, 1, trellis_term);
	basic_decoder(Lext12, Lext21, outbits_x, outbits_y2, LC, 2, trellis_term);
	
}
result(LC, outbits_x, decisionbits, Lext12, Lext21);
}


void basic_decoder(float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y,  float LC, int which_decoder, int trellis_term)
{
extern int inv_permutation_bits[DATASIZE];
float maxf(float a, float b);
void guard_alpha(float *start_alpha);
void guard_beta(float *start_beta, int which_decoder, int trellis_term);
int cpu_alpha_state_0[8] 	= 	{0, 3, 4, 7, 1, 2, 5, 6};
int cpu_alpha_state_1[8] 	= 	{1, 2, 5, 6, 0, 3, 4, 7};
int cpu_beta_state_0[8] 	= 	{0, 4, 5, 1, 2, 6, 7, 3};
int cpu_beta_state_1[8] 	= 	{4, 0, 1, 5, 6, 2, 3, 7};
int cpu_alpha_encbit_0[8] 	= 	{-1, 1, 1, -1, -1, 1, 1, -1};
int cpu_alpha_encbit_1[8] 	= 	{1, -1, -1, 1, 1, -1, -1, 1};
int cpu_beta_encbit_0[8] 	= 	{-1, -1, 1, 1, 1, 1, -1, -1};
int cpu_beta_encbit_1[8] 	= 	{1, 1, -1, -1, -1, -1, 1, 1};
float alpha[(DATASIZE+1)*8];
float Lext_input[DATASIZE], channel_x[DATASIZE], channel_y[DATASIZE];
float beta[8], temp_beta[8];
float start_alpha[8], start_beta[8];
float gamma_0, gamma_1;
float alpha_0, alpha_1;
float beta_0, beta_1;
float splus, sminus;
float Lexternal;
int i, j;


for(i=0; i<DATASIZE; i++)
{
	channel_y[i] = outbits_y[i];
}
if(which_decoder == 1)
{
	for(i=0; i<DATASIZE; i++)
	{
		Lext_input[i] = Lext_in[i];
		channel_x[i] = outbits_x[i];
	}
}
else if(which_decoder == 2)
{
	for(i=0; i<DATASIZE; i++)
	{
		Lext_input[i] = Lext_in[inv_permutation_bits[i]];
		channel_x[i] = outbits_x[inv_permutation_bits[i]];
	}
}

guard_alpha(start_alpha);
guard_beta(start_beta, which_decoder, trellis_term);

for(j=0; j<8; j++)
{
	alpha[j] = start_alpha[j];
	beta[j] = start_beta[j];
}

//*******************************************Alpha evaluation**********************************************
for(i=0; i<DATASIZE; i++)
{
	for(j=0; j<8; j++)
	{
		gamma_0 = 0.5*(Lext_input[i] + LC*channel_x[i])*(-1) + 0.5*LC*channel_y[i]*cpu_alpha_encbit_0[j];
		gamma_1 = 0.5*(Lext_input[i] + LC*channel_x[i])*(1) + 0.5*LC*channel_y[i]*cpu_alpha_encbit_1[j];
		alpha_0 = alpha[i*8 + cpu_alpha_state_0[j] ] + gamma_0;
		alpha_1 = alpha[i*8 + cpu_alpha_state_1[j] ] + gamma_1;
		alpha[(i+1)*8 + j] = maxf(alpha_0, alpha_1);
	}
}
//*******************************************Alpha evaluation**********************************************
/************Print alpha test block*******************
FILE *fp1, *fp2, *fp3, *fp4;
fp1 = fopen("alpha_check.dat", "w");
fp2 = fopen("beta_check.dat", "w");
fp3 = fopen("Lext_check.dat", "w");
fp4 = fopen("Full_Lext_check.dat", "w");
for(i=0; i<DATASIZE; i++)
{
	for(j=0; j<8; j++)
	{
		fprintf(fp1, "%f\t", alpha[(i+1)*8+j]);
	}
	fprintf(fp1, "\n");
}
fclose(fp1);
************Print alpha test block*******************/


//*******************************************Beta and Lext evaluation**********************************************
for(i=DATASIZE-1; i>=0; i--)
{
	//***************************************Beta evaluation***********************************************
	if(i == DATASIZE-1)
	{
		for(j=0; j<8; j++)
		{
			beta[j] = EQUAL_GUARD_VALUE;
		}
	}
	else
	{
		for(j=0; j<8; j++)
		{
			gamma_0 = 0.5*(Lext_input[i+1] + LC*channel_x[i+1])*(-1) + 0.5*LC*channel_y[i+1]*cpu_beta_encbit_0[j];
			gamma_1 = 0.5*(Lext_input[i+1] + LC*channel_x[i+1])*(1) + 0.5*LC*channel_y[i+1]*cpu_beta_encbit_1[j];
			beta_0 = beta[cpu_beta_state_0[j] ] + gamma_0;
			beta_1 = beta[cpu_beta_state_1[j] ] + gamma_1;
			temp_beta[j] = maxf(beta_0, beta_1);
		}	
		for(j=0; j<8; j++)
		{
			beta[j] = temp_beta[j];
		}
	}
	/********Print beta block**********
	for(j=0; j<8; j++)
	{
		fprintf(fp2,"%f\t", beta[j]);
	}
	fprintf(fp2, "\n");
	********Print beta block**********/
	//***************************************Beta evaluation***********************************************

	//***************************************Lexternal evaluation******************************************	
	sminus = 0;
	splus =  0;
	for(j=0; j<8; j++)
	{
		gamma_0 = 0.5*(Lext_input[i] + LC*channel_x[i])*(-1) + 0.5*LC*channel_y[i]*cpu_beta_encbit_0[j];
		gamma_1 = 0.5*(Lext_input[i] + LC*channel_x[i])*(1) + 0.5*LC*channel_y[i]*cpu_beta_encbit_1[j];
		sminus = maxf(sminus, alpha[i*8+j] + beta[cpu_beta_state_0[j]] + gamma_0);
		splus = maxf(splus, alpha[i*8+j] + beta[cpu_beta_state_1[j]] + gamma_1); 
	}
	Lexternal = splus - sminus - LC*channel_x[i]- Lext_input[i];
	//Lexternal = Lext_input[i];
	//Lexternal = LC*channel_x[i];
	/*************************************************Mega Lext check*****************************************************************
	fprintf(fp4, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", channel_x[i], Lexternal,sminus, splus,splus - sminus, LC*channel_x[i], Lext_input[i]);
	*************************************************Mega Lext check*****************************************************************/
	//***************************************Lexternal evaluation******************************************
	if(which_decoder == 1)
	{
		Lext_out[i] = Lexternal;
	}
	else if(which_decoder == 2)
	{
		Lext_out[inv_permutation_bits[i]] = Lexternal;
	}
}
//*******************************************Beta and Lext evaluation**********************************************
/************Printing Lexternal*****************
for(i=0; i<DATASIZE; i++)
{
	fprintf(fp3, "%f\n", Lext_out[i]);
}
fclose(fp2);
fclose(fp3);
fclose(fp4);
************Printing Lexternal*****************/
}


void guard_alpha(float *start_alpha)
{
int i;
start_alpha[0] = 0;
for(i=1; i<8; i++)
{	
	start_alpha[i] = MINUS_INFINITY;
}
}

void guard_beta(float *start_beta, int which_decoder, int trellis_term)
{
int i;
if(which_decoder == 1 && trellis_term == 1)
{
	start_beta[0] = 0;
	for(i=1; i<8; i++)
	{
		start_beta[i] = MINUS_INFINITY;
	}
}
else
{
	for(i=0; i<8; i++)
	{
		start_beta[i] = EQUAL_GUARD_VALUE;
	}
}
}

float maxf(float a, float b)
{
float myabs(float a);
float mymax(float a, float b);
return mymax(a,b);
//return mymax(a,b) + log(1+exp((-1)*myabs(a-b)));
}




float myabs(float a)
{
if(a>=0) return a;
else return (-1)*a;
}

float mymax(float a, float b)
{
if(a>=b)
	return a;
else 
	return b;
}


