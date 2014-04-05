#include <stdio.h>
#include <stdlib.h>
#include "hashdefined.h"

void just_decoder_parallel_close(float snr, int iter_num, int blocks, int guard_size, float *outbits_x, float *outbits_y1, float *outbits_y2, int *decisionbits, int trellis_term, int guarding_type)
{
void result(float LC, float *outbits_x, int *decisionbits, float *Lext12, float *Lext21);
void parallel_decoder_close(int blocks, int guard_size, float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y,  float LC, int which_decoder);
float Lext12[DATASIZE], Lext21[DATASIZE];
extern int inv_permutation_bits[DATASIZE];
float LC;
int iter;
int i;
LC = (4.0/3.0)*pow(10,snr/10.0);
for(i=0; i<DATASIZE; i++)
{
	Lext21[i] = 0;
	Lext12[i] = 0;
}
for(iter =1; iter<= iter_num; iter++)
{
	parallel_decoder_close(blocks, guard_size, Lext21, Lext12, outbits_x, outbits_y1, LC, 1);
	parallel_decoder_close(blocks, guard_size, Lext12, Lext21, outbits_x, outbits_y2, LC, 2);	
}
/********Lext12 and Lext21 check block**************
FILE *fp;
int pi;
fp = fopen("Lext_host.dat", "w");
for(pi=0; pi<DATASIZE; pi++)
{
	fprintf(fp, "%f\t%f\t%f\t%f\n", LC*outbits_x[pi], Lext12[pi], Lext21[pi], LC*outbits_x[pi] + Lext12[pi] + Lext21[pi]);
}
fclose(fp);
********Lext12 and Lext21 check block**************/
result(LC, outbits_x, decisionbits, Lext12, Lext21);
}


void parallel_decoder_close(int blocks, int guard_size, float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y,  float LC, int which_decoder)
{
void figment_decoder_close(int figment, int blocks, int guard_size, float *Lext_in, float *Lext_out, float *channel_x, float *channel_y,  float LC, int which_decoder);
int i;
//******************************************Using the figment decoder to do parallel decoding******************************************
for(i=0; i<blocks; i++)
{
	figment_decoder_close(i, blocks, guard_size, Lext_in, Lext_out, outbits_x, outbits_y, LC, which_decoder);
}
//******************************************Using the figment decoder to do parallel decoding******************************************
}

void figment_decoder_close(int figment, int blocks, int guard_size, float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y,  float LC, int which_decoder)
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
int block_size;
float *alpha, *Lext, *channel_x, *channel_y;
float *lambda_0, *lambda_1;
block_size = DATASIZE/blocks;
alpha = (float*)malloc((block_size+1)*8*sizeof(float));
Lext = (float*)malloc((block_size+1)*sizeof(float));
channel_x = (float*)malloc((block_size+1)*sizeof(float));
channel_y = (float*)malloc((block_size+1)*sizeof(float));
lambda_0 = (float*)malloc(8*8*sizeof(float));
lambda_1 = (float*)malloc(8*8*sizeof(float));

float beta[8], temp_beta[8];
float start_alpha[8], start_beta[8];
float gamma_0, gamma_1;
float alpha_0, alpha_1;
float beta_0, beta_1;
float Lext_0, Lext_1;
float Lexternal;
int i, j, k, l, tid, index;
int block_end;
block_end = (figment == blocks-1) ? block_size : block_size+1;
for(i=0; i<block_end; i++)
{
	channel_y[i] = outbits_y[figment*block_size+i];
}
if(which_decoder == 1)
{
	for(i=0; i<block_end; i++)
	{
		Lext[i] = Lext_in[figment*block_size+i];
		channel_x[i] = outbits_x[figment*block_size+i];
	}
}
else if(which_decoder == 2)
{
	for(i=0; i<block_end; i++)
	{
		Lext[i] = Lext_in[inv_permutation_bits[figment*block_size+i]];
		channel_x[i] = outbits_x[inv_permutation_bits[figment*block_size+i]];
	}
}

guard_alpha(start_alpha);
guard_beta(start_beta, 0, 1);

//*******************************************Alpha evaluation**********************************************
for(tid=0; tid<8; tid++)
{
	alpha[tid]=start_alpha[tid];
}
for(i=0; i<block_size; i++)
{
	for(tid=0; tid<8; tid++)
	{
		gamma_0 = (0.5f)*(Lext[i] + LC*channel_x[i])*(-1) + (0.5f)*LC*channel_y[i]*cpu_alpha_encbit_0[tid];
		gamma_1 = (0.5f)*(Lext[i] + LC*channel_x[i])*(1) + (0.5f)*LC*channel_y[i]*cpu_alpha_encbit_1[tid];
		alpha_0 = alpha[i*8 + cpu_alpha_state_0[tid] ] + gamma_0;
		alpha_1 = alpha[i*8 + cpu_alpha_state_1[tid] ] + gamma_1;
		alpha[(i+1)*8 + tid] = maxf(alpha_0, alpha_1);
	}
}
//*******************************************Alpha evaluation**********************************************
/************Print alpha test block*******************
int pj;
FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6;
fp1 = fopen("alpha_check_2.dat", "w");
fp2 = fopen("beta_check_2.dat", "w");
fp3 = fopen("Lext_check_2.dat", "w");
fp4 = fopen("Full_Lext_check_2.dat", "w");
fp5 = fopen("lambda_0_check.dat", "w");
fp6 = fopen("lambda_1_check.dat",  "w");
for(i=0; i<block_size; i++)
{
	for(pj=0; pj<8; pj++)
	{
		fprintf(fp1, "%f\t", alpha[(i+1)*8+pj]);
	}
	fprintf(fp1, "\n");
}
fclose(fp1);
************Print alpha test block*******************/
//******************************Beta and Lext evaluation*******************************
for(tid=0; tid<8; tid++)
{
	beta[tid] = start_beta[tid];
}
for(k=block_size/8; k>0; k--)
{
	for(j=0; j<8; j++)
	{
		i = k*8-j-1;
		if((figment == blocks-1) && (i == block_size-1))
		{
			for(tid=0; tid<8; tid++)
			{
				temp_beta[tid] = EQUAL_GUARD_VALUE;
				beta[tid] = EQUAL_GUARD_VALUE;
			}
		}
		else
		{
			for(tid=0; tid<8; tid++)
			{
				//**********************************Beta evaluation****************************
				gamma_0 = (0.5f)*(LC*channel_x[i+1]+Lext[i+1])*(-1) + (0.5f)*LC*channel_y[i+1]*cpu_beta_encbit_0[tid];
				gamma_1 = (0.5f)*(LC*channel_x[i+1]+Lext[i+1])*(1) + (0.5f)*LC*channel_y[i+1]*cpu_beta_encbit_1[tid];
				beta_0 = beta[cpu_beta_state_0[tid]] + gamma_0;
				beta_1 = beta[cpu_beta_state_1[tid]] + gamma_1;
				temp_beta[tid] = maxf(beta_0, beta_1);
			}
			for(tid=0; tid<8; tid++)
			{
				beta[tid] = temp_beta[tid];
			}
		}
		for(tid=0; tid<8; tid++)
		{
			//**********************************Lambda evaluation**************************
			gamma_0 = (0.5f)*(LC*channel_x[i]+Lext[i])*(-1) + (0.5f)*LC*channel_y[i]*cpu_beta_encbit_0[tid];
			gamma_1 = (0.5f)*(LC*channel_x[i]+Lext[i])*(1) + (0.5f)*LC*channel_y[i]*cpu_beta_encbit_1[tid];
			lambda_0[j*8+tid] = alpha[i*8+tid] + temp_beta[cpu_beta_state_0[tid]] + gamma_0;
			lambda_1[j*8+tid] = alpha[i*8+tid] + temp_beta[cpu_beta_state_1[tid]] + gamma_1;
			//**********************************Lambda evaluation**************************/
		}
		/**************Print lambda block**************
		for(pj=0; pj<8; pj++)
		{
			fprintf(fp5, "%f\t", lambda_0[j*8+pj]);
			fprintf(fp6, "%f\t", lambda_1[j*8+pj]);
		}
		fprintf(fp5, "\n");
		fprintf(fp6, "\n");
		**************Print lambda block**************/	
		/********Print beta block**********
		for(pj=0; pj<8; pj++)
		{
			fprintf(fp2,"%f\t", beta[pj]);
		}
		fprintf(fp2, "\n");
		********Print beta block************/
	}
	for(tid=0; tid<8; tid++)
	{
		//****************************Lexternal evaluation****************************
		Lext_0 = lambda_0[tid*8];
		Lext_1 = lambda_1[tid*8];
		for(l=1; l<=7; l++)
		{
				Lext_0 = maxf(Lext_0, lambda_0[tid*8 + l]);
				Lext_1 = maxf(Lext_1, lambda_1[tid*8 + l]);
		}
		Lexternal = (Lext_1 - Lext_0) - LC*channel_x[k*8 - tid - 1] - Lext[k*8 - tid - 1];
		/****************Full Lext print check**********************
		fprintf(fp4, "%f\t%f\t%f\t%f\t%f\n", channel_x[i+7-tid], Lexternal,Lext_0, Lext_1, Lext[i+7-tid]);
		****************Full Lext print check**********************/
		index = figment*block_size + k*8 - tid - 1;
		if(which_decoder == 1)
		{
			Lext_out[index] = Lexternal;
		}
		else if(which_decoder == 2)
		{
			Lext_out[inv_permutation_bits[index]] = Lexternal;
		}
		//****************************Lexternal evaluation****************************/	
	}	
}
/************Printing Lexternal*****************
for(i=0; i<DATASIZE; i++)
{
	fprintf(fp3, "%f\n", Lext_out[i]);
}
fclose(fp2);
fclose(fp3);
fclose(fp4);
fclose(fp5);
fclose(fp6);
************Printing Lexternal*****************/
free(alpha);
free(Lext);
free(channel_x);
free(channel_y);
free(lambda_0);
free(lambda_1);
}
