#include <stdio.h>
#include <stdlib.h>
#include "hashdefined.h"

void just_decoder_parallel(float snr, int iter_num, int blocks, int guard_size, float *outbits_x, float *outbits_y1, float *outbits_y2, int *decisionbits, int trellis_term, int guarding_type)
{
void result(float LC, float *outbits_x, int *decisionbits, float *Lext12, float *Lext21);
void parallel_decoder(int blocks, int guard_size, float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y, float *read_ag, float *read_bg , float *write_ag, float *write_bg, float LC, int which_decoder, int trellis_term, int guarding_type);
float Lext12[DATASIZE], Lext21[DATASIZE];
float *read_alpha_guard_1, *read_beta_guard_1, *read_alpha_guard_2, *read_beta_guard_2;
float *write_alpha_guard_1, *write_beta_guard_1, *write_alpha_guard_2, *write_beta_guard_2;
extern int inv_permutation_bits[DATASIZE];
float LC;
int iter;
int i;

read_alpha_guard_1 = (float*)malloc(blocks*8*sizeof(float));
read_beta_guard_1 = (float*)malloc(blocks*8*sizeof(float));
read_alpha_guard_2 = (float*)malloc(blocks*8*sizeof(float));
read_beta_guard_2 = (float*)malloc(blocks*8*sizeof(float));
write_alpha_guard_1 = (float*)malloc(blocks*8*sizeof(float));
write_beta_guard_1 = (float*)malloc(blocks*8*sizeof(float));
write_alpha_guard_2 = (float*)malloc(blocks*8*sizeof(float));
write_beta_guard_2 = (float*)malloc(blocks*8*sizeof(float));
for(i=0; i<blocks*8; i++)
{
	read_alpha_guard_1[i] = EQUAL_GUARD_VALUE;
	read_alpha_guard_2[i] = EQUAL_GUARD_VALUE;
	read_beta_guard_1[i] = EQUAL_GUARD_VALUE;
	read_beta_guard_2[i] = EQUAL_GUARD_VALUE;
}

LC = (4.0/3.0)*pow(10,snr/10.0);
for(i=0; i<DATASIZE; i++)
{
	Lext21[i] = 0;
	Lext12[i] = 0;
}
for(iter =1; iter<= iter_num; iter++)
{
	parallel_decoder(blocks, guard_size, Lext21, Lext12, outbits_x, outbits_y1, read_alpha_guard_1, read_beta_guard_1, write_alpha_guard_1, write_beta_guard_1, LC, 1, trellis_term, guarding_type);
	parallel_decoder(blocks, guard_size, Lext12, Lext21, outbits_x, outbits_y2, read_alpha_guard_2, read_beta_guard_2, write_alpha_guard_2, write_beta_guard_2, LC, 2, trellis_term, guarding_type);
}
result(LC, outbits_x, decisionbits, Lext12, Lext21);
}


void parallel_decoder(int blocks, int guard_size, float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y, float *read_ag, float *read_bg , float *write_ag, float *write_bg, float LC, int which_decoder, int trellis_term, int guarding_type)
{
void figment_decoder(int figment, int blocks, int guard_size, float *Lext_in, float *Lext_out, float *channel_x, float *channel_y, float *read_ag, float *read_bg, float *write_ag, float *write_bg, float LC, int which_decoder, int trellis_term, int guarding_type);
int i;
//******************************************Using the figment decoder to do parallel decoding******************************************
for(i=0; i<blocks; i++)
{
	figment_decoder(i, blocks, guard_size, Lext_in, Lext_out, outbits_x, outbits_y, read_ag, read_bg, write_ag, write_bg, LC, which_decoder, trellis_term, guarding_type);
}
//Synchronizing read, write : alphas and betas
for(i=0; i<blocks*8; i++)
{
	read_ag[i] = write_ag[i];
	read_bg[i] = write_bg[i];
}
//Synchronizing read, write : alphas and betas
//******************************************Using the figment decoder to do parallel decoding******************************************
}

void figment_decoder(int figment, int blocks, int guard_size, float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y, float *read_ag, float *read_bg, float *write_ag, float *write_bg, float LC, int which_decoder, int trellis_term, int guarding_type)
{
extern int inv_permutation_bits[DATASIZE];
float maxf(float a, float b);
void p_guard_alpha(float *start_alpha, float *guard_alpha, int figment, int guard_size, float *channel_x_before, float *channel_y_before, float *Lext_before, float LC, int guarding_type);
void p_guard_beta(float *start_beta, float *guard_beta, int figment, int guard_size, float *channel_x_after, float *channel_y_after, float *Lext_after, float LC, int blocks, int which_decoder, int trellis_term, int guarding_type);
int cpu_alpha_state_0[8] 	= 	{0, 3, 4, 7, 1, 2, 5, 6};
int cpu_alpha_state_1[8] 	= 	{1, 2, 5, 6, 0, 3, 4, 7};
int cpu_beta_state_0[8] 	= 	{0, 4, 5, 1, 2, 6, 7, 3};
int cpu_beta_state_1[8] 	= 	{4, 0, 1, 5, 6, 2, 3, 7};
int cpu_alpha_encbit_0[8] 	= 	{-1, 1, 1, -1, -1, 1, 1, -1};
int cpu_alpha_encbit_1[8] 	= 	{1, -1, -1, 1, 1, -1, -1, 1};
int cpu_beta_encbit_0[8] 	= 	{-1, -1, 1, 1, 1, 1, -1, -1};
int cpu_beta_encbit_1[8] 	= 	{1, 1, -1, -1, -1, -1, 1, 1};
int block_size;
float *alpha, *Lext_input, *channel_x, *channel_y;
float *Lext_before, *Lext_after, *channel_x_before, *channel_x_after, *channel_y_before, *channel_y_after;
block_size = DATASIZE/blocks;
alpha = (float*)malloc((block_size+1)*8*sizeof(float));
Lext_input = (float*)malloc((block_size+1)*sizeof(float));
channel_x = (float*)malloc((block_size+1)*sizeof(float));
channel_y = (float*)malloc((block_size+1)*sizeof(float));

channel_x_before = (float*)malloc(guard_size*sizeof(float));
channel_y_before = (float*)malloc(guard_size*sizeof(float));
Lext_before = (float*)malloc(guard_size*sizeof(float));
channel_x_after = (float*)malloc(guard_size*sizeof(float));
channel_y_after = (float*)malloc(guard_size*sizeof(float));
Lext_after = (float*)malloc(guard_size*sizeof(float));


float beta[8], temp_beta[8];
float start_alpha[8], start_beta[8];
float gamma_0, gamma_1;
float alpha_0, alpha_1;
float beta_0, beta_1;
float splus, sminus, lambda_0, lambda_1;
float Lexternal;
int i, j;
int block_end;
//************Fetching channel and Lext values for decoding*****************
block_end = (figment == blocks-1) ? block_size : block_size+1;
for(i=0; i<block_end; i++)
{
	channel_y[i] = outbits_y[figment*block_size+i];
}
if(which_decoder == 1)
{
	for(i=0; i<block_end; i++)
	{
		Lext_input[i] = Lext_in[figment*block_size+i];
		channel_x[i] = outbits_x[figment*block_size+i];
	}
}
else if(which_decoder == 2)
{
	for(i=0; i<block_end; i++)
	{
		Lext_input[i] = Lext_in[inv_permutation_bits[figment*block_size+i]];
		channel_x[i] = outbits_x[inv_permutation_bits[figment*block_size+i]];
	}
}
//************Fetching channel and Lext values for decoding*****************

//**********************Fetching values for performing guarding**********************
//**********************Before the window fetching**********************
if(figment != 0)
{
	for(i=0; i<guard_size; i++)
	{
		channel_y_before[i] = outbits_y[figment*block_size - guard_size + i];
	}
	if(which_decoder == 1)
	{
		for(i=0; i<guard_size; i++)
		{
			channel_x_before[i] = outbits_x[figment*block_size - guard_size + i];
			Lext_before[i] = Lext_in[figment*block_size - guard_size + i];
		}
	}
	else if(which_decoder == 2)
	{
		for(i=0; i<guard_size; i++)
		{
			channel_x_before[i] = outbits_x[inv_permutation_bits[figment*block_size - guard_size + i]];
			Lext_before[i] = Lext_in[inv_permutation_bits[figment*block_size - guard_size + i]];
		}
	}
}
//**********************Before the window fetching**********************

//**********************After the window fetching**********************
if(figment != blocks-1)
{
	for(i=0; i<guard_size; i++)
	{
		channel_y_after[i] = outbits_y[(figment+1)*block_size +1 + i];
	}
	if(which_decoder == 1)
	{
		for(i=0; i<guard_size; i++)
		{
			channel_x_after[i] = outbits_x[(figment+1)*block_size +1 + i];
			Lext_after[i] = Lext_in[(figment+1)*block_size +1 + i];			
		}
	}
	else if(which_decoder == 2)
	{
		for(i=0; i<guard_size; i++)
		{
			channel_x_after[i] = outbits_x[inv_permutation_bits[(figment+1)*block_size +1 + i]];
			Lext_after[i] = Lext_in[inv_permutation_bits[(figment+1)*block_size +1 + i]];
			
		}
	}
}
//**********************After the window fetching**********************
//**********************Fetching values for performing guarding**********************


p_guard_alpha(start_alpha, read_ag, figment, guard_size, channel_x_before, channel_y_before, Lext_before, LC, guarding_type);
p_guard_beta(start_beta, read_bg, figment, guard_size, channel_x_after, channel_y_after, Lext_after, LC, blocks, which_decoder, trellis_term, guarding_type);

for(j=0; j<8; j++)
{
	alpha[j] = start_alpha[j];
	beta[j] = start_beta[j];
}


//*******************************************Alpha evaluation**********************************************
for(i=0; i<block_size; i++)
{
	for(j=0; j<8; j++)
	{
		gamma_0 = (0.5f)*(Lext_input[i] + LC*channel_x[i])*(-1) + (0.5f)*LC*channel_y[i]*cpu_alpha_encbit_0[j];
		gamma_1 = (0.5f)*(Lext_input[i] + LC*channel_x[i])*(1) + (0.5f)*LC*channel_y[i]*cpu_alpha_encbit_1[j];
		alpha_0 = alpha[i*8 + cpu_alpha_state_0[j] ] + gamma_0;
		alpha_1 = alpha[i*8 + cpu_alpha_state_1[j] ] + gamma_1;
		alpha[(i+1)*8 + j] = maxf(alpha_0, alpha_1);
	}
}
//***************************Writing end alpha for use in the next iteration***************************
if(guarding_type == 1)
{
	for(j=0; j<8; j++)
	{
		write_ag[figment*8+j] = alpha[block_size*8+j];
	}
}
else if(guarding_type == 3)
{
	for(j=0; j<8; j++)
	{
		write_ag[figment*8+j] = alpha[(block_size-guard_size)*8+j];
	}
}
//***************************Writing end alpha for use in the next iteration***************************

//*******************************************Alpha evaluation**********************************************
/************Print alpha test block*******************
FILE *fp1, *fp2, *fp3, *fp4;
fp1 = fopen("alpha_check_1.dat", "w");
fp2 = fopen("beta_check_1.dat", "w");
fp3 = fopen("Lext_check_1.dat", "w");
fp4 = fopen("Full_Lext_check_1.dat", "w");
for(i=0; i<block_size; i++)
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
for(i=block_size-1; i>=0; i--)
{
	//***************************************Beta evaluation***********************************************
	if((figment == blocks-1) && (i == block_size-1))
	{
		for(j=0; j<8; j++)
		{
			beta[j] = EQUAL_GUARD_VALUE;
			temp_beta[j] = EQUAL_GUARD_VALUE;
		}
	}
	else
	{
		for(j=0; j<8; j++)
		{
			gamma_0 = (0.5f)*(Lext_input[i+1] + LC*channel_x[i+1])*(-1) + (0.5f)*LC*channel_y[i+1]*cpu_beta_encbit_0[j];
			gamma_1 = (0.5f)*(Lext_input[i+1] + LC*channel_x[i+1])*(1) + (0.5f)*LC*channel_y[i+1]*cpu_beta_encbit_1[j];
			beta_0 = beta[cpu_beta_state_0[j] ] + gamma_0;
			beta_1 = beta[cpu_beta_state_1[j] ] + gamma_1;
			temp_beta[j] = maxf(beta_0, beta_1);
		}	
		for(j=0; j<8; j++)
		{
			beta[j] = temp_beta[j];
		}
	}
	//***********************Writing beta for initialisation in the next iteration*************************
	if(guarding_type == 1)
	{
		if(i == 0)
		{
			for(j=0; j<8; j++)
			{
				write_bg[figment*8+j] = beta[j];
			}
		}
	}
	else if(guarding_type == 3)
	{
		if(i == guard_size)
		{
			for(j=0; j<8; j++)
			{
				write_bg[figment*8+j] = beta[j];
			}
		}
	}
	//***********************Writing beta for initialisation in the next iteration*************************
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
		gamma_0 = (0.5f)*(Lext_input[i] + LC*channel_x[i])*(-1) + (0.5f)*LC*channel_y[i]*cpu_beta_encbit_0[j];
		gamma_1 = (0.5f)*(Lext_input[i] + LC*channel_x[i])*(1) + (0.5f)*LC*channel_y[i]*cpu_beta_encbit_1[j];
		lambda_0 = alpha[i*8+j] + beta[cpu_beta_state_0[j]] + gamma_0;
		lambda_1 = alpha[i*8+j] + beta[cpu_beta_state_1[j]] + gamma_1;
		sminus = maxf(sminus, lambda_0);
		splus = maxf(splus, lambda_1);
	}
	Lexternal = splus - sminus - LC*channel_x[i]- Lext_input[i];
	/*************************************************Mega Lext check*****************************************************************
	fprintf(fp4, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", channel_x[i], Lexternal,sminus, splus,splus - sminus, LC*channel_x[i], Lext_input[i]);
	*************************************************Mega Lext check*****************************************************************/
	//***************************************Lexternal evaluation******************************************
	if(which_decoder == 1)
	{
		Lext_out[figment*block_size + i] = Lexternal;
	}
	else if(which_decoder == 2)
	{
		Lext_out[inv_permutation_bits[figment*block_size + i]] = Lexternal;
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
free(alpha);
free(channel_x);
free(channel_y);
free(Lext_input);
}



void p_guard_alpha(float *start_alpha, float *guard_alpha, int figment, int guard_size, float *channel_x_before, float *channel_y_before, float *Lext_before, float LC, int guarding_type)
{
float maxf(float a, float b);
int i,j,k;
float gamma_0, gamma_1, alpha_0, alpha_1;
float temp_alpha[8];
int cpu_alpha_state_0[8] 	= 	{0, 3, 4, 7, 1, 2, 5, 6};
int cpu_alpha_state_1[8] 	= 	{1, 2, 5, 6, 0, 3, 4, 7};
int cpu_alpha_encbit_0[8] 	= 	{-1, 1, 1, -1, -1, 1, 1, -1};
int cpu_alpha_encbit_1[8] 	= 	{1, -1, -1, 1, 1, -1, -1, 1};
if(figment == 0)
{
	start_alpha[0] = 0;
	for(i=1; i<8; i++)
	{
		start_alpha[i] = MINUS_INFINITY;
	}
		
}
else
{
	if(guarding_type == 1)
	{
		for(i=0; i<8; i++)
		{
			start_alpha[i] = guard_alpha[(figment-1)*8 + i];
		}
	
	}
	else if(guarding_type == 2 || guarding_type == 3)
	{
		if(guarding_type == 2)
		{
			for(j=0; j<8; j++)
			{
				start_alpha[j] = EQUAL_GUARD_VALUE;
			}
		}
		else if(guarding_type == 3)
		{
			for(j=0; j<8; j++)
			{
				start_alpha[j] = guard_alpha[(figment-1)*8 + j];
			}
		}
		for(k=0; k<guard_size; k++)
		{
			for(j=0; j<8; j++)
			{
				gamma_0 = (0.5f)*(Lext_before[k] + LC*channel_x_before[k])*(-1) + (0.5f)*LC*channel_y_before[k]*cpu_alpha_encbit_0[j];
				gamma_1 = (0.5f)*(Lext_before[k] + LC*channel_x_before[k])*(1) + (0.5f)*LC*channel_y_before[k]*cpu_alpha_encbit_1[j];
				alpha_0 = start_alpha[cpu_alpha_state_0[j] ] + gamma_0;
				alpha_1 = start_alpha[cpu_alpha_state_1[j] ] + gamma_1;
				temp_alpha[j] = maxf(alpha_0, alpha_1);
			}
			for(j=0; j<8; j++)
			{
				start_alpha[j] = temp_alpha[j];
			}
		}
	}
}
}

void p_guard_beta(float *start_beta, float *guard_beta, int figment, int guard_size, float *channel_x_after, float *channel_y_after, float *Lext_after, float LC, int blocks, int which_decoder, int trellis_term, int guarding_type)
{
float maxf(float a, float b);
int i,j,k;
float gamma_0, gamma_1, beta_0, beta_1;
float temp_beta[8];
int cpu_beta_state_0[8] 	= 	{0, 4, 5, 1, 2, 6, 7, 3};
int cpu_beta_state_1[8] 	= 	{4, 0, 1, 5, 6, 2, 3, 7};
int cpu_beta_encbit_0[8] 	= 	{-1, -1, 1, 1, 1, 1, -1, -1};
int cpu_beta_encbit_1[8] 	= 	{1, 1, -1, -1, -1, -1, 1, 1};
if(figment == blocks-1)
{
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
else
{
	if(guarding_type == 1)
	{
		
		for(i=0; i<8; i++)
		{
			start_beta[i] = guard_beta[(figment+1)*8 + i];
		}
		
	}
	else if(guarding_type == 2 || guarding_type == 3)
	{
		if(guarding_type == 2)
		{
			for(j=0; j<8; j++)
			{
				start_beta[j] = EQUAL_GUARD_VALUE;
			}
		}
		else if(guarding_type == 3)
		{
			for(j=0; j<8; j++)
			{
				start_beta[j] = guard_beta[(figment+1)*8 + j];
			}
		}
		for(k=guard_size-1; k>=0; k--)
		{
			for(j=0; j<8; j++)
			{
				gamma_0 = (0.5f)*(Lext_after[k] + LC*channel_x_after[k])*(-1) + (0.5f)*LC*channel_y_after[k]*cpu_beta_encbit_0[j];
				gamma_1 = (0.5f)*(Lext_after[k] + LC*channel_x_after[k])*(1) + (0.5f)*LC*channel_y_after[k]*cpu_beta_encbit_1[j];
				beta_0 = start_beta[cpu_beta_state_0[j] ] + gamma_0;
				beta_1 = start_beta[cpu_beta_state_1[j] ] + gamma_1;
				temp_beta[j] = maxf(beta_0, beta_1);
			}	
			for(j=0; j<8; j++)
			{
				start_beta[j] = temp_beta[j];
			}
		}	
	}
}
}


