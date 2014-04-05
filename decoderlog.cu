#include <stdio.h>
#include <stdlib.h>
#include "hashdefined.h"
//Most of the decoding would be done on the device
//Gamma evaluation on the device
__constant__ int deinter[DATASIZE];
__constant__ int alpha_state_0[8], alpha_state_1[8], alpha_encbit_0[8], alpha_encbit_1[8];
__constant__ int beta_state_0[8], beta_state_1[8], beta_encbit_0[8], beta_encbit_1[8];

void just_decoder_gpu(float snr, int iter_num, int blocks, int guard_size, float *outbits_x, float *outbits_y1, float *outbits_y2, int *decisionbits, int trellis_term, int guarding_type)
{
void result(float LC, float *outbits_x, int *decisionbits,float *Lext12, float *Lext21);
void setup_perm_bits();
void transition_array_creator();
__global__ void gpu_decoder(int blocks, int guard_size, float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y, float *alpha_guard_read, float *beta_guard_read, float *alpha_guard_write, float *beta_guard_write,  float LC, int which_decoder, int trellis_term, int guarding_type);
__global__ void Lext_initialiser(float *Lext12, float *Lext21);
__global__ void guard_initialiser(int blocks, float *alpha_guard_read_1, float *beta_guard_read_1,  float *alpha_guard_read_2,  float *beta_guard_read_2);
__global__ void sync_guard_values(int blocks, float *alpha_guard_read_1, float *alpha_guard_write_1, float *beta_guard_read_1, float *beta_guard_write_1, float *alpha_guard_read_2, float *alpha_guard_write_2, float *beta_guard_read_2, float *beta_guard_write_2 );
__global__ void result_ker(float LC, float *dev_outbits_x, int *dev_decisionbits, float *dev_Lext12, float *dev_Lext21);
extern int inv_permutation_bits[DATASIZE];
float LC;
int iter;
LC = (4.0/3.0)*pow(10,snr/10.0);
//Allocate memory on the device
//and transferring the required values
float *dev_outbits_x,  *dev_outbits_y1, *dev_outbits_y2;
float *dev_Lext12, *dev_Lext21;
int *dev_decisionbits;
float *dev_alpha_guard_read_1, *dev_alpha_guard_read_2, *dev_beta_guard_read_1, *dev_beta_guard_read_2;
float *dev_alpha_guard_write_1, *dev_alpha_guard_write_2, *dev_beta_guard_write_1, *dev_beta_guard_write_2;

//Allocating memory on the device for outbits and transferring data
cudaMalloc((void **)&dev_outbits_x,DATASIZE*sizeof(float));
cudaMalloc((void **)&dev_outbits_y1,DATASIZE*sizeof(float));
cudaMalloc((void **)&dev_outbits_y2,DATASIZE*sizeof(float));
cudaMemcpy(dev_outbits_x,outbits_x,DATASIZE*sizeof(float), cudaMemcpyHostToDevice);
cudaMemcpy(dev_outbits_y1,outbits_y1,DATASIZE*sizeof(float),cudaMemcpyHostToDevice);
cudaMemcpy(dev_outbits_y2,outbits_y2,DATASIZE*sizeof(float),cudaMemcpyHostToDevice);
//Allocating memory on the device for Lext12 and Lext21
cudaMalloc((void **)&dev_Lext12, DATASIZE*sizeof(float));
cudaMalloc((void **)&dev_Lext21, DATASIZE*sizeof(float));
//Allocating memory for decisionbits to be set at the end of all this
cudaMalloc((void**)&dev_decisionbits, DATASIZE*sizeof(int));
//Allocating memory on the device for alpha and beta guards
cudaMalloc((void**)&dev_alpha_guard_read_1, blocks*8*sizeof(float));
cudaMalloc((void**)&dev_alpha_guard_read_2, blocks*8*sizeof(float));
cudaMalloc((void**)&dev_beta_guard_read_1, blocks*8*sizeof(float));
cudaMalloc((void**)&dev_beta_guard_read_2, blocks*8*sizeof(float));
cudaMalloc((void**)&dev_alpha_guard_write_1, blocks*8*sizeof(float));
cudaMalloc((void**)&dev_alpha_guard_write_2, blocks*8*sizeof(float));
cudaMalloc((void**)&dev_beta_guard_write_1, blocks*8*sizeof(float));
cudaMalloc((void**)&dev_beta_guard_write_2, blocks*8*sizeof(float));
//Finished allocating memory on the device
Lext_initialiser<<<DATASIZE/64,64>>>(dev_Lext12, dev_Lext21);
setup_perm_bits();
transition_array_creator();
guard_initialiser<<<(blocks*8)/64, 64>>>(blocks, dev_alpha_guard_read_1, dev_beta_guard_read_1,  dev_alpha_guard_read_2,  dev_beta_guard_read_2);
//The iterative decoding is done here
//With repeated calls to kernels on the gpu
//cudaPrintfInit();
for(iter =1; iter<= iter_num; iter++)
{	
	gpu_decoder<<<blocks,8>>>(blocks, guard_size, dev_Lext21, dev_Lext12, dev_outbits_x, dev_outbits_y1, dev_alpha_guard_read_1, dev_beta_guard_read_1, dev_alpha_guard_write_1, dev_beta_guard_write_1, LC, 1, trellis_term, guarding_type);
	gpu_decoder<<<blocks,8>>>(blocks, guard_size, dev_Lext12, dev_Lext21, dev_outbits_x, dev_outbits_y2, dev_alpha_guard_read_2, dev_beta_guard_read_2, dev_alpha_guard_write_2, dev_beta_guard_write_2, LC, 2, trellis_term, guarding_type);
	if(guarding_type != 2)
	{
		sync_guard_values<<<(blocks*8)/64, 64>>>(blocks, dev_alpha_guard_read_1, dev_alpha_guard_write_1, dev_beta_guard_read_1, dev_beta_guard_write_1, dev_alpha_guard_read_2, dev_alpha_guard_write_2, dev_beta_guard_read_2, dev_beta_guard_write_2);
	}
	//cudaPrintfDisplay(stdout, true);
}
//cudaPrintfEnd();
result_ker<<<48,128>>>(LC, dev_outbits_x, dev_decisionbits, dev_Lext12, dev_Lext21);
cudaMemcpy(decisionbits, dev_decisionbits, DATASIZE*sizeof(int), cudaMemcpyDeviceToHost);

/********Lext12 and Lext21 check block**************
FILE *fp;
int pi;
fp = fopen("Lext_device.dat", "w");
for(pi=0; pi<DATASIZE; pi++)
{
	fprintf(fp, "%f\t%f\t%f\t%f\n", LC*outbits_x[pi], Lext12[pi], Lext21[pi], LC*outbits_x[pi] + Lext12[pi] + Lext21[pi]);
}
fclose(fp);
********Lext12 and Lext21 check block**************/


//Freeing the allocated memory on the device
cudaFree(dev_Lext12);
cudaFree(dev_Lext21);
cudaFree(dev_outbits_x);
cudaFree(dev_outbits_y1);
cudaFree(dev_outbits_y2);
cudaFree(dev_alpha_guard_read_1);
cudaFree(dev_alpha_guard_read_2);
cudaFree(dev_beta_guard_read_1);
cudaFree(dev_beta_guard_read_2);
cudaFree(dev_alpha_guard_write_1);
cudaFree(dev_alpha_guard_write_2);
cudaFree(dev_beta_guard_write_1);
cudaFree(dev_beta_guard_write_2);
cudaFree(dev_decisionbits);
//Done freeing memory on the device
}

__global__ void gpu_decoder(int blocks, int guard_size, float *Lext_in, float *Lext_out, float *outbits_x, float *outbits_y, float *alpha_guard_read, float *beta_guard_read, float *alpha_guard_write, float *beta_guard_write, float LC, int which_decoder, int trellis_term, int guarding_type)
{
__device__ void guard_alpha_gpu(int blocks, int guard_size, float *start_alpha, float *alpha_guard, float *channel_x_before, float *channel_y_before, float *Lext_before, float LC, int guarding_type);
__device__ void guard_beta_gpu(int blocks, int guard_size, float *start_beta, float *beta_guard, float *channel_x_after, float *channel_y_after, float *Lext_after, float LC, int which_decoder, int trellis_term, int guarding_type);
__device__ int inter(int k);
__device__ float maxstar(float a,float b);
int tid = threadIdx.x;
int bid = blockIdx.x;
int i,n, j, k, l, index, n_eq;
float gamma_0, gamma_1, alpha_0, alpha_1;
float gamma_00, gamma_11, gamma_000, gamma_111;
float beta_0, beta_1;
float Lext_0, Lext_1, Lexternal;


int block_size = DATASIZE/blocks;

//int block_start = bid*block_size;
//int block_end = (bid+1)*block_size - 1;

__shared__ float channel_x[MAX_SUB_BLOCK_SIZE+1], channel_y[MAX_SUB_BLOCK_SIZE+1];
__shared__ float alpha[(MAX_SUB_BLOCK_SIZE+1)*8], Lext[MAX_SUB_BLOCK_SIZE+1];  //520+65*3 715 floats
__shared__ float lambda_0[8*8], lambda_1[8*8]; //128 floats
__shared__ float start_alpha[8], start_beta[8], temp_beta[8], beta[8]; //32 floats
__shared__ float channel_x_before[MAX_GUARD_SIZE], channel_x_after[MAX_GUARD_SIZE], channel_y_before[MAX_GUARD_SIZE], channel_y_after[MAX_GUARD_SIZE]; //64 floats
__shared__ float  Lext_before[MAX_GUARD_SIZE], Lext_after[MAX_GUARD_SIZE]; //32 floats
__shared__ float alpha_guard_copy[8], beta_guard_copy[8]; //16 float

__shared__ float gamma_e1[8], gamma_e2[8];



//************Fecthing into shared memory for the purpose of guarding************
if(guarding_type != 2)
{
	if(bid != 0)
	{
		alpha_guard_copy[tid] = alpha_guard_read[(bid-1)*8+tid];
	}
	if(bid != blocks-1)
	{
		beta_guard_copy[tid] = beta_guard_read[(bid+1)*8+tid];
	}
}
if(guarding_type != 1)
{
	if(bid != 0)
	{
		for(i=0; i<MAX_GUARD_SIZE/8; i++)
		{
			channel_y_before[i*8+tid] = outbits_y[bid*block_size - MAX_GUARD_SIZE + i*8+tid];
		}
		if(which_decoder == 1)
		{
			for(i=0; i<MAX_GUARD_SIZE/8; i++)
			{
				channel_x_before[i*8+tid] = outbits_x[bid*block_size - MAX_GUARD_SIZE + i*8+tid];
				Lext_before[i*8+tid] = Lext_in[bid*block_size - MAX_GUARD_SIZE + i*8+tid];
			}
		}
		else if(which_decoder == 2)
		{
			for(i=0; i<MAX_GUARD_SIZE/8; i++)
			{
				channel_x_before[i*8+tid] = outbits_x[deinter[bid*block_size - MAX_GUARD_SIZE + i*8+tid]];
				Lext_before[i*8+tid] = Lext_in[deinter[bid*block_size - MAX_GUARD_SIZE + i*8+tid]];
			}	
		}
	}
	if(bid != blocks-1)
	{
		for(i=0; i<MAX_GUARD_SIZE/8; i++)
		{
			channel_y_after[i*8+tid] = outbits_y[(bid+1)*block_size +1 + i*8+tid];
		}
		if(which_decoder == 1)
		{
			for(i=0; i<MAX_GUARD_SIZE/8; i++)
			{
				channel_x_after[i*8+tid] = outbits_x[(bid+1)*block_size +1 + i*8+tid];
				Lext_after[i*8+tid] = Lext_in[(bid+1)*block_size +1 + i*8+tid];
			}
		}
		else if(which_decoder == 2)
		{
			for(i=0; i<MAX_GUARD_SIZE/8; i++)
			{
				channel_x_after[i*8+tid] = outbits_x[deinter[(bid+1)*block_size +1 + i*8+tid]];
				Lext_after[i*8+tid] = Lext_in[deinter[(bid+1)*block_size +1 + i*8+tid]];
			}	
		}
	}
}
//************Fecthing into shared memory for the purpose of guarding************
guard_alpha_gpu(blocks, guard_size, start_alpha, alpha_guard_copy, channel_x_before, channel_y_before, Lext_before, LC, guarding_type);
guard_beta_gpu(blocks, guard_size, start_beta, beta_guard_copy, channel_x_after, channel_y_after, Lext_after, LC, which_decoder, trellis_term, guarding_type);
alpha[tid]=start_alpha[tid];
beta[tid] = start_beta[tid];

//*******Performing coalesced memory access to fetch channel values into shared memory******
for(i=0; i<block_size/8; i++)
{
	index = bid*block_size+i*8+tid;
	channel_y[i*8+tid] = outbits_y[index];
}
//*******Performing coalesced memory access to fetch channel values into shared memory******

//********Loading Lext_in appropriately, depending on the which_decoder*********************
for(i=0; i<block_size/8; i++)
{
	index = bid*block_size+i*8+tid;
	if(which_decoder == 1)
	{
		channel_x[i*8+tid] = outbits_x[index];	
		Lext[i*8+tid] = Lext_in[index];
	}
	else if(which_decoder == 2)
	{
		channel_x[i*8+tid] = outbits_x[deinter[index]];
		Lext[i*8+tid] = Lext_in[deinter[index]];
	}
}
//********Loading Lext_in appropriately, depending on the which_decoder*********************

//***********Fetching the values corresponding to first bit of the next block****************
if(bid != blocks-1)
{
	index = (bid+1)*block_size;
	channel_y[block_size] = outbits_y[index];
	if(which_decoder == 1)
	{
		channel_x[block_size] = outbits_x[index];
		Lext[block_size] = Lext_in[index];
	}
	else if(which_decoder == 2)
	{
		channel_x[block_size] = outbits_x[deinter[index]];
		Lext[block_size] = Lext_in[deinter[index]];
	}
}
//***********Fetching the values corresponding to first bit of the next block****************

//*************************************Alpha evaluation*******************************


/*for(i=0; i<block_size; i++)
{	
	
	gamma_0 = (0.5f)*(LC*channel_x[i]+Lext[i])*(-1) + (0.5f)*LC*channel_y[i]*alpha_encbit_0[tid];
	alpha_0 = alpha[ i*8 + alpha_state_0[tid] ] + gamma_0;
	gamma_1 = (0.5f)*(LC*channel_x[i]+Lext[i])*(1) + (0.5f)*LC*channel_y[i]*alpha_encbit_1[tid];
	alpha_1 = alpha[ i*8 + alpha_state_1[tid] ] + gamma_1;
	alpha[ (i+1)*8 + tid ] = maxstar(alpha_0, alpha_1);
	//__syncthreads();
}*/

for(k=0; k<block_size/8; k++)
{
	gamma_e1[tid] = (0.5f)*(LC*channel_x[k*8+tid]+Lext[k*8+tid]);
	gamma_e2[tid] = (0.5f)*LC*channel_y[k*8+tid];
	for(j=0; j<8; j++)
	{	
		i = k*8+j;
		gamma_0 = gamma_e1[j]*(-1) + gamma_e2[j]*alpha_encbit_0[tid];
		gamma_1 = gamma_e1[j] + gamma_e2[j]*alpha_encbit_1[tid];
		alpha_0 = alpha[ i*8 + alpha_state_0[tid] ] + gamma_0;
		alpha_1 = alpha[ i*8 + alpha_state_1[tid] ] + gamma_1;
		alpha[ (i+1)*8 + tid ] = maxstar(alpha_0, alpha_1);
	//__syncthreads();
	}
}

//***************************Writing end alpha for use in the next iteration***************************
if(guarding_type == 1)
{
	alpha_guard_write[bid*8+tid] = alpha[block_size*8+tid];
}
else if(guarding_type == 3)
{
	alpha_guard_write[bid*8+tid] = alpha[(block_size-guard_size)*8+tid];
}

//***************************Writing end alpha for use in the next iteration***************************
//********************Alpha print check********************
//cuPrintf("%f\n", alpha[(2)*8+tid]);
//********************Alpha print check********************
//*************************************Alpha evaluation*******************************

//******************************Beta and Lext evaluation*******************************
if( bid != blocks-1)
{
	gamma_000 = (0.5f)*(LC*channel_x[block_size]+Lext[block_size])*(-1) + (0.5f)*LC*channel_y[block_size]*beta_encbit_0[tid];
	gamma_111 = (0.5f)*(LC*channel_x[block_size]+Lext[block_size])*(1) + (0.5f)*LC*channel_y[block_size]*beta_encbit_1[tid];
	for(k=block_size/8; k>0; k--)
	{
		gamma_e1[tid] = (0.5f)*(LC*channel_x[k*8-tid-1]+Lext[k*8-tid-1]);
		gamma_e2[tid] = (0.5f)*LC*channel_y[k*8-tid-1];
		for(j=0; j<=7; j++)
		{
			n = k*8-j-1;
			//**********************************Beta evaluation****************************
			//gamma_000 = (0.5f)*(LC*channel_x[n+1]+Lext[n+1])*(-1) + (0.5f)*LC*channel_y[n+1]*beta_encbit_0[tid];
			//gamma_111 = (0.5f)*(LC*channel_x[n+1]+Lext[n+1])*(1) + (0.5f)*LC*channel_y[n+1]*beta_encbit_1[tid];
			beta_0 = beta[beta_state_0[tid]] + gamma_000;
			beta_1 = beta[beta_state_1[tid]] + gamma_111;
			temp_beta[tid] = maxstar(beta_0, beta_1);
			beta[tid] = temp_beta[tid];			 
				
			//**********************************Lambda evaluation**************************		
			//gamma_000 = (0.5f)*(LC*channel_x[n]+Lext[n])*(-1) + (0.5f)*LC*channel_y[n]*beta_encbit_0[tid];
			//gamma_111 = (0.5f)*(LC*channel_x[n]+Lext[n])*(1) + (0.5f)*LC*channel_y[n]*beta_encbit_1[tid];
			gamma_000 = gamma_e1[j]*(-1) + gamma_e2[j]*beta_encbit_0[tid];
			gamma_111 = gamma_e1[j] + gamma_e2[j]*beta_encbit_1[tid];
			lambda_0[j*8+tid] = gamma_000 + beta[beta_state_0[tid]] + alpha[n*8+tid];
			lambda_1[j*8+tid] = gamma_111 + beta[beta_state_1[tid]] + alpha[n*8+tid];
			//**********************************Lambda evaluation***************************
			
			//***********************Writing beta for initialisation in the next iteration*************************
			if(guarding_type == 3)
			{
				if(n == guard_size)
				{
					beta_guard_write[bid*8+tid] = beta[tid];
				}
			}	
			//***********************Writing beta for initialisation in the next iteration*************************
			//**********************************Beta evaluation****************************
		
			//****************Print lambda check**********************
			//if(n==block_size-1) cuPrintf("%f\t%f\n", lambda_0[j*8+tid], lambda_1[j*8+tid]);
			//****************Print lambda check**********************
			//********************Beta print check********************
			//if(n== 0) cuPrintf("%f\n", beta[tid]);
			//********************Beta print check********************		
		}
	
		//****************************Lexternal evaluation****************************
		
		Lext_0 = lambda_0[tid*8];
		Lext_1 = lambda_1[tid*8];
		//Lext_0 = 0;
		//Lext_1 = 0;
		for(l=1; l<8; l++)
		{
			//index = tid*8+(l+tid)%8;
			//index = tid*8+(l+tid)&7;
			index = tid*8 + l;
			Lext_0 = maxstar(Lext_0, lambda_0[index]);
			Lext_1 = maxstar(Lext_1, lambda_1[index]);
		}
		n_eq = k*8 - tid - 1;
		Lexternal = (Lext_1 - Lext_0) - LC*channel_x[n_eq] - Lext[n_eq];
		
		//****************Full Lext print check**********************
		//if(k==1) cuPrintf("%f\t%f\t%f\t%f\t%f\n", channel_x[k*8 - tid - 1],Lexternal, Lext_0, Lext_1, Lext[k*8 - tid - 1]);
		//if(k==8) cuPrintf("%f\n", Lexternal);
		//****************Full Lext print check**********************
		
		index = bid*block_size + n_eq;
		if(which_decoder == 1)
		{
			Lext_out[index] = Lexternal; //30 cycles
		}
		else if(which_decoder == 2)
		{
			Lext_out[deinter[index]] = Lexternal; //40 cycles
		}
		//****************************Lexternal evaluation****************************
		//__syncthreads();
		
	}
	//***********************Writing beta for initialisation in the next iteration*************************
	if(guarding_type == 1)
	{
		beta_guard_write[bid*8+tid] = beta[tid];
	}
	//***********************Writing beta for initialisation in the next iteration*************************
	//******************************Beta and Lext evaluation*******************************
}
else
{
	//gamma_00 = (0.5f)*(LC*channel_x[block_size-1]+Lext[block_size-1])*(-1) + (0.5f)*LC*channel_y[block_size-1]*beta_encbit_0[tid];
	//gamma_11 = (0.5f)*(LC*channel_x[block_size-1]+Lext[block_size-1])*(1) + (0.5f)*LC*channel_y[block_size-1]*beta_encbit_1[tid];
	for(k=block_size/8; k>0; k--)
	{
		for(j=0; j<=7; j++)
		{
			n = k*8-j-1;
			if(n == block_size-1)
			{
				temp_beta[tid] = EQUAL_GUARD_VALUE;
				beta[tid] = EQUAL_GUARD_VALUE;
			}
			else
			{
				//**********************************Beta evaluation****************************
				gamma_00 = (0.5f)*(LC*channel_x[n+1]+Lext[n+1])*(-1) + (0.5f)*LC*channel_y[n+1]*beta_encbit_0[tid];
				gamma_11 = (0.5f)*(LC*channel_x[n+1]+Lext[n+1])*(1) + (0.5f)*LC*channel_y[n+1]*beta_encbit_1[tid];
				beta_0 = beta[beta_state_0[tid]] + gamma_00;
				beta_1 = beta[beta_state_1[tid]] + gamma_11;
				temp_beta[tid] = maxstar(beta_0, beta_1);
				beta[tid] = temp_beta[tid];			 
			}
						
			//**********************************Lambda evaluation**************************		
			gamma_00 = (0.5f)*(LC*channel_x[n]+Lext[n])*(-1) + (0.5f)*LC*channel_y[n]*beta_encbit_0[tid];
			gamma_11 = (0.5f)*(LC*channel_x[n]+Lext[n])*(1) + (0.5f)*LC*channel_y[n]*beta_encbit_1[tid];
			lambda_0[j*8+tid] = gamma_00 + beta[beta_state_0[tid]] + alpha[n*8+tid];
			lambda_1[j*8+tid] = gamma_11 + beta[beta_state_1[tid]] + alpha[n*8+tid];
			//**********************************Lambda evaluation***************************
			
		
			//***********************Writing beta for initialisation in the next iteration*************************
			if(guarding_type == 3)
			{
				if(n == guard_size)
				{
					beta_guard_write[bid*8+tid] = beta[tid];
				}
			}	
			//***********************Writing beta for initialisation in the next iteration*************************
			//**********************************Beta evaluation****************************
		
			//****************Print lambda check**********************
			//if(n==block_size-1) cuPrintf("%f\t%f\n", lambda_0[j*8+tid], lambda_1[j*8+tid]);
			//****************Print lambda check**********************
			//********************Beta print check********************
			//if(n== 0) cuPrintf("%f\n", beta[tid]);
			//********************Beta print check********************		
		}
	
		//****************************Lexternal evaluation****************************
		
		Lext_0 = lambda_0[tid*8];
		Lext_1 = lambda_1[tid*8];
		//Lext_0 = 0;
		//Lext_1 = 0;
		for(l=1; l<8; l++)
		{
			//index = tid*8+(l+tid)%8;
			//index = tid*8+(l+tid)&7;
			index = tid*8 + l;
			Lext_0 = maxstar(Lext_0, lambda_0[index]);
			Lext_1 = maxstar(Lext_1, lambda_1[index]);
		}
		n_eq = k*8 - tid - 1;
		Lexternal = (Lext_1 - Lext_0) - LC*channel_x[n_eq] - Lext[n_eq];
		
		//****************Full Lext print check**********************
		//if(k==1) cuPrintf("%f\t%f\t%f\t%f\t%f\n", channel_x[k*8 - tid - 1],Lexternal, Lext_0, Lext_1, Lext[k*8 - tid - 1]);
		//if(k==8) cuPrintf("%f\n", Lexternal);
		//****************Full Lext print check**********************
		
		index = bid*block_size + n_eq;
		if(which_decoder == 1)
		{
			Lext_out[index] = Lexternal; //30 cycles
		}
		else if(which_decoder == 2)
		{
			Lext_out[deinter[index]] = Lexternal; //40 cycles
		}
		//****************************Lexternal evaluation****************************
		//__syncthreads();
		
	}
	//***********************Writing beta for initialisation in the next iteration*************************
	if(guarding_type == 1)
	{
		beta_guard_write[bid*8+tid] = beta[tid];
	}
	//***********************Writing beta for initialisation in the next iteration*************************
	//******************************Beta and Lext evaluation*******************************
}
}

__device__ void guard_alpha_gpu(int blocks, int guard_size, float *start_alpha, float *alpha_guard, float *channel_x_before, float *channel_y_before, float *Lext_before, float LC, int guarding_type)
{
__device__ float maxstar(float a,float b);
int tid = threadIdx.x;
int bid = blockIdx.x;
__shared__ float temp_alpha[8];
int k;
float gamma_0, gamma_1, alpha_0, alpha_1;
if( bid != 0)
{
	if(guarding_type == 1)
	{
		start_alpha[tid] = alpha_guard[tid];
	}
	else if(guarding_type == 2 || guarding_type == 3)
	{
		if(guarding_type == 2)
		{
			start_alpha[tid] = EQUAL_GUARD_VALUE;
		}
		else if(guarding_type == 3)
		{
			start_alpha[tid] = alpha_guard[tid];
		}
		for(k=MAX_GUARD_SIZE-guard_size; k<MAX_GUARD_SIZE; k++)
		{
			
			gamma_0 = (0.5f)*(Lext_before[k] + LC*channel_x_before[k])*(-1) + (0.5f)*LC*channel_y_before[k]*alpha_encbit_0[tid];
			gamma_1 = (0.5f)*(Lext_before[k] + LC*channel_x_before[k])*(1) + (0.5f)*LC*channel_y_before[k]*alpha_encbit_1[tid];
			alpha_0 = start_alpha[alpha_state_0[tid] ] + gamma_0;
			alpha_1 = start_alpha[alpha_state_1[tid] ] + gamma_1;
			temp_alpha[tid] = maxstar(alpha_0, alpha_1);
			//__syncthreads();
			start_alpha[tid] = temp_alpha[tid];
			//__syncthreads();
		}
	}
}
else
{
	if(tid == 0)
	{
		start_alpha[0] = 0;
	}
	else
	{
		start_alpha[tid] = MINUS_INFINITY;
	}
}
}


__device__ void guard_beta_gpu(int blocks, int guard_size, float *start_beta, float *beta_guard, float *channel_x_after, float *channel_y_after, float *Lext_after, float LC, int which_decoder, int trellis_term, int guarding_type)
{
__device__ float maxstar(float a,float b);
int tid = threadIdx.x;
int bid = blockIdx.x;
int k;
float gamma_0, gamma_1, beta_0, beta_1;
__shared__ float temp_beta[8];
if(bid != blocks-1)
{
	if(guarding_type == 1)
	{
		start_beta[tid] = beta_guard[tid];
		
	}
	else if(guarding_type == 2 || guarding_type == 3)
	{
		if(guarding_type == 2)
		{
			start_beta[tid] = EQUAL_GUARD_VALUE;
		}
		else if(guarding_type == 3)
		{
			start_beta[tid] = beta_guard[tid];
		}
		for(k=guard_size-1; k>=0; k--)
		{
			gamma_0 = (0.5f)*(Lext_after[k] + LC*channel_x_after[k])*(-1) + (0.5f)*LC*channel_y_after[k]*beta_encbit_0[tid];
			gamma_1 = (0.5f)*(Lext_after[k] + LC*channel_x_after[k])*(1) + (0.5f)*LC*channel_y_after[k]*beta_encbit_1[tid];
			beta_0 = start_beta[beta_state_0[tid]] + gamma_0;
			beta_1 = start_beta[beta_state_1[tid]] + gamma_1;
			temp_beta[tid] = maxstar(beta_0, beta_1);
			//__syncthreads();
			start_beta[tid] = temp_beta[tid];
			//__syncthreads();
		}	
	}	
}
else
{
	if(which_decoder == 1 && trellis_term == 1)
	{
		if(tid == 0)
		{
			start_beta[0] = 0;
		}
		else
		{
			start_beta[tid] = MINUS_INFINITY;
		}
	}
	else
	{
		start_beta[tid] = EQUAL_GUARD_VALUE;
	}	
}
}



void transition_array_creator()
{
	int cpu_alpha_state_0[8] 	= 	{0, 3, 4, 7, 1, 2, 5, 6};
	int cpu_alpha_state_1[8] 	= 	{1, 2, 5, 6, 0, 3, 4, 7};
	int cpu_beta_state_0[8] 	= 	{0, 4, 5, 1, 2, 6, 7, 3};
	int cpu_beta_state_1[8] 	= 	{4, 0, 1, 5, 6, 2, 3, 7};
	int cpu_alpha_encbit_0[8] 	= 	{-1, 1, 1, -1, -1, 1, 1, -1};
	int cpu_alpha_encbit_1[8] 	= 	{1, -1, -1, 1, 1, -1, -1, 1};
	int cpu_beta_encbit_0[8] 	= 	{-1, -1, 1, 1, 1, 1, -1, -1};
	int cpu_beta_encbit_1[8] 	= 	{1, 1, -1, -1, -1, -1, 1, 1};
	cudaMemcpyToSymbol(alpha_state_0, cpu_alpha_state_0, 8*sizeof(int));
	cudaMemcpyToSymbol(alpha_state_1, cpu_alpha_state_1, 8*sizeof(int));
	cudaMemcpyToSymbol(beta_state_0, cpu_beta_state_0, 8*sizeof(int));
	cudaMemcpyToSymbol(beta_state_1, cpu_beta_state_1, 8*sizeof(int));
	cudaMemcpyToSymbol(alpha_encbit_0, cpu_alpha_encbit_0, 8*sizeof(int));
	cudaMemcpyToSymbol(alpha_encbit_1, cpu_alpha_encbit_1, 8*sizeof(int));
	cudaMemcpyToSymbol(beta_encbit_0, cpu_beta_encbit_0, 8*sizeof(int));
	cudaMemcpyToSymbol(beta_encbit_1, cpu_beta_encbit_1, 8*sizeof(int));	
}

void result(float LC, float *outbits_x, int *decisionbits, float *Lext12, float *Lext21)
{
float L1uk;
int k;
for(k=0; k<DATASIZE; k++)
{
	L1uk = LC*outbits_x[k]+ Lext12[k] + Lext21[k];
	//L1uk =Lext12[k];
	//L1uk =Lext12[k];
	//L1uk = Lext12[k] + Lext21[k];
	decisionbits[k] = L1uk > 0 ? 1 : 0;	
}
}


void setup_perm_bits()
{
	extern int inv_permutation_bits[DATASIZE];
	cudaMemcpyToSymbol(deinter, inv_permutation_bits, DATASIZE*sizeof(int));
}



__device__ int inter(int n)
{
	return (((1+6*n)%DATASIZE)*n)%DATASIZE;
}


__device__ float maxstar(float a,float b)
{
	return max(a,b);
	//return max(a,b) + __logf(1+__expf((-1)*abs(a-b)));
}

__global__ void sync_guard_values(int blocks, float *alpha_guard_read_1, float *alpha_guard_write_1, float *beta_guard_read_1, float *beta_guard_write_1, float *alpha_guard_read_2, float *alpha_guard_write_2, float *beta_guard_read_2, float *beta_guard_write_2 )
{
//Copy write values to the read values
int index = threadIdx.x + blockIdx.x*blockDim.x;
if(index < blocks*8)
{
	alpha_guard_read_1[index] = alpha_guard_write_1[index];
	beta_guard_read_1[index] = beta_guard_write_1[index];
	alpha_guard_read_2[index] = alpha_guard_write_2[index];
	beta_guard_read_2[index] = beta_guard_write_2[index];
}
}


__global__ void guard_initialiser(int blocks, float *alpha_guard_read_1, float *beta_guard_read_1,  float *alpha_guard_read_2,  float *beta_guard_read_2)
{
//Copy write values to the read values
int index = threadIdx.x + blockIdx.x*blockDim.x;
if(index < blocks*8)
{
	alpha_guard_read_1[index] = EQUAL_GUARD_VALUE;
	beta_guard_read_1[index] = EQUAL_GUARD_VALUE;
	alpha_guard_read_2[index] = EQUAL_GUARD_VALUE;
	beta_guard_read_2[index] = EQUAL_GUARD_VALUE;
}
}


__global__ void Lext_initialiser(float *Lext12, float *Lext21)
{
int index = threadIdx.x + blockIdx.x*blockDim.x;
if(index < DATASIZE)
{
	Lext12[index] = 0;
	Lext21[index] = 0;
}
}
__global__ void result_ker(float LC, float *dev_outbits_x, int *dev_decisionbits, float *dev_Lext12, float *dev_Lext21)
{
float L1uk;
int index = threadIdx.x + blockIdx.x*blockDim.x;
if(index < DATASIZE)
{
	L1uk = LC*dev_outbits_x[index] + dev_Lext12[index] + dev_Lext21[index];
	dev_decisionbits[index] = L1uk > 0 ? 1 : 0;	
}
}
