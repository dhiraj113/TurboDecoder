//#include "cuPrintf.cu"
#include "decoderlog.cu"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hashdefined.h"
int inv_permutation_bits[DATASIZE];
main()
{
void decoder_types();
void simulations();
decoder_types();
//simulations();
system("rm decoder_output.dat encoded_bits.dat channel_out.dat input_bits.dat");
}

void simulations()
{
void simulation_iter(int numb_of_bits, int decoder_type, int guarding_type, int guard_size, int trellis_term_enable, int only_plotting);
void simulation_guard(int numb_of_bits, int decoder_type, int iter, int trellis_term_enable, int only_plotting);
void simulation_blocks(int numb_of_bits, int decoder_type, int iter, int guarding_type, int guard_size, int trellis_term_enable, int only_plotting);
int numb_of_bits = 1e8;
int decoder_type = 4;
int guarding_type = 3;
int guard_size = 8;
int iter = 5;
int trellis_term_enable = 0;
simulation_iter(numb_of_bits, decoder_type, guarding_type, guard_size, trellis_term_enable, 0);
simulation_guard(numb_of_bits, decoder_type, iter, trellis_term_enable, 0);
simulation_blocks(numb_of_bits, decoder_type, iter, guarding_type, guard_size, trellis_term_enable, 0);
}


void decoder_types()
{
void permuter_bits();
void encoder_and_noise(int numb_of_bits, int trellis_termination, float snr, int noise_enable);
void decode_and_analyse(int numb_of_bits, float snr, int iter_num, int blocks, int guard_size, int trellis_term_enable, int decoder_kind, int guarding_type, int *bit_error_count, int *frame_error_count);
float snr;
int numb_of_bits, iter, blocks, guard_size;
int trellis_term_enable, noise_enable;
int bit_error_count, frame_error_count;
trellis_term_enable = 0;
noise_enable = 1;
iter = 5;
snr = 1.5;
numb_of_bits = 1e7;
blocks = NO_OF_BLOCKS;
guard_size = 5;
permuter_bits();
encoder_and_noise(numb_of_bits, trellis_term_enable, snr, noise_enable);
decode_and_analyse(numb_of_bits, snr, iter, blocks, guard_size, trellis_term_enable, 4,1, &bit_error_count, &frame_error_count);
}



void decode_and_analyse(int numb_of_bits, float snr, int iter_num, int blocks, int guard_size, int trellis_term_enable, int decoder_kind, int guarding_type, int *bit_error_count, int *frame_error_count)
{
void full_decoder(int numb_of_bits, float snr, int iter_num, int blocks, int guard_size, int trellis_term, int decoder_kind, int guarding_type);
void analysis(int numb_of_bits, int *bit_error_count, int *frame_error_count);
full_decoder(numb_of_bits, snr, iter_num, blocks, guard_size, trellis_term_enable, decoder_kind, guarding_type);
analysis(numb_of_bits, bit_error_count, frame_error_count);
}



void full_decoder(int numb_of_bits, float snr, int iter_num, int blocks, int guard_size, int trellis_term, int decoder_kind, int guarding_type)
{
void just_decoder_basic(float snr, int iter_num, float *outbits_x, float *outbits_y1, float *outbits_y2, int *decision_bits, int trellis_term);
void just_decoder_parallel(float snr, int iter_num, int blocks, int guard_size, float *outbits_x, float *outbits_y1, float *outbits_y2, int *decision_bits, int trellis_term, int guarding_type);
void just_decoder_parallel_close(float snr, int iter_num, int blocks, int guard_size, float *outbits_x, float *outbits_y1, float *outbits_y2, int *decisionbits, int trellis_term, int guarding_type);
void just_decoder_gpu(float snr, int iter_num, int blocks, int guard_size, float *outbits_x, float *outbits_y1, float *outbits_y2, int *decision_bits, int trellis_term, int guarding_type);
clock_t start = clock(), hold1, hold2;
//*************Loading up decoder inputs to memory**************
numb_of_bits = ((numb_of_bits+DATASIZE-1)/DATASIZE)*DATASIZE;
FILE *fp, *gp;
fp = fopen("channel_out.dat", "r");
gp = fopen("decoder_output.dat", "w");
float *mega_outbits_x, *mega_outbits_y1, *mega_outbits_y2;
int *mega_decision_bits;
mega_outbits_x = (float*)malloc((numb_of_bits)*sizeof(float));
mega_outbits_y1 = (float*)malloc((numb_of_bits)*sizeof(float));
mega_outbits_y2 = (float*)malloc((numb_of_bits)*sizeof(float));
mega_decision_bits = (int*)malloc((numb_of_bits)*sizeof(int));
int numb;
int i;
for(i=0; i<numb_of_bits; i++)
{
	fscanf(fp,"%f\t%f\t%f\n", &mega_outbits_x[i], &mega_outbits_y1[i], &mega_outbits_y2[i]);
}
fclose(fp);
//*************Done Loading up decoder inputs to memory***************
//printf("Time elapsed for loading into memory: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
hold1 = clock();
//**************************Using cudahostalloc so speed up cudamemcpy**************************
float *outbits_x, *outbits_y1, *outbits_y2;
int *decisionbits;
cudaHostAlloc((void**)&outbits_x, DATASIZE*sizeof(float), cudaHostAllocDefault);
cudaHostAlloc((void**)&outbits_y1, DATASIZE*sizeof(float), cudaHostAllocDefault);
cudaHostAlloc((void**)&outbits_y2, DATASIZE*sizeof(float), cudaHostAllocDefault);
cudaHostAlloc((void**)&decisionbits, DATASIZE*sizeof(int), cudaHostAllocDefault);
//**************************Using cudahostalloc so speed up cudamemcpy**************************
for(numb = 0; numb< numb_of_bits/DATASIZE; numb++)
{
	for(i=0; i<DATASIZE; i++)
	{
		outbits_x[i] = mega_outbits_x[numb*DATASIZE + i];
		outbits_y1[i] = mega_outbits_y1[numb*DATASIZE + i];
		outbits_y2[i] = mega_outbits_y2[numb*DATASIZE + i];
	}
	//***********DATASIZE decoder******************
	switch(decoder_kind)
	{
		case 1 :{	just_decoder_basic(snr, iter_num, outbits_x, outbits_y1, outbits_y2, decisionbits, trellis_term);
				break;
			}
		case 2 :{	just_decoder_parallel(snr, iter_num, blocks, guard_size, outbits_x, outbits_y1, outbits_y2, decisionbits, trellis_term, guarding_type);
				break;
			}
		case 3 :{	just_decoder_parallel_close(snr, iter_num, blocks, guard_size, outbits_x, outbits_y1, outbits_y2, decisionbits, trellis_term, guarding_type);
				break;
			}
		case 4 :{	just_decoder_gpu(snr, iter_num, blocks, guard_size, outbits_x, outbits_y1, outbits_y2, decisionbits, trellis_term, guarding_type);
				break;
			}
		default:{
				printf("Wrong value give to decoder_kind\ncheck, exiting abruptly\n");
				exit(1);
			}
	}
	//***********DATASIZE decoder******************
	for(i=0; i<DATASIZE; i++)
	{
		mega_decision_bits[numb*DATASIZE+i] = decisionbits[i];
	}
}
hold2 = clock();
printf("Speed of decoding = %d Kbits/sec\n", (int)((numb_of_bits/((double)(hold2-hold1)/CLOCKS_PER_SEC))/1e3));
//printf("Time elapsed for decoding: %f\n", ((double)clock() - hold1) / CLOCKS_PER_SEC);
//printf("%d\n",(int)((numb_of_bits/((double)(hold2-hold1)/CLOCKS_PER_SEC))/1e3));
//**************************Freeing up cudahostalloc allocated memory**************************
cudaFreeHost(outbits_x);
cudaFreeHost(outbits_y1);
cudaFreeHost(outbits_y2);
cudaFreeHost(decisionbits);
//**************************Freeing up cudahostalloc allocated memory**************************
//**************Storing back the decoded bits to file*******************
for(i=0; i<numb_of_bits; i++)
{
	fprintf(gp, "%d\n",mega_decision_bits[i]);
}
free(mega_outbits_x);
free(mega_outbits_y1);
free(mega_outbits_y2);
free(mega_decision_bits);
fclose(gp);
//**************Done storing back the decoded bits to file**************
//printf("Time elapsed for writing to file: %f\n", ((double)clock() - hold2) / CLOCKS_PER_SEC);
}

void analysis(int numb_of_bits, int *bit_error_count, int *frame_error_count)
{
FILE *fp, *gp;
fp = fopen("input_bits.dat", "r");
gp = fopen("decoder_output.dat", "r");
int i,k;
int actual_bit, decision_bit;
int frame_error_flag, no_of_frames;
numb_of_bits = ((numb_of_bits+DATASIZE-1)/DATASIZE)*DATASIZE;
no_of_frames = numb_of_bits/DATASIZE;
*bit_error_count = 0;
*frame_error_count = 0;
for(k=0; k<no_of_frames; k++)
{
	frame_error_flag = 0;
	for(i=0; i<DATASIZE; i++)
	{
		fscanf(fp, "%d\n", &actual_bit);
		fscanf(gp, "%d\n", &decision_bit);
		if(actual_bit != decision_bit)
		{
			(*bit_error_count)++;
			frame_error_flag = 1;
			//printf("%d----%d\n", i%DATASIZE, i/DATASIZE);
		}
	}
	if(frame_error_flag == 1)
	{
		(*frame_error_count)++;
	}
}
fclose(fp);
fclose(gp);
printf("Bit_Error_count = %d\nFrame_Error_count=%d\n", *bit_error_count, *frame_error_count);
printf("BER = %g\nFER = %f\n", ((*bit_error_count)*1.0)/numb_of_bits, ((*frame_error_count)*1.0)/no_of_frames);
//printf("percentage correct = %f\n", ((numb_of_bits-error_count)*100.0)/numb_of_bits);
}


