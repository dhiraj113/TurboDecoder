#include <stdio.h>
#include <stdlib.h>
#include "hashdefined.h"
void simulation_iter(int numb_of_bits, int decoder_type, int guarding_type, int guard_size, int trellis_term_enable, int only_plotting)
{
//Plots BER vs SNR in logscale for different iterations
//Input information to be given is type of guarding, guard_size
void permuter_bits();
void encoder_and_noise(int numb_of_bits, int trellis_termination, float snr, int noise_enable);
int decode_and_analyse(int numb_of_bits, float snr, int iter_num, int blocks, int guard_size, int trellis_term_enable, int decoder_kind, int guarding_type);
void plotting_1(int size, int *iter_array, int decoder_type, int guarding_type, int guard_size);
float snr;
float snr_start = 0.0;
float snr_end = 2.005;
float resolution = 0.1;
int iter_array[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
int size = 10;
int i,k, error_count, numb_snrs = 0;
int *errordata;
FILE *fp, *gp;
if(only_plotting != 1)
{
	fp = fopen("ber_vs_snr_iters.dat", "w");
	gp = fopen("ber_vs_snr_iters_linear.dat", "w");
	permuter_bits();
	for(snr = snr_start; snr<=snr_end; snr += resolution)
	{
		encoder_and_noise(numb_of_bits, 0, snr, 1);
		fprintf(fp,"%f\t", snr);
		numb_snrs++;
		for(i = 0; i<size; i++)
		{
			error_count = decode_and_analyse(numb_of_bits, snr, iter_array[i], NO_OF_BLOCKS, guard_size, trellis_term_enable, decoder_type, guarding_type);
			fprintf(fp, "%d\t", error_count);
			fprintf(gp, "%d\n", error_count);
		}
		fprintf(fp, "\n");
		printf("Done for snr = %f\n", snr);
	}
	fclose(fp);
	fclose(gp);
	//Converting the data into a format that can be used by gnuplot
	gp = fopen("ber_vs_snr_iters_linear.dat", "r");
	errordata = (int*)malloc(numb_snrs*size*sizeof(int));
	for(i=0; i<numb_snrs*size; i++)
	{
		fscanf(gp, "%d\n", &errordata[i]);
	}
	fclose(gp);
	fp = fopen("ber_vs_snr_iters_gnu.dat", "w");
	for(k=0; k<size; k++)
	{
		fprintf(fp, "#####Iteration %d#####\n", iter_array[k]);
		for(i=0; i<numb_snrs; i++)
		{
			fprintf(fp,"%f\t%d\t%g\n", snr_start + i*resolution, errordata[k+i*size], errordata[k+i*size]/(numb_of_bits*(1.0)));
		}
		fprintf(fp, "\n\n\n\n");
	}
	fclose(fp);
	free(errordata);
}
plotting_1(size, iter_array, decoder_type, guarding_type, guard_size);
}


void simulation_blocks(int numb_of_bits, int decoder_type, int iter, int guarding_type, int guard_size, int trellis_term_enable, int only_plotting)
{
void permuter_bits();
void encoder_and_noise(int numb_of_bits, int trellis_termination, float snr, int noise_enable);
int decode_and_analyse(int numb_of_bits, float snr, int iter_num, int blocks, int guard_size, int trellis_term_enable, int decoder_kind, int guarding_type);
void plotting_3(int size, int *blocks_array, int decoder_type, int iter, int guarding_type, int guard_size);
float snr;
float snr_start = 0.0;
float snr_end = 2.005;
float resolution = 0.1;
int blocks_array[5] = {96, 128, 192};
int size = 3;
int i,k, error_count, numb_snrs = 0;
int *errordata;
FILE *fp, *gp;
if(only_plotting != 1)
{
	fp = fopen("ber_vs_snr_blocks.dat", "w");
	gp = fopen("ber_vs_snr_blocks_linear.dat", "w");
	permuter_bits();
	for(snr = snr_start; snr<=snr_end; snr += resolution)
	{
		encoder_and_noise(numb_of_bits, 0, snr, 1);
		fprintf(fp,"%f\t", snr);
		numb_snrs++;		
		for(i = 0; i<size; i++)
		{
			error_count = decode_and_analyse(numb_of_bits, snr, iter, blocks_array[i], guard_size, trellis_term_enable, decoder_type, guarding_type);
			fprintf(fp, "%d\t", error_count);
			fprintf(gp, "%d\n", error_count);
		}
		error_count = decode_and_analyse(numb_of_bits, snr, iter, NO_OF_BLOCKS, 0, trellis_term_enable, 1, 1); //on the cpu :(
		fprintf(fp, "%d\t", error_count);
		fprintf(gp, "%d\n", error_count);
		fprintf(fp, "\n");
		printf("Done for snr = %f\n", snr);
	}
	fclose(fp);
	fclose(gp);
	//Converting the data into a format that can be used by gnuplot
	gp = fopen("ber_vs_snr_blocks_linear.dat", "r");
	errordata = (int*)malloc(numb_snrs*size*sizeof(int));
	for(i=0; i<numb_snrs*size; i++)
	{
		fscanf(gp, "%d\n", &errordata[i]);
	}
	fclose(gp);
	fp = fopen("ber_vs_snr_blocks_gnu.dat", "w");
	for(k=0; k<size+1; k++)
	{
		if(k == size) fprintf(fp,"#####blocks 1#####\n");
		else fprintf(fp, "#####blocks %d#####\n", blocks_array[k]);
		for(i=0; i<numb_snrs; i++)
		{
			fprintf(fp,"%f\t%d\t%g\n", snr_start + i*resolution, errordata[k+i*size], errordata[k+i*size]/(numb_of_bits*(1.0)));
		}
		fprintf(fp, "\n\n\n\n");
	}
	fclose(fp);
	free(errordata);
}
plotting_3(size, blocks_array, decoder_type,iter, guarding_type, guard_size);
}

void simulation_guard(int numb_of_bits, int decoder_type, int iter, int trellis_term_enable, int only_plotting)
{
//Plots BER vs SNR in logscale for different types of guarding and corresponding different guard sizes
//Input information to be given is the number of iterations to be be performed
void permuter_bits();
void encoder_and_noise(int numb_of_bits, int trellis_termination, float snr, int noise_enable);
int decode_and_analyse(int numb_of_bits, float snr, int iter_num, int blocks, int guard_size, int trellis_term_enable, int decoder_kind, int guarding_type);
void plotting_2(int size_2, int size_3, int *guard_size_2, int *guard_size_3, int decoder_type, int iter);
float snr;
float snr_start = 0.0;
float snr_end = 2.005;
float resolution = 0.1;
int i,k, error_count, numb_snrs = 0;
int size_2 = 4;
int size_3 = 3;
int guard_size_2[4] = {3, 5, 8, 10};
int guard_size_3[3] = {3, 5, 8};
int *errordata;
FILE *fp, *gp;
if(only_plotting != 1)
{
	fp = fopen("ber_vs_snr_guard.dat", "w");
	gp = fopen("ber_vs_snr_guard_linear.dat", "w");
	permuter_bits();
	for(snr = snr_start; snr<=snr_end; snr += resolution)
	{
		encoder_and_noise(numb_of_bits, 0, snr, 1);
		fprintf(fp,"%f\t", snr);
		numb_snrs++;
		error_count = decode_and_analyse(numb_of_bits, snr, iter, NO_OF_BLOCKS, 0, trellis_term_enable, decoder_type, 1);
		fprintf(fp, "%d\t", error_count);
		fprintf(gp, "%d\n", error_count);
		for(i = 0; i<size_2; i++)
		{
			error_count = decode_and_analyse(numb_of_bits, snr, iter, NO_OF_BLOCKS, guard_size_2[i], trellis_term_enable, decoder_type, 2);
			fprintf(fp, "%d\t", error_count);
			fprintf(gp, "%d\n", error_count);
		
		}
		for(i = 0; i<size_3; i++)
		{
			error_count = decode_and_analyse(numb_of_bits, snr, iter, NO_OF_BLOCKS, guard_size_3[i], trellis_term_enable, decoder_type, 3);
			fprintf(fp, "%d\t", error_count);
			fprintf(gp, "%d\n", error_count);
		}
	
		error_count = decode_and_analyse(numb_of_bits, snr, iter, NO_OF_BLOCKS, 0, trellis_term_enable, decoder_type, 2);
		fprintf(fp, "%d\t", error_count);
		fprintf(gp, "%d\n", error_count);
		fprintf(fp, "\n");
		printf("Done for snr = %f\n", snr);
	}
	fclose(fp);
	fclose(gp);
	//Converting the data into a format that can be used by gnuplot
	gp = fopen("ber_vs_snr_guard_linear.dat", "r");
	errordata = (int*)malloc(numb_snrs*(2 + size_2 + size_3)*sizeof(int));
	for(i=0; i<numb_snrs*(2 + size_2 + size_3); i++)
	{
		fscanf(gp, "%d\n", &errordata[i]);
	}
	fclose(gp);
	fp = fopen("ber_vs_snr_guard_gnu.dat", "w");
	for(k=0; k<2 + size_2 + size_3; k++)
	{
		if(k == 0) fprintf(fp, "#####Previous value initialisation#####\n");
		else if(k>0 && k < 1+size_2) fprintf(fp, "#####Only training, guard size = %d#####\n", guard_size_2[k-1]);
		else if(k < 2 + size_2 + size_3 - 1) fprintf(fp, "#####Only training, guard size = %d#####\n", guard_size_2[k-1-size_2]);
		else fprintf(fp, "####No guarding at all, equal value initialisation#####\n");
		for(i=0; i<numb_snrs; i++)
		{
			fprintf(fp,"%f\t%d\t%g\n",snr_start+i*resolution,errordata[k+i*(2+size_2+size_3)], errordata[k+i*(2+size_2+size_3)]/(numb_of_bits*(1.0)));
		}
		fprintf(fp, "\n\n\n\n");
	}
	fclose(fp);
	free(errordata);
}
plotting_2(size_2, size_3, guard_size_2, guard_size_3, decoder_type, iter);
}
