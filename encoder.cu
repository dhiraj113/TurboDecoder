#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hashdefined.h"

void encoder_and_noise(int numb_of_bits, int trellis_termination, float snr, int noise_enable)
//Output of this is a file by name "bits_to_be_decoded.dat"
{
void bit_generator(int numb_of_bits, int trellis_termination);
void just_encoder(int numb_of_bits);
void awgnchannel(int numb_of_bits, float snr, int noise_enable);
bit_generator(numb_of_bits, trellis_termination);
just_encoder(numb_of_bits);
awgnchannel(numb_of_bits, snr, noise_enable);
}

void bit_generator(int numb_of_bits, int trellis_termination)
//Generates random bits and stores them in the file named "input_bits.dat"
//Takes as input number of bits to generate
//Expects '0' or '1' as the trellis termination bit
//'0' --> No trellis termination
//'1' --> Trellis termination is present
//Encoding is necessary to do the trellis termination
{
int next_state(int current_state, int input);
int i, k, numb = numb_of_bits/DATASIZE;
int input_bit;
int currentstate = 0;
int trellis_term_table[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,1,1},{1,1,1},{1,0,1},{0,0,1}};
numb_of_bits = ((numb_of_bits+DATASIZE-1)/DATASIZE)*DATASIZE;
if(trellis_termination == 0 || trellis_termination == 1)
{
	FILE *fp;
	fp = fopen("input_bits.dat", "w");
	if(trellis_termination == 0)
	{
		srand((unsigned int)time((time_t *)NULL)); //modifying the seed for rand based on time
		//srand(time(NULL));
		for(i = 1; i <= numb_of_bits; i++)
		{
			fprintf(fp, "%d\n", rand() % 2);
		}
	}
	else if(trellis_termination == 1)
	{
		for(k=1; k<= numb; k++)
		{
			srand((unsigned int)time((time_t *)NULL));
			for(i=1; i<=DATASIZE-3; i++)
			{
				input_bit = rand() % 2;
				currentstate = next_state(currentstate, input_bit);
				fprintf(fp, "%d\n", input_bit);
			}
			fprintf(fp, "%d\n", trellis_term_table[currentstate][0]);
			fprintf(fp, "%d\n", trellis_term_table[currentstate][1]);
			fprintf(fp, "%d\n", trellis_term_table[currentstate][2]);
		}
	}
	fclose(fp);
}
else
{
	printf("Improper input to decide upon trellis termination\n");
	printf("Exiting abruptly\n");
	exit(1);
}	
}


void just_encoder(int numb_of_bits)
//encodes the bits stored in the file "input_bits.dat"
//Stores the encoded bits in a new file "encoded_bits.dat"
{
int interleaver(int index);
int next_state(int current_state, int input);
int output(int current_state, int input);
FILE *fp, *gp;
int inputbits[DATASIZE], deinter_bits[DATASIZE], encbits_1[DATASIZE], encbits_2[DATASIZE];
int currentstate;
int i,k, numb;
fp = fopen("input_bits.dat","r");
gp = fopen("encoded_bits.dat","w");
numb_of_bits = ((numb_of_bits+DATASIZE-1)/DATASIZE)*DATASIZE;
numb = numb_of_bits/DATASIZE;
for(k = 0; k<numb; k++)
{
	for(i = 0; i < DATASIZE; i++)
	{
		fscanf(fp, "%d\n", &inputbits[i]);
		deinter_bits[interleaver(i)] = inputbits[i];
	}
	//*****************Encoder 1********************
	currentstate = 0;
	for(i = 0; i<DATASIZE; i++)
	{
		encbits_1[i] = output(currentstate, inputbits[i]);
		currentstate = next_state(currentstate, inputbits[i]);
	}
	//*****************Encoder 1********************
	//*****************Encoder 2********************
	currentstate = 0;
	for( i = 0; i<DATASIZE; i++)
	{
		encbits_2[i] = output(currentstate, deinter_bits[i]);
		currentstate = next_state(currentstate, deinter_bits[i]);
	}
	//*****************Encoder 2********************
	for(i = 0; i < DATASIZE; i++)
	{
		fprintf(gp, "%d\t%d\t%d\n", inputbits[i]==0 ? -1 : 1, encbits_1[i], encbits_2[i]);
	}
}
fclose(fp);
fclose(gp);
}



void awgnchannel(int numb_of_bits, float snr, int noise_enable)
//Takes as input the encoded bits file by name encoded_bits.dat
//Adds noise to the encoded bits of that file
//Takes as input the snr of the noise to be added
{
float gasdev(long *idum);
int bit, enc1_bit, enc2_bit, i;
float outbit_x, outbit_y1, outbit_y2;
FILE *fp, *gp;
fp = fopen("encoded_bits.dat", "r");
gp = fopen("channel_out.dat", "w");
long int seedvalue = ((unsigned int)time((time_t *)NULL));//seed needs to be negative for ran1() which gasdev() calls
seedvalue = -seedvalue;
float sigma, rate = 1.0/3.0;
sigma = pow(10,-snr/20)/sqrt(2*rate); //The SNR value input is in db
numb_of_bits = ((numb_of_bits+DATASIZE-1)/DATASIZE)*DATASIZE;
if(noise_enable == 1)
{
	for(i=0; i<numb_of_bits; i++)
	{
		//***********************with noise****************************
		fscanf(fp, "%d\t%d\t%d\n", &bit, &enc1_bit, &enc2_bit);
		outbit_x = bit + sigma*gasdev(&seedvalue);
		outbit_y1 = enc1_bit + sigma*gasdev(&seedvalue);
		outbit_y2 = enc2_bit + sigma*gasdev(&seedvalue);
		fprintf(gp, "%f\t%f\t%f\n", outbit_x, outbit_y1, outbit_y2);
		//***********************with noise****************************/
	}
}
else if(noise_enable == 0)
{
	for(i=0; i<numb_of_bits; i++)
	{
		//***********************without noise****************************
		fscanf(fp, "%d\t%d\t%d\n", &bit, &enc1_bit, &enc2_bit);
		outbit_x = bit;
		outbit_y1 = enc1_bit;
		outbit_y2 = enc2_bit;
		fprintf(gp, "%f\t%f\t%f\n", outbit_x, outbit_y1, outbit_y2);
		//***********************without noise****************************/
	}
}
fclose(fp);
fclose(gp);
}

void permuter_bits()
{
int interleaver(int n);
extern int inv_permutation_bits[DATASIZE];
int k;
for(k=0; k<DATASIZE; k++)
{
	inv_permutation_bits[(((1+6*k)%DATASIZE)*k)%DATASIZE] = k;
}
/**********************Deinterleaver check block****************
FILE *fp;
fp = fopen("deinter.bits.dat", "w");
for(k=0; k<DATASIZE; k++)
{

	fprintf(fp,"%d\t%d\t%d\n",k,interleaver(k), inv_permutation_bits[k]);
}
**********************Deinterleaver check block****************/
}

int interleaver(int n)
{
return (((1+6*n)%DATASIZE)*n)%DATASIZE;
}


int deinterleaver(int index)
{
//performs de-interleaving using the already created permutation_bits array
extern int inv_permutation_bits[DATASIZE];
return inv_permutation_bits[index];
}
