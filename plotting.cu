#include <stdio.h>
#include <stdlib.h>
#include "hashdefined.h"
void ber_plotting_1(int size, int *iter_array, int decoder_type, int guarding_type, int guard_size)
{
//Calling gnuplot and setting up various things in it
int i;
FILE *fp;
printf("Calling gnuplot\n");
fp = fopen("ber_gnucommands_1.txt","w");
fprintf(fp,"set term postscript eps enhanced color\n");
fprintf(fp,"set output \"BER_vs_SNR_plot_iter.ps\"\n");
fprintf(fp,"set xlabel \"SNR(in dB)\"\n");
fprintf(fp,"set ylabel \"BER(logscale)\"\n");
fprintf(fp,"set title \"BER vs SNR plot for different iterations\"\n");
switch(guarding_type)
{
	case 1 : {
			fprintf(fp,"set title \"Guarding type - only Previous value Initialisation\"\n");
		 	break;
		 }
	case 2 : {
			fprintf(fp,"set title \"Guarding type - only training window\"\tGuard size = %d\n", guard_size);
			break;
		 }
	case 3 : {
			fprintf(fp,"set title \"Guarding type - Previous value initialisation + training window, Guard size = %d\"\n", guard_size);
			break;
		 }
}
fprintf(fp,"set logscale y\n");
fprintf(fp,"plot");
for(i = 0; i<size-1; i++)
{
	fprintf(fp,"\"ber_vs_snr_iters_gnu.dat\" index %d:%d using 1:3 title \"Iter %d\" with lines,\\\n",i,i,iter_array[i]);
}
fprintf(fp,"\"ber_vs_snr_iters_gnu.dat\" index %d:%d using 1:3 title \"Iter %d\" with lines\n", size-1, size-1, iter_array[size-1]);
fclose(fp);
system("gnuplot ber_gnucommands_1.txt"); //using the script file gnucommands to plot
printf("Done plotting, go to the folder to see the plot\n");
}

void ber_plotting_3(int size, int *blocks_array, int decoder_type, int iter, int guarding_type, int guard_size)
{
//Calling gnuplot and setting up various things in it
int i;
FILE *fp;
printf("Calling gnuplot\n");
fp = fopen("ber_gnucommands_3.txt","w");
fprintf(fp,"set term postscript eps enhanced color\n");
fprintf(fp,"set output \"BER_vs_SNR_plot_block.ps\"\n");
fprintf(fp,"set xlabel \"SNR(in dB)\"\n");
fprintf(fp,"set ylabel \"BER(logscale)\"\n");
fprintf(fp,"set title \"BER vs SNR plot for different no. of blocks\"\n");
switch(guarding_type)
{
	case 1 : {
			fprintf(fp,"set title \"Guarding type - only Previous value Initialisation\"\n");
		 	break;
		 }
	case 2 : {
			fprintf(fp,"set title \"Guarding type - Only training window, Guard size = %d\"\n", guard_size);
			break;
		 }
	case 3 : {
			fprintf(fp,"set title \"Guarding type - Previous value initialisation + training window, Guard size = %d\"\n", guard_size);
			break;
		 }
}
fprintf(fp,"set logscale y\n");
fprintf(fp,"plot");
for(i = 0; i<size-1; i++)
{
	fprintf(fp,"\"ber_vs_snr_blocks_gnu.dat\" index %d:%d using 1:3 title \"# of blocks %d\" with lines,\\\n",i,i,blocks_array[i]);
}
fprintf(fp,"\"ber_vs_snr_blocks_gnu.dat\" index %d:%d using 1:3 title \"# of blocks %d\" with lines\n", size-1, size-1, blocks_array[size-1]);
fclose(fp);
system("gnuplot ber_gnucommands_3.txt"); //using the script file gnucommands to plot
printf("Done plotting, go to the folder to see the plot\n");
}

void ber_plotting_2(int size_2, int size_3, int *guard_size_2, int *guard_size_3, int decoder_type, int iter)
{
//Calling gnuploat and setting up various things in it
int i;
FILE *fp;
printf("Calling gnuplot\n");
fp = fopen("ber_gnucommands_2.txt","w");
fprintf(fp,"set term postscript eps enhanced color\n");
fprintf(fp,"set output \"BER_vs_SNR_plot_guard.ps\"\n");
fprintf(fp,"set xlabel \"SNR(in dB)\"\n");
fprintf(fp,"set ylabel \"BER(logscale)\"\n");
fprintf(fp,"set title \"BER vs SNR plot for different types of guarding, No. of iterations = %d\"\n", iter);
fprintf(fp,"set logscale y\n");
fprintf(fp,"plot");
fprintf(fp,"\"ber_vs_snr_guard_gnu.dat\" index 0:0 using 1:3 title \"Prev val init\" with lines,\\\n");
for(i = 1; i<size_2+1; i++)
{
	fprintf(fp,"\"ber_vs_snr_guard_gnu.dat\" index %d:%d using 1:3 title \"Only guard, size=%d\" with lines,\\\n",i,i,guard_size_2[i-1]);
}
for(i = size_2+1; i<size_3+size_2+1; i++)
{
	fprintf(fp,"\"ber_vs_snr_guard_gnu.dat\" index %d:%d using 1:3 title \"Prev val init and guard, size=%d\" with lines,\\\n",i,i,guard_size_3[i-1-size_2]);
}
fprintf(fp,"\"ber_vs_snr_guard_gnu.dat\" index %d:%d using 1:3 title \"Without guarding at all\" with lines\n", size_3+size_2+1, size_3+size_2+1);
//fprintf(fp,"\"ber_vs_snr_guard_gnu.dat\" index %d:%d using 1:3 title \"Max-log map with no parallelisation\" with lines\n", size_3+size_2+2, size_3+size_2+2);
fclose(fp);
system("gnuplot ber_gnucommands_2.txt"); //using the script file gnucommands to plot
printf("Done plotting, go to the folder to see the plot\n");
}

void fer_plotting_1(int size, int *iter_array, int decoder_type, int guarding_type, int guard_size)
{
//Calling gnuplot and setting up various things in it
int i;
FILE *fp;
printf("Calling gnuplot\n");
fp = fopen("fer_gnucommands_1.txt","w");
fprintf(fp,"set term postscript eps enhanced color\n");
fprintf(fp,"set output \"FER_vs_SNR_plot_iter.ps\"\n");
fprintf(fp,"set xlabel \"SNR(in dB)\"\n");
fprintf(fp,"set ylabel \"FER(logscale)\"\n");
fprintf(fp,"set title \"FER vs SNR plot for different iterations\"\n");
switch(guarding_type)
{
	case 1 : {
			fprintf(fp,"set title \"Guarding type - only Previous value Initialisation\"\n");
		 	break;
		 }
	case 2 : {
			fprintf(fp,"set title \"Guarding type - only training window\"\tGuard size = %d\n", guard_size);
			break;
		 }
	case 3 : {
			fprintf(fp,"set title \"Guarding type - Previous value initialisation + training window, Guard size = %d\"\n", guard_size);
			break;
		 }
}
fprintf(fp,"set logscale y\n");
fprintf(fp,"plot");
for(i = 0; i<size-1; i++)
{
	fprintf(fp,"\"fer_vs_snr_iters_gnu.dat\" index %d:%d using 1:3 title \"Iter %d\" with lines,\\\n",i,i,iter_array[i]);
}
fprintf(fp,"\"fer_vs_snr_iters_gnu.dat\" index %d:%d using 1:3 title \"Iter %d\" with lines\n", size-1, size-1, iter_array[size-1]);
fclose(fp);
system("gnuplot fer_gnucommands_1.txt"); //using the script file gnucommands to plot
printf("Done plotting, go to the folder to see the plot\n");
}




void fer_plotting_3(int size, int *blocks_array, int decoder_type, int iter, int guarding_type, int guard_size)
{
//Calling gnuplot and setting up various things in it
int i;
FILE *fp;
printf("Calling gnuplot\n");
fp = fopen("fer_gnucommands_3.txt","w");
fprintf(fp,"set term postscript eps enhanced color\n");
fprintf(fp,"set output \"FER_vs_SNR_plot_block.ps\"\n");
fprintf(fp,"set xlabel \"SNR(in dB)\"\n");
fprintf(fp,"set ylabel \"FER(logscale)\"\n");
fprintf(fp,"set title \"FER vs SNR plot for different no. of blocks\"\n");
switch(guarding_type)
{
	case 1 : {
			fprintf(fp,"set title \"Guarding type - only Previous value Initialisation\"\n");
		 	break;
		 }
	case 2 : {
			fprintf(fp,"set title \"Guarding type - Only training window, Guard size = %d\"\n", guard_size);
			break;
		 }
	case 3 : {
			fprintf(fp,"set title \"Guarding type - Previous value initialisation + training window, Guard size = %d\"\n", guard_size);
			break;
		 }
}
fprintf(fp,"set logscale y\n");
fprintf(fp,"plot");
for(i = 0; i<size-1; i++)
{
	fprintf(fp,"\"fer_vs_snr_blocks_gnu.dat\" index %d:%d using 1:3 title \"# of blocks %d\" with lines,\\\n",i,i,blocks_array[i]);
}
fprintf(fp,"\"fer_vs_snr_blocks_gnu.dat\" index %d:%d using 1:3 title \"# of blocks %d\" with lines\n", size-1, size-1, blocks_array[size-1]);
fclose(fp);
system("gnuplot fer_gnucommands_3.txt"); //using the script file gnucommands to plot
printf("Done plotting, go to the folder to see the plot\n");
}

void fer_plotting_2(int size_2, int size_3, int *guard_size_2, int *guard_size_3, int decoder_type, int iter)
{
//Calling gnuploat and setting up various things in it
int i;
FILE *fp;
printf("Calling gnuplot\n");
fp = fopen("fer_gnucommands_2.txt","w");
fprintf(fp,"set term postscript eps enhanced color\n");
fprintf(fp,"set output \"FER_vs_SNR_plot_guard.ps\"\n");
fprintf(fp,"set xlabel \"SNR(in dB)\"\n");
fprintf(fp,"set ylabel \"FER(logscale)\"\n");
fprintf(fp,"set title \"FER vs SNR plot for different types of guarding, No. of iterations = %d\"\n", iter);
fprintf(fp,"set logscale y\n");
fprintf(fp,"plot");
fprintf(fp,"\"fer_vs_snr_guard_gnu.dat\" index 0:0 using 1:3 title \"Prev val init\" with lines,\\\n");
for(i = 1; i<size_2+1; i++)
{
	fprintf(fp,"\"fer_vs_snr_guard_gnu.dat\" index %d:%d using 1:3 title \"Only guard, size=%d\" with lines,\\\n",i,i,guard_size_2[i-1]);
}
for(i = size_2+1; i<size_3+size_2+1; i++)
{
	fprintf(fp,"\"fer_vs_snr_guard_gnu.dat\" index %d:%d using 1:3 title \"Prev val init and guard, size=%d\" with lines,\\\n",i,i,guard_size_3[i-1-size_2]);
}
fprintf(fp,"\"fer_vs_snr_guard_gnu.dat\" index %d:%d using 1:3 title \"Without guarding at all\" with lines\n", size_3+size_2+1, size_3+size_2+1);
//fprintf(fp,"\"fer_vs_snr_guard_gnu.dat\" index %d:%d using 1:3 title \"Max-log map with no parallelisation\" with lines\n", size_3+size_2+2, size_3+size_2+2);
fclose(fp);
system("gnuplot fer_gnucommands_2.txt"); //using the script file gnucommands to plot
printf("Done plotting, go to the folder to see the plot\n");
}


