#include <stdio.h>
//The following trellis functions are written specifically for
//g1(D)/g2(D) = (1 + D + D3 )/(1 + D2 + D3 )
//Need to be generalised
int next_state(int current_state, int input)
{
int cpu_beta_state_0[8] 	= 	{0, 4, 5, 1, 2, 6, 7, 3};
int cpu_beta_state_1[8] 	= 	{4, 0, 1, 5, 6, 2, 3, 7};
int temp=5;
if(input == 0 || input == -1)
{
	temp = cpu_beta_state_0[current_state];
}
else if(input == 1)
{
	temp = cpu_beta_state_1[current_state];
}
return temp;
}



int output(int current_state, int input)
{
int cpu_beta_encbit_0[8] 	= 	{-1, -1, 1, 1, 1, 1, -1, -1};
int cpu_beta_encbit_1[8] 	= 	{1, 1, -1, -1, -1, -1, 1, 1};
int temp=5;
if(input == 0 || input == -1)
{
	temp = cpu_beta_encbit_0[current_state];
}
else if(input == 1)
{
	temp = cpu_beta_encbit_1[current_state];
}
return temp;
}


