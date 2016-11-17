#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NB_CHANNELS  4
#define NB_SAMPLES   110

int main (void){

int i;

int *buffer = (double*)malloc(NB_CHANNELS*NB_SAMPLES*sizeof(double));


for (i=0;i<NB_CHANNELS*NB_SAMPLES;i++){
	buffer[i] = i;
	//printf(" %i ",buffer[i]);
	
	}
memcpy(&buffer[0],&buffer[(NB_CHANNELS*NB_SAMPLES/5)-1],(NB_CHANNELS*NB_SAMPLES));

for (i=0;i<NB_SAMPLES*NB_CHANNELS;i++){
	printf(" %i ",buffer[i]);
	
	
	}

}
