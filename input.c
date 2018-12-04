#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"

inputStruct* commandLineParser( int argc, char* argv[] )
{
	int i, length;
	char *in="-i", *c="-c", *o="-o", *d="-d", *frechet="Frechet", *DTW="DTW", *inputFilePath, *configurationFilePath, *outputFilePath, *metric;
	inputStruct *input;

	input = malloc( sizeof(inputStruct) );

	for( i=1; i<argc; i+=2 )
	{
		if( strcmp(argv[i], in) == 0 )
		{
			length = strlen(argv[i+1]);
			input -> inputFilePath = malloc( (length+1) * sizeof(char) );
			strcpy(input -> inputFilePath, argv[i+1]);
		}
		else if( strcmp(argv[i], c) == 0 )
		{
			length = strlen(argv[i+1]);
			input -> configurationFilePath = malloc( (length+1) * sizeof(char) );
			strcpy(input -> configurationFilePath, argv[i+1]);
		}
		else if( strcmp(argv[i], o) == 0 )
		{
			length = strlen(argv[i+1]);
			input -> outputFilePath = malloc( (length+1) * sizeof(char) );
			strcpy(input -> outputFilePath, argv[i+1]);
		}
		else if(  strcmp(argv[i], d) == 0 )
		{
			length = strlen(argv[i+1]);
			input -> metric = malloc( (length+1) * sizeof(char) );
			strcpy(input -> metric, argv[i+1]);
		}
		else
		{
			printf("Error: Command line arguments not valid\n");
			exit(EXIT_FAILURE);
		}
	}
	
	return input;
}

inputStruct_user* userParser()
{
	char userLine[SIZE];
	char *split;
	inputStruct_user *input;

	input = malloc( sizeof(inputStruct_user) );

	printf("Choose initialization { 1:random, 2:kmeans++ }\n");
	printf("> ");
	fgets(userLine, SIZE, stdin);

	split = strtok(userLine, "\r\n");
	if (strcmp(split, "1") == 0)
		input -> initializationMethod = 1;
	else if (strcmp(split, "2") == 0)
		input -> initializationMethod = 2;
	else
		printf("Error: only {1, 2} are acceptable values\n");

	memset(userLine, 0, SIZE);

	printf("Choose assignment { 1:Lloyd's, 2:Range search }\n");
	printf("> ");
	fgets(userLine, SIZE, stdin);

	split = strtok(userLine, "\r\n");
	if (strcmp(split, "1") == 0)
		input -> assignmentMethod = 1;
	else if (strcmp(split, "2") == 0)
		input -> assignmentMethod = 2;
	else
		printf("Error: only {1, 2} are acceptable values\n");

	return input;
}