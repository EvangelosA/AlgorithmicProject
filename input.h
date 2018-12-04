#ifndef INPUT__H__
#define INPUT__H__

#define SIZE 1000

typedef struct inputStruct
{
	char *inputFilePath;
	char *outputFilePath;
	char *configurationFilePath;
	char *metric;
}inputStruct;

typedef struct inputStruct_user
{
	int initializationMethod;
	int assignmentMethod;
}inputStruct_user;

inputStruct* commandLineParser( int , char** );
inputStruct_user* userParser();

#endif