#ifndef __LOG_H__
#define __LOG_H__

#include <stdlib.h>


#define LOG(_msg, ...)  	fprintf(stderr, "log: " _msg "\n", ##__VA_ARGS__);

#define DIE(_msg, ...)		{ fprintf(stderr, "die: " _msg "\n", ##__VA_ARGS__); exit(1); }

#ifdef _DEBUG_
#define DEBUG(_msg, ...)	fprintf(stderr, "debug: " _msg "\n", ##__VA_ARGS__);
#else
#define DEBUG(_msg, ...)	do {} while (0) ;
#endif


#endif 
