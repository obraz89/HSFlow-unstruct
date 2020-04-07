// HSFlow-unstruct.cpp : Defines the entry point for the console application.
//


//#include <sys/types.h>
//#include <sys/stat.h>
//#include <string.h>   // strcmp

#include "common_data.h"
#include "common_procs.h"

#include "flow_common.h"

//--------------------------------< Functions >---------------------------------

/**
* Make std::string with printf-style formatting
*
* @param fmt format string like "%s"
* @return
*/


TState G_State;    // solving state

// Computational domain composed of blocks with own grid and field
t_Domain G_Domain;
