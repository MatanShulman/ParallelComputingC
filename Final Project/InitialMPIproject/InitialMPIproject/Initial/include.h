struct NODE { int parent; float prob;  float func; }typedef node;
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>  
#include <cstdlib>
#include <ctime>
#include "iostream"
#include <time.h>
#define _CRT_SECURE_NO_DEPRECATE
#define observSize 30000
#define observRows 30000
#define StatesCols 1000
#define ABPATH "D:\\data\\AB.txt"
#define ObservationPath "D:\\data\\Observation.txt"
#define transtionPath "D:\\data\\Transition.txt"
#define transtionPathDone "D:\\data\\TransitionPath.txt"
