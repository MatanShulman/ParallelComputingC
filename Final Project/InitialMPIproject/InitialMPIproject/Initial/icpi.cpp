/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
*  (C) 2001 by Argonne National Laboratory.
*      See COPYRIGHT in top-level directory.
*/

/* This is an interactive version of cpi */
#include "include.h"
#include "Function.h"

using namespace std;

int main(int argc, char *argv[])
{
	node	**trellis;
	node	** nodeMat;
	float	number, sum;
	float	*observation, *slaveObserv;
	float	**emission, **transtion;
	double t1,t2;
	int		end = 0, process = 0, start = 0, amount = 0, simple = 1, flag = 1;
	int		i, j, z, mallocNode, cols, rows, newPath, counter, parent, namelen, numprocs, myid;
	int		* finalPath, *temp, *zeroVector;
	
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	MPI_Get_processor_name(processor_name, &namelen);

	MPI_Status status;


	//memory allocation function to ab matrix and transtion matrix
	memInitEmission(&emission);
	memInitTranstion(&transtion);
	
	//loading data from file to ab matrix and transtion matrix
	if (myid == 0){
		loadMatrixFromFile(emission, StatesCols, 2, ABPATH);
		loadMatrixFromFile(transtion, StatesCols, StatesCols, transtionPath);
		}


	//broadcast ab matrix and transition matrix to all processes
	for (i = 0; i < StatesCols; i++)
		MPI_Bcast(emission[i], 2, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
	for (i = 0; i < StatesCols; i++)
		MPI_Bcast(transtion[i], StatesCols, MPI_FLOAT, 0, MPI_COMM_WORLD);
	

	
	for (i = 0; i < StatesCols; i++)
		#pragma omp parallel for private(j)										//!!!!!!!!!!!!!!!!!!
			for (j = 0; j < StatesCols; j++)
				transtion[i][j] = log(transtion[i][j]);
	
	
	//array for master to print in the end after all procsess are done
	finalPath = (int*)malloc(sizeof(int)*observSize);	
	
	//temporary array for all slave to work on and  copy to matser final path array
	temp = (int*)malloc(sizeof(int)*observSize);			

	// initial temp array for validtion in the "print path" function 
	#pragma omp parallel for private(i)	
		for (i = 0; i < observSize; i++)
				temp[i] = observRows + 1;
	
	if (myid == 0){
		t1 = MPI_Wtime();
		//memory allocation for ab matrix
		memInitObservation(&observation);
		
		//loading data for ab matrix
		loadArrayFromFile(observation, observRows, ObservationPath);
		
		//vector for observtion zero's;
		zeroVector = (int*)malloc(sizeof(int)* 1000);			
		counter = 1;
		zeroVector[0] = 0;
		for (i = 1; i < observRows; i++){
			if (observation[i] == 0){
				zeroVector[counter] = i;
				counter++;
			}
		}

		cout << " THE NUMBER OF ZERO'S IN VECTOR \"ZERO ARE\" ARE :  " << counter << "\n";
		fflush(stdout);

		zeroVector[counter] = observRows - 1;

		for (int i = 1; i < numprocs; i++)
			MPI_Send(&counter, 1, MPI_INT, i, 0, MPI_COMM_WORLD);			//Send  number of zero's  to all procesess
		

		if (numprocs - 1 >= counter){
			for (int i = 1; i < numprocs; i++)
				MPI_Send(&simple, 1, MPI_INT, i, 0, MPI_COMM_WORLD);		//Send mode number to all  procesess or complex or simple mode
		}
		else{
			simple = 0;
			for (int i = 1; i < numprocs; i++){
				MPI_Send(&simple, 1, MPI_INT, i, 0, MPI_COMM_WORLD);		
			}
		}
	
		if (simple){									
			cout << " MASTER INSIDE SINPLE MODE  " << "\n";
			fflush(stdout);
			//Send  Data to  all procesess that enter simple mode 
			for (i = 1; i < counter + 1; i++){
				MPI_Send(&zeroVector[i - 1], 1, MPI_INT, i, 0, MPI_COMM_WORLD);								   //Send start signal for allocation matrix
				MPI_Send(&zeroVector[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);									   //Send end signal for allocation matrix
				MPI_Send(&amount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);										   //amount-> point to next part of observtion to send to slave
				MPI_Send(observation + amount, (zeroVector[i]) - amount, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
				amount = zeroVector[i] + 1;
			}
			//All procesess finish
			for (i = 1; i < counter + 1; i++){
				//Receive  the temp array from all procesess to copy it by order to final path
				MPI_Recv(&process, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);						
				MPI_Recv(&start, 1, MPI_INT, process, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&end, 1, MPI_INT, process, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(temp, observRows, MPI_FLOAT, process, 0, MPI_COMM_WORLD, &status);
				if (start == 0){
					for (z = start; z <= end; z++){
						finalPath[z] = temp[z];
					}
				}
				else{
					for (z = start + 1; z <= end; z++){
						finalPath[z] = temp[z];
					}
				}
			}
			writePathToFile(finalPath, transtionPathDone);

			cout << "\n MASTER IS FINISH    " << "\n";
			fflush(stdout);
		}

		else{
			cout << " MASTER INSIDE COMPLEX MODE   " << "\n";
			fflush(stdout);
			for (i = 1; i < counter + 1; i++){
				if (i + 1 < counter + 1){																//Check for last itiration
					MPI_Recv(&process, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);			
					flag = 1;																			//Signal for more work left
					MPI_Send(&flag, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
				}
				else{																					//Last itiration !!!
					flag = 0;
					for (z = 1; z < numprocs - 1; z++){
						MPI_Recv(&process, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);		//Waitng to oene last slave to give work.
						MPI_Send(&flag, 1, MPI_INT, process, 0, MPI_COMM_WORLD);						//Signal him to Continue
					}
					flag = 2;
					MPI_Recv(&process, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);			//Free all slave from work
					MPI_Send(&flag, 1, MPI_INT, process, 0, MPI_COMM_WORLD);							//Signal  to slave and master that work is done 
				}

				//Send data to slave's in complex mode
				MPI_Send(&zeroVector[i - 1], 1, MPI_INT, process, 0, MPI_COMM_WORLD);					//send start signal
				MPI_Send(&zeroVector[i], 1, MPI_INT, process, 0, MPI_COMM_WORLD);						//send end signal
				MPI_Send(&amount, 1, MPI_INT, process, 0, MPI_COMM_WORLD);								//amount-> point to next part of observtion to send to slave
				MPI_Send(observation + amount, (zeroVector[i]) - amount, MPI_FLOAT, process, 0, MPI_COMM_WORLD);
				amount = zeroVector[i] + 1;
				
				if (flag == 2){
					cout << " MASTER IS FINISH IN COMPLEX MODE    " << "\n";
					fflush(stdout);
					break;
				}
			}
		}
	}
	else{
		MPI_Recv(&counter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);						//Receive number for simple mode to restrict invalid slaves
		MPI_Recv(&simple, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);						//Mode selection to all slave's
		if (simple){
			cout << "simple mode " << myid << "\n";
			fflush(stdout);
			if (myid <= counter){
				MPI_Recv(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);				//start signal for rows
				MPI_Recv(&mallocNode, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);			//End signal for rows
				MPI_Recv(&amount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);				//Receive number for  observtion memory alocation 
				observation = (float*)malloc(sizeof(float)*(mallocNode - amount));
				end = mallocNode - amount;													//Final end signal for rows
				MPI_Recv(observation, end, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);		
								
				//malloc and initial for temp matrix
				nodeMat = (node**)malloc(sizeof(node*)*((mallocNode - amount) + 1));
				for (i = 0; i < (mallocNode - amount) + 1; i++){
					*(nodeMat + i) = (node*)malloc(sizeof(node)* StatesCols);
				}
				initTrellis(&nodeMat, end + 1, observation, emission);
								
				calculate(nodeMat, end + 1, transtion, myid);		
				//matrixPrint(nodeMat, StatesCols, end + 1);								//Optional
				findPath(nodeMat, end + 1, start, finalPath, mallocNode);					//Return an array of path for nodemat matrix
				
				//Return the final path back to master
				MPI_Send(&myid, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);							
				MPI_Send(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(&mallocNode, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(finalPath, observRows, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
				for (int i = 0; i < end + 1; i++)
					free(nodeMat[i]);
				free(nodeMat);
				free(observation);
			}
		}
		else{

			flag = 1;
			//Complex mode 
			while (flag){
				MPI_Send(&myid, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
				if (flag == 0)																			//All process stop working exept one for last itiration 
					break;
				MPI_Recv(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);							//start signal
				MPI_Recv(&mallocNode, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);						//End signal for rows
				MPI_Recv(&amount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);							//Receive number for  observtion memory alocation
				observation = (float*)malloc(sizeof(float)*(mallocNode - amount));
				end = mallocNode - amount;																//Final end signal for rows
				MPI_Recv(observation, end, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
				
				//malloc and initial for temp matrix
				nodeMat = (node**)malloc(sizeof(node*)*((mallocNode - amount) + 1));
				for (i = 0; i < (mallocNode - amount) + 1; i++){
					*(nodeMat + i) = (node*)malloc(sizeof(node)* StatesCols);
				}
				initTrellis(&nodeMat, end + 1, observation, emission);
				
				calculate(nodeMat, end + 1, transtion, myid);
				//matrixPrint(nodeMat, StatesCols, end + 1);											//Optional
				findPath(nodeMat, end + 1, start, temp, mallocNode);
				for (int i = 0; i < end + 1; i++)	
					free(nodeMat[i]);
				free(nodeMat);
				free(observation);
				if (flag == 2)
					break;
				}
		}
}
	if (!simple){																				//Master printing to  file final path
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 0){
		for (i = 1; i < numprocs; i++){
			MPI_Recv(temp, observRows, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);	//Receive all the path from all procesess
			for (z = 0; z < observRows; z++){
				if (temp[z] != observRows + 1){													//Avoid override of other procesess path
					finalPath[z] = temp[z];
				}
			}
		}
	writePathToFile(finalPath, transtionPathDone);
	}
	else {
		MPI_Send(temp, observRows, MPI_INT, 0, 0, MPI_COMM_WORLD);								//All the sla
	}
}
	
	if (myid == 0) {
		t2 = MPI_Wtime();
		printf("\nthe time of process number %d is : %.2f\n", myid, t2 - t1);
	}
	free(emission);
	free(transtion);
	MPI_Finalize();
	return 0;


}

// Allocatoin  for slave's matrix
void memInitMatrix(node ***trellis){
	int i;
	//trellis for transtion matrix
	*trellis = (node**)malloc(sizeof(node*)*observRows);
	for (i = 0; i < observRows; i++){
		*(*trellis + i) = (node*)malloc(sizeof(node)* StatesCols);
	}
}
//allocatoin  for slaves Observation vector
void memInitObservation(float **obsrv){
	*obsrv = (float*)malloc(sizeof(float)*observRows);
}
//allocatoin  for slaves AB vector
void memInitEmission(float ***ab){
	int i;
	*ab = (float**)malloc(sizeof(float*)*StatesCols);
	for (i = 0; i < StatesCols; i++){
		*((*ab) + i) = (float*)malloc(sizeof(float)* 2);
	}

}
//allocatoin  for slaves transtion matrix
void memInitTranstion(float ***trans){
	int i = 0;
	*trans = (float**)malloc(sizeof(float*)*StatesCols);
	for (i = 0; i < StatesCols; i++){
		*((*trans) + i) = (float*)malloc(sizeof(float)*StatesCols);
	}

}
// algorhtem implantation
void calculate(node**nodeArr, int end, float **transtion,int myid){
	int i, j, z, tid;
	float number;
	float observNum;
	cout << "slave  number [" << myid <<"] IS START INSIDE CLAULCATE , end =  " <<end<< "\n";
	fflush(stdout);
	for (i = 0; i < end - 1; i++){
		#pragma omp parallel for private(number,j,z)	
			for (z = 0; z < StatesCols; z++){
				for (j = 0; j < StatesCols; j++){
					
					number = (nodeArr[i][j].prob + transtion[j][z] + nodeArr[i][j].func);
						if (number > nodeArr[i + 1][z].prob || nodeArr[i + 1][z].prob == 0){
						nodeArr[i + 1][z].prob = number;
						nodeArr[i + 1][z].parent = j;
					}
				}
	     	 }
    	 }
     }
//(Optional) Matrix print
void matrixPrint(node **matrix,int cols,int rows){
	int i, j;
	printf("original matrix \n");
	printf("---------------\n");
	for (i = 0; i < rows; i++){
		for (j = 0; j < cols; j++){
			printf("[%.20f, %d]  ", matrix[i][j].prob, matrix[i][j].parent);
		}
		printf("\n\n\n");

	}


}
//(Optional) Print observtion vector
void printOobservtion(float *observation,int size){
	int i;

	printf("observation array \n");
	printf("----------------- \n");
	for (i = 0; i < size; i++){
		printf("%.2f  ", observation[i]);
	}
	printf("\n\n\n");
}
//(Optional) Normalizstion for matrix vector
void normalizstion(float **transtion){
	int i, j, sum;
	for (i = 0; i < StatesCols; i++){
		sum = 0;
		for (j = 0; j < StatesCols; j++){
			sum += transtion[i][j];
		}

		for (j = 0; j < StatesCols; j++){
			transtion[i][j] = transtion[i][j] / sum;
		}


	}
}
//(Optional) Print transtion matrix
void printTranstion(float **transtion){
	int i, j;
	printf("transtion array \n");
	printf("--------------- \n");
	for (i = 0; i < StatesCols; i++){
		for (j = 0; j < StatesCols; j++)
			printf("%.2f  ", transtion[i][j]);
		printf("\n");
	}
	printf("\n\n");
}
//(Optional) Print AB matrix
void printEmission(float **emission){
	int i, j;
	printf("emission array \n");
	printf("--------------- \n");
	for (i = 0; i < StatesCols; i++){
		for (j = 0; j < 2; j++)
			printf("%.2f  ", emission[i][j]);
		printf("\n");
	}

	printf("\n\n");
}
//initial slave's matrix in simple or complex mode
void initTrellis(node ***trellis,int rows,float * observ, float ** ab){
	int i, j;
	float ob;

	for (i = 0; i < rows; i++){

		ob = observ[i];
	#pragma omp parallel for  private(j)	
		for (j = 0; j < StatesCols; j++){
			if (i == 0){
				(*trellis)[i][j].prob = 0; //log(1)
				(*trellis)[i][j].parent = -1;
				(*trellis)[i][j].func = log(ab[j][0] * exp((float)-1 * (pow(ob- ab[j][1], 2))));
			}

			else{
				(*trellis)[i][j].prob = 0;
				(*trellis)[i][j].parent = -1;
				(*trellis)[i][j].func = log(ab[j][0] * exp((float)-1 * (pow(ob - ab[j][1], 2))));
			}
		}
	}
}
//Find a path of the given matrix.
void findPath(node **trellis,int observrows,int start,int * findPath,int end){
	//cout << "FIND PATH START TO WORK " << "\n";
	//fflush(stdout);
	int  parent, rows, cols,i,j;
	float number = trellis[observrows - 1][0].prob;


	
	//printf("Finding Path: \n");
	//printf("------------ \n");

	//finding the biggset probbelaty in the last line 
	for (j = 0; j < StatesCols; j++){
		if (trellis[observrows - 1][j].prob >= number){
			number = trellis[observrows - 1][j].prob;
			parent = trellis[observrows - 1][j].parent;
			rows = observrows - 1;
			cols = j;

		}
	}

	// init the main array of the path that will go back to master.
	findPath[end] = -(cols + 1);
	end--;

	//cout << "START ANOTER PATH " << "\n";
	//fflush(stdout);
	//printf("Final Path : \n", cols);
	//printf(" state:%d ->", cols);


	//start go back from last through parent 
	for (i = observrows - 1; i > 0; i--){
		if (i != 1){
		//	printf(" state:%d  - > ", parent, i - 1);
			findPath[end] = parent;
			end--;
			parent = trellis[i - 1][parent].parent;
		//	printf("  parent is : %d  ", parent);
			
		}
		else {
		//	printf(" state:%d  ", parent, i - 1);
			parent = trellis[i][parent].parent;
			findPath[end] = parent;
		
		}
			
	}
		//cout << "number of rows is   " << observrows << "\n";
		//fflush(stdout);
		printf("\n\n\n");
	}
//Set random values to transtion matirx and ab matrix
void testAllValuesTE(float **trans, float *ab[]){
	float r3 ,a,b,c;
	int i, j;



	
	
	for (i = 0; i < StatesCols; i++){
		for (j = 0; j < StatesCols; j++){
			r3 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			if (i == 0 && j == 0){
				
				

				trans[i][j] = 0.5;
			}
			else if (i == StatesCols - 1 && j == StatesCols - 1){
			
				trans[i][j] = 0.5;
			}
			else
			
				trans[i][j] = r3;
		}
	}
	for (i = 0; i < StatesCols; i++){
		for (j=0; j < 2; j++){
			r3 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			ab[i][j] = r3;
		}
	}
	


}
//Set random values observion array
void testAllValuesOB(float obsrv[]){
	float r3;
	int i;
	
	for (i = 0; i < observRows; i++){
		r3 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		obsrv[i] = r3;
		if (i == 0 ){
			obsrv[i] = 0.5;
		}
		 if (i == observRows - 1){
			obsrv[i] = 0.5;
		}
		 if (i % 8 == 0){
	
		 obsrv[i] = 0;
		 }
		 

	}
}
//Loading data from file to Oobservtion vector
bool loadArrayFromFile(float arr[], int size, char fpath[])
{
	FILE* f = fopen(fpath, "r+");

	if (f == NULL)
	{
		printf("\nFailed opening the file..\n");
		return false;
	}

	for (int i = 0; i < size; i++)
	{
		fscanf(f, "%f", &arr[i]);
	//	printf("read number %d", i);
	}

	fclose(f);
	return true;
}
//Master path  printing  to file 
void writePathToFile(int arr[],char fpath[])
{
	int i = 0;
	FILE* f = fopen(fpath, "w");
	if (f == NULL){
		printf("\nFailed opening the file..\n");
		return;
	}
	else {
		fprintf(f, "NEW PATH \n-----------\n ");
		for (i = 0; i < observRows; i++){
			if (arr[i] < 0){
				fprintf(f, "\n\n\n ");
				fprintf(f, "NEW PATH \n-----------\n ", arr[i]);
				fprintf(f, "%d, ", -arr[i]);
				i++;
			}
			if (i == observRows)
				break;
			if (i != observRows - 1){
				
				fprintf(f, "%d, ", arr[i]);
			}
			else{
				fprintf(f, "%d", arr[i]);
			}
		}
	}
	fprintf(f, "\n\n\n Amount of total path parents's is : %d",i);
	fclose(f);
	return;
}
//load data from file to Ab
bool loadMatrixFromFile(float *mat[], int rows, int cols, char fpath[])
{
	FILE* f = fopen(fpath, "r+");

	if (f == NULL)
	{
		printf("\nFailed opening the file..\n");
		return false;
	}

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			fscanf(f, "%f", &mat[i][j]);
		}
	}

	fclose(f);
	return true;
}
//(Optional) Master print final path to consle
void masterPrintPath(float * finalPath){
	int i;
	for (i = 0; i < observRows; i++){
	printf(" %.2f ", finalPath[i]);
	}


	printf("\nStart to find a path :  \n-------------------\n");
	for (int i = 0; i < observRows; i++)
	{
		if (finalPath[i] >= 0){
			printf("%.2f  ", finalPath[i]);
		}

		if (finalPath[i] < 0){
			if (finalPath[i] + 1 == 0){
				printf("%.2f  ", (finalPath[i] + 1));
				printf("\n\na new path was found \n-------------------\n");
			}
			else{
				printf("%.2f  ", -(finalPath[i] + 1));
				printf("\n\na new path was found \n-------------------\n");
			}
		}
	}
}