void writePathToFile(int arr[], char fpath[]);
void memInitMatrix(node ***trellis);
void normalizstion(float **transtion);												//Optional
void matrixPrint(node **matrix, int cols, int rows);								//Optional
void printOobservtion(float *observation, int size);								//Optional
void printTranstion(float **transtion);												//Optional
void printEmission(float **emission);												//Optional
void initTrellis(node ***trellis, int rows, float * observ, float ** ab);
void memInitTranstion(float ***trans);
void memInitObservation(float **obsrv);
void memInitEmission(float ***ab);
void calculate(node**nodeArr, int end, float **transtion, int myid);
void findPath(node **trellis, int observrows, int start, int * findPath, int end);
void masterPrintPath(float * finalPath);
bool loadArrayFromFile(float arr[], int size, char fpath[]);
bool loadMatrixFromFile(float *mat[], int rows, int cols, char fpath[]);

