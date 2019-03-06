#include "pch.h"			//pre-compiled headers
#include "mpi.h"			//MPI library
#include "stdio.h"			//basic library
#include <iostream>			//for example for cout
#include <math.h>			//Math class for PI, fabs
#include <vector>			//Vector class
#include "stdlib.h"			//basic library
#include "matrix.h"			//Matrix class, inspired on c++ course
#include <chrono>			//Measure time in serialized program
#include <string>			//For read argv
#include <fstream>			//Save file

using namespace std;


void ThomasAlgorithm(vector<double> &mainVector, int sizeX, double deltaT, double deltaX, const double u)
{
	double alpha = 0, m = 0;
	vector<double> A(sizeX);
	vector<double> B(sizeX);
	vector<double> C(sizeX);

	alpha = u * (deltaT/deltaX);
	A[0] = 0.0;
	B[sizeX-1] = 1 + alpha;
	C[sizeX-1] = 0.0;

	for(int x = 0; x < sizeX - 1; x++)
	{
		A[x+1] = (-0.5) * alpha;
		B[x] = 1;
		C[x] = 0.5*alpha;
	}

	for(int x = 1; x < sizeX - 1; x++)
	{
		m = A[x]/B[x-1];
		B[x] = B[x] - m * C[x-1];
		mainVector[x] = mainVector[x] - m * mainVector[x-1];
	}

	mainVector[sizeX-1] = mainVector[sizeX-1]/B[sizeX-1];
	
	for(int x = sizeX - 2; x>=1; x--)
	{
		mainVector[x] = (mainVector[x] - C[x] * mainVector[x+1])/B[x];
	}		

}


int main(int argc, char* argv[])
{
	//Variables for MPI
	double t1, t2;		//variable for time performance
	int process_id;		//process_id: label for process, 
	int processes; 		//processes: number of processes from console
	MPI_Init(&argc, &argv);
	MPI_Comm communicator;
	communicator = MPI_COMM_WORLD;
	MPI_Comm_size(communicator, &processes);
	MPI_Comm_rank(communicator, &process_id);
	MPI_Status status;

	//Variables for Implicit
	string s = argv[1];
	double deltaT = stod(s);
	const double PI = atan(1) * 4; //Constant double value PI = atan(1) * 4
	const double u = 250;	//Constant double value from excercise 250 m/s
	double deltaX = 0.5;	//The same for each deltaT
	const double x = 400;	//Domain is xE[0;400] meters
	const double t = 0.5;	//t goes from 0.0 to 0.5sec
	int sizeX = 0; 
	int sizeT = 0;

	//compute sizeX (columns) and sizeT (rows), this is size of matrix
	sizeX = (x / deltaX);
	sizeT = (t / deltaT);
	
	//The final matrix with results from Implicit Scheme
	Matrix matrixGlobal(sizeT,sizeX);
	
	//inital condition in vector
	vector<double> mainVector(sizeX);

	if (process_id == 0)
	{
		
		//Initial Condition
		for (int x = 0; x < sizeX; x++)
		{
			if (x >= 0 && x <= 100)
			{
				mainVector[x] = 0;
			}
			if (x >= 100 && x <= 220)
			{
				mainVector[x] = 100 * (sin((PI)*((((double)x * 0.5) - 50) / 60)));
			}
			if (x >= 220 && x < sizeX)
			{
				mainVector[x] = 0;
			}
		}
		for(int i = 0; i<sizeX;i++)
				matrixGlobal[0][i] = mainVector[i];
		//Send vector to next process
		MPI_Send(&mainVector[0], sizeX, MPI_DOUBLE, process_id + 1, 0, communicator);
		MPI_Send(&matrixGlobal[0][0], sizeX*sizeT, MPI_DOUBLE, process_id + 1, 1, communicator);

	}

	
	
	t1 = MPI_Wtime();	//start measure time

	//Implicit Scheme FTBS
	if(process_id>0)
	{
		//P1 receives data from P0 (only once).
		if(process_id==1)
		{
			MPI_Recv(&mainVector[0], sizeX, MPI_DOUBLE, 0, 0, communicator, &status);	
			MPI_Recv(&matrixGlobal[0][0],sizeX*sizeT, MPI_DOUBLE, 0, 1, communicator, &status);
		}
		//Current process (different than 0 and 1) receives data from previous process.
		else
		{
			
			MPI_Recv(&matrixGlobal[0][0],sizeX*sizeT, MPI_DOUBLE, process_id-1, 2, communicator, &status);
			MPI_Recv(&mainVector[0], sizeX, MPI_DOUBLE, process_id-1, 3, communicator, &status);	
		}
		
		
		//How many iterations can each process do?
		int k = (sizeT-1) / (processes-1);
		int tmp = (process_id * (k))-(k)+1;
		
		if(process_id==processes-1)
		{
			k = (sizeT - 1) - tmp;
			
		}

		//Solve problem with Thomas Algorithm k times
		for(int i = 0; i < k; i++)
		{
			ThomasAlgorithm(mainVector, sizeX, deltaT, deltaX, u);
			for(int j = 0; j<sizeX;j++)
				matrixGlobal[tmp][j] = mainVector[j];
			tmp++;
		}
		
		//Current process(different than 0 and 1) send data to next process.
		if(process_id!=processes-1)
		{
			MPI_Send(&matrixGlobal[0][0], sizeX*sizeT, MPI_DOUBLE, process_id+1, 2, communicator);
			MPI_Send(&mainVector[0], sizeX, MPI_DOUBLE, process_id+1, 3, communicator);
		}
	}

	t2 = MPI_Wtime();	//end measure time
	cout << "\n\nElapsed time is "<< t2 - t1 <<" for process: "<<process_id;
	
	//Save Matrix after calculations	
	if (process_id == processes-1)
	{		

		//equation computing step for each size matrix 
		int step = ((sizeT ) / 0.5)*0.1;
		ofstream f;
		f.open("ExplicitResult");
		f << "Steps [from 0.0 sec to 0.5 sec] \n";
		for (int i = 0; i < sizeT; i = i +step)
		{
			f << "\n";
			f << "[" << i * deltaT << "]";
			for (int j = 0; j < sizeX; j++)
			{
				f<<" "<<matrixGlobal[i][j];
			}
			f << "\n";
		}
		f << "\n";
		f.close();
	}			


	if(process_id==processes-1)
	{

		//Analytical and norms
		Matrix analytical(sizeT,sizeX);

		//Initial Condition
			for (int x = 0; x < sizeX; x++)
			{
				if (x >= 0 && x <= 100)
				{
					analytical[0][x] = 0;
				}
				if (x >= 100 && x <= 220)
				{
					analytical[0][x] = 100 * (sin((PI)*((((double)x * 0.5) - 50) / 60)));
				}
				if (x >= 220 && x < sizeX)
				{
					analytical[0][x] = 0;
				}
			}

			//Boundry Condition
			for (int t = 0; t < sizeT; t++)
			{
				analytical[t][0] = 0;
				analytical[t][sizeX - 1] = 0;
			}

		double T;
		for (int i = 0; i < sizeT; i++)
		{
			T = i * deltaT;
			for (int x = 0; x < sizeX; x++)
			{
				double Y = 5 * x;
				if (Y <= 110 + (250 * T) && Y >= 50 + 250 * T)
				{
					analytical[i][x] = 100 * (sin((PI)*
						((Y - 50 - (250 * T)) / 60)));
				}
			}
		}

		//Subtraction: m1- m2
		//analytical solution - numerical solution
		Matrix result = Matrix(analytical.getNrows(), matrixGlobal.getNcols());
		result = analytical - matrixGlobal;
		double size = analytical.getNcols()*matrixGlobal.getNrows();

		
		//Print result
		cout << "\n--NORMS--";
		cout<<"\n One norm: "<<result.one_norm()/size;
		cout<<"\n Second norm: "<<result.two_norm()/size;
		cout<<"\n Uniform norm: "<<result.uniform_norm()/size;
		cout << "\n";	


	}
	
	MPI_Finalize();
}




