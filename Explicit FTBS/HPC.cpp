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

//Create vectors of coordinates for each sub-domains in matrix
void mapping(vector<int> &x_start, vector<int> &t_start, vector<int> &x_end, vector<int> &t_end, int cell_x, int cell_t, int x, int t)
{

	//Compute the points for space dimension - for each colum.
	if(x==1)
	{
		x_start[x-1] = 1;
		x_end[x-1]=cell_x-2;
	}
	if(x==2)
	{
		x_start[0] = 1;
		x_end[0] = cell_x-1;

		x_start[1] = cell_x;
		x_end[1] = x_start[1]+cell_x - 1;
	}
	else
	{
		int k;
		int tmp2 = 0;
		x_start[0] = 1;
		x_end[0] = cell_x-1;
		for(k=1;k<x-1;k++)
		{	
			tmp2+=cell_x;
			x_start[k] = tmp2;
			x_end[k] = x_start[k]+cell_x - 1;
					
		}
		tmp2+=cell_x;
		x_start[x-1] = tmp2;
		x_end[x-1] = x_start[x-1]+cell_x-2; 
	}
	
	//Compute the points for time dimension - for each row.
	int i;
	int tmp =0;
	t_start[0]=1;
	t_end[0]=cell_t-1;
	for(i=1;i<t;i++)
	{
		tmp+=cell_t;
		t_start[i] = tmp;
		t_end[i]=t_start[i]+cell_t-1;
		
	}
}

// Update matrix: add new computed values
void merge(Matrix &matrix, Matrix &prevMatrix, int x_s, int x_e, int t_s, int t_e)
{
	for (int t = t_s; t <= t_e; t++)
	{ 				
			for (int x = x_s; x <= x_e; x++)
			{
				(matrix)[t][x] = prevMatrix[t][x];
			}
	}
}



int main(int argc, char* argv[])
{
	//MPI part
	double t1, t2;		//variable for time performance
	int process_id;		//process_id: label for process, 
	int processes; 		//processes: number of processes from console
	MPI_Init(&argc, &argv);
	
	MPI_Comm communicator;
	communicator = MPI_COMM_WORLD;
	MPI_Comm_size(communicator, &processes);
	MPI_Comm_rank(communicator, &process_id);
	MPI_Status status;


	string s = argv[1];
	double deltaT = stod(s);
	const double PI = atan(1) * 4; //Constant double value PI = atan(1) * 4
	const double u = 250;	//Constant double value from excercise 250 m/s
	double deltaX = 0.5;	//The same for each deltaT
	const double x = 400;	//Domain is xE[0;400] meters
	const double t = 0.5;	//t goes from 0.0 to 0.5sec

	int sizeX; 
	int sizeT;
	

	int periods[2];
	int dims[2]; 
	int reorder = 1;
	int dimension = 2;
	int tmp_dims[2];

	//Create Cartesian grid
	periods[0] = 0;
	periods[1] = 0;
	MPI_Dims_create(processes, dimension, dims);
	
	//swap
	tmp_dims[0] = dims[0];
	tmp_dims[1] = dims[1];
	dims[0] = tmp_dims[1];
	dims[1] = tmp_dims[0];
		

	//compute sizeX (columns) and sizeT (rows), this is size of matrix
	sizeX = (x / deltaX);
	sizeT = (t / deltaT);
	Matrix matrix(sizeT,sizeX);

	//size for each cell in grid
	int cell_x = (sizeX/dims[0]);
	int cell_t = (sizeT/dims[1]);


	//Process to map: compute start and end point for each cell
	//return vectors with points start and end
	vector<int> x_start, x_end, t_start, t_end;
	x_start.resize(dims[0]);
	x_end.resize(dims[0]);
	t_start.resize(dims[1]);
	t_end.resize(dims[1]);
	mapping(x_start, t_start, x_end, t_end, cell_x, cell_t, dims[0], dims[1]);


	// Resize Matrix, fill by zeros, inital and boundry conditions
	if (process_id == 0)
	{

		//Initial Condition
		for (int x = 0; x < sizeX; x++)
		{
			if (x >= 0 && x <= 100)
			{
				matrix[0][x] = 0;
			}
			if (x >= 100 && x <= 220)
			{
				matrix[0][x] = 100 * (sin((PI)*((((double)x * 0.5) - 50) / 60)));
			}
			if (x >= 220 && x < sizeX)
			{
				matrix[0][x] = 0;
			}
		}

		//Boundry Condition
		for (int t = 0; t < sizeT; t++)
		{
			matrix[t][0] = 0;
			matrix[t][sizeX - 1] = 0;
		}

	}

	t1 = MPI_Wtime();	//start measure time

	//P0 solve zero part of matrix.
	if(process_id == 0)
	{
		for (int t = t_start[0]; t <= t_end[0]; t++)
		{
			for (int x = x_start[0]; x <= x_end[0]; x++)
			{

				matrix[t][x] = (1 - ((u*deltaT) / deltaX))*matrix[t - 1][x] 
				+((u*deltaT) / deltaX)*matrix[t - 1][x - 1];

			}		
		}
	}

	//P0 - root process - sends and receives matrix. 
	//Start from first part of matrix.
	if(process_id == 0)
	{	
		//temporary matrix
		Matrix prevMatrix(sizeT,sizeX);
		prevMatrix = matrix;

		//vector with coordinates
		vector<pair<int,int>> cells;
		
		int j = 0, i = 1, itmp = 1, jtmp = 0;
		int currentProcessSend = 1;
		int currentProcessRecv = 1;
		
		//while each cells/parts of the matrix have been solved.
	loop:	while((i+1)*(j+1)<=((dims[0]*dims[1])))
		{
			cells.push_back(make_pair(itmp,jtmp));
			itmp--;

			if(itmp < 0)
			{
				for(int s = 0; s < cells.size(); s++)
				{
					//Send the matrix to current process
					MPI_Send(&(matrix[0][0]),sizeT*sizeX,MPI_DOUBLE,
					currentProcessSend, 1234, communicator);
					//Send the coordinates to current process
					MPI_Send(&(cells)[s],2,MPI_INT,
					currentProcessSend, 12345, communicator);

					currentProcessSend++;
				}
				for(int r = 0; r < cells.size(); r++)
				{
					
					MPI_Recv(&(matrix[0][0]),sizeT*sizeX,MPI_DOUBLE,
					currentProcessRecv,666,communicator,&status);

					merge(prevMatrix, matrix,
					x_start[cells[r].first],x_end[cells[r].first],
					t_start[cells[r].second],t_end[cells[r].second]);

					currentProcessRecv++; 
					matrix=prevMatrix;
				}
				

				cells.clear();
				i++;		
				if(i < dims[0])
				{
					itmp = i;
					jtmp = 0;
				}
				else
				{
					i = dims[0] - 1;
					itmp = dims[0] - 1;
					j++;
					jtmp=j;	
				}
				goto loop;
			}
			else
			{
				jtmp++;
				if(jtmp > dims[1] - 1)
				{
				for(int s = 0; s < cells.size(); s++)
				{
						
					MPI_Send(&(matrix[0][0]),sizeT*sizeX,MPI_DOUBLE,
					currentProcessSend, 1234, communicator);
					MPI_Send(&(cells)[s],2,MPI_INT,currentProcessSend,
					12345, communicator);
					currentProcessSend++;
				}
				for(int r = 0; r < cells.size(); r++)
				{

					MPI_Recv(&(matrix[0][0]),sizeT*sizeX,MPI_DOUBLE,
					currentProcessRecv,666,communicator,&status);

					//merge function: update the matrix with newly calculated parts
					merge(prevMatrix, matrix,x_start[cells[r].first],x_end[cells[r].first],
					t_start[cells[r].second],t_end[cells[r].second]);
					currentProcessRecv++;
					matrix=prevMatrix;
				}

				cells.clear();
				j++;
				jtmp = j;
				itmp = dims[0] - 1;
				goto loop;
				}
			}
				
		}

		//SAVE RESULTS INTO FILE		

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
				f<<" "<<matrix[i][j];
			}
			f << "\n";
		}
		f << "\n";
		f.close();

	}


	//Explicit Scheme FTBS
	//Processes different than root process
	if(process_id>0)
	{
		
		int t_tmp = 0;
		int x_tmp = 0;
		int nonZeroCount=0;
		vector<pair<int,int>> cells(1);

		//Receive the matrix and coordinates from process P0
		MPI_Recv(&(matrix[0][0]),sizeT*sizeX,MPI_DOUBLE,0,1234,communicator,&status);
		MPI_Recv(&(cells)[0],2,MPI_INT,0,12345,communicator,&status);


		for (int t = t_start[cells[0].second]; t <= t_end[cells[0].second]; t++)
		{ 				
			for (int x = x_start[cells[0].first]; x <= x_end[cells[0].first]; x++)
			{
				//solve
				matrix[t][x] = (1 - ((u*deltaT) / deltaX))*matrix[t - 1][x] 
				+((u*deltaT) / deltaX)*matrix[t - 1][x - 1];
			}

				
		}

		//Send solved part to root process (P0)
		MPI_Send(&(matrix[0][0]), sizeT*sizeX, MPI_DOUBLE, 0, 666, communicator);	
		
	}

	t2 = MPI_Wtime();	//end measure time
	cout << "\n\nElapsed time is "<< t2 - t1 <<" for process: "<<process_id;
	

	if(process_id==0)
	{

		//For Analytical and norms
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
				double Y = 0.5 * x;
				if (Y <= 110 + (250 * T) && Y >= 50 + 250 * T)
				{
					analytical[i][x] = 100 * (sin((PI)*
						((Y - 50 - (250 * T)) / 60)));
				}
			}
		}

		//Subtraction: m1- m2
		//analytical solution - numerical solution
		Matrix result = Matrix(analytical.getNrows(), matrix.getNcols());
		result = analytical - matrix;
		double size = analytical.getNcols()*matrix.getNrows();

		
		//Print result
		cout << "\n--NORMS--";
		cout<<"\n One norm: "<<result.one_norm()/size;
		cout<<"\n Second norm: "<<result.two_norm()/size;
		cout<<"\n Uniform norm: "<<result.uniform_norm()/size;
		cout << "\n";	


	}

	MPI_Finalize();
}




