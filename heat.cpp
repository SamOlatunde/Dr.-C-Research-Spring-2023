//*****************************************************************
// End of Semester Project 
// Name: Samuel Olatunde , and Sunil Rasaily 
// GPU Programming Date: Date of Submission (11/28/2022)
//******************************************************************
// Computes the temparture distrubtion of a square metal sheet with 
// edge temperatures known
//******************************************************************
#include<stdio.h>
#include<mpi.h>
#include<timer.h>
#include<stdlib.h>



//error tolerance 
const float eT  = 0.00001;
-
// data size
#define N 4

// limit for the max number of iterations
#define limit 100

double checkSum(float * h)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            sum += h[i*N+j];
        }
    }

    return sum;
}

void print(float * a)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f\t", a[i* N + j]);
        }

        printf("\n\n");
    }
}

struct
{
  int p; // total number of processes
  MPI_Comm comm; //communicator for entire grid
  MPI_Comm row_comm; //Communicator for my row
  MPI_Comm col_comm;//communicator for my col

 //Edge point flags 
  int first_row;
  int last_row;
  int first_column;
  int last_column;

  int q; // order of grid
  int my_row; //my row number 
  int my_col// my column number
  int my_rank;  // my rank in the grid communicator 
  MPI_Comm first_row;
  MPI_Comm last_row;
  MPI_Comm first_column;
  MPI_Comm last_column;
}

void Setup_grid(GRID_INFO_TYPE * grid);
//Protypes 
void initMetalPlate(float *h, float * g,float edgeTemp, int N);
void calcIntTempDistribution(float * h,float *g);
int converged (float newValue, float oldValue );

int main(int argc, char ** argv)
{
    MPI_Init(&argc,&argv);
   /* MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);*/
    

    GRID_INFO_TYPE grid;
    Setup_grid(&grid);

    int source;
    int dest;
    int tag=23;
    MPI_Status status;
    int i = grid.my_row, j = grid.my_col;

    // variable declarations
   float /** h,*/ *g ;
   float w =0.0, x = 0.0, y = 0.0, z = 0.0;
  
   // Allocate Dynamic Memory in host 
  // h = (float *) malloc ((N*N) * sizeof(float));
   g = (float *) malloc((grid.p) * sizeof(float));

   float edgeTemp =70.5;//300;
   double tStart = 0.0, tStop = 0.0, tElapsed = 0.0;
    
   //initialize matrix
   initMetalPlate(h,g, edgeTemp, grid.q);
   
   // time computation
   GET_TIME(tStart);

   // need to find a way to test for convergence

    
    if (gird.first_row) 
    {
      w = edgeTemp;//top_value;
    }
    else 
    {
      w = g[i-1][j];
    }

    if (gird.last_row) 
    {
      x = edgeTemp;//bottom_value;
    }
    else
    {
      x = g[i+1][j];
    }

    if (gird.first_column) 
    {
      y = edgeTemp;//left_value;
    }
    else 
    {
      y = g[i][j-1];
    }

    if (gird.last_column) 
    {
      z = edgeTemp;//right_value;
    }
    else 
    {
      z = g[i][j+1];
    }

    iteration = 0; /* for process i,j */

    do
     {
        iteration++;
        g[i][j] = 0.25 * (w + x+ y + z);
        if !(gird.first_row) send(&g[i][j],1, MPI_FLOAT, P_i-1,tag, grid.comm);
        if !(gird.last_row) send(&g[i][j],1, MPI_FLOAT,P_i+1,tag, grid.comm);
        if !(gird.first_column) send(&g[i][j],1, MPI_FLOAT, P_i,j-1,tag, grid.comm);
        if !(gird.last_column) send(&g[i][j],1, MPI_FLOAT,P_i,j+1, tag, grid.comm);
        if !(gird.first_row) recv(&w,1, MPI_FLOAT, P_i-1,j, tag, grid.comm, status);
        if !(gird.last_row) recv(&x,1, MPI_FLOAT, P_i+1,j, tag, grid.comm,status);
        if !(gird.first_column) recv(&y,1, MPI_FLOAT, P_i,j-1, tag, grid.comm,status), ;
        if !(gird.last_column) recv(&z,1, MPI_FLOAT, P_i,j+1, tag, grid.comm,status);
      } while (/*(!converged(i,j)) &&*/ (iteration < limit));
    //send(&g, &i, &j, &iteration, P_Master);
   
  //Need help with partitioning 
   

    //  printf("g: \n");
    //  print(g);

    //  printf("h: \n");
    //  print(h);
    
   //calcIntTempDistribution(h,g);
   GET_TIME(tStop);
   
   // Compute how long it took
   tElapsed = tStop - tStart;
   
   printf("The code to be timed took %e seconds\n", tElapsed);
   //printf("checkSum: %f\n", checkSum(h));
   //print(h);

   // Dallocate dynamic memory
   
   MPI_Finalize();

    return 0;
}


//*******************************************************************
// Name::Setup_grid()
// Parameters: 
// 
//********************************************************************
void Setup_grid(GRID_INFO_TYPE * grid)
{
  int old_rank;
   int dimensions[2];
   int periods[2];
   int reorder =1;
   int coordinates[2];
   
   //Set up Global Grid Information
   MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
   MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);
   grid->q = (int) sqrt((double)grid->p);
   dimensions[0] = dimensions[1] = grid->q;
   periods[0] = periods[1] = 0;
   MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,periods,reorder,&(grid->comm));
   MPI_Comm_rank(grid->comm, &(grid->my_rank));
   MPI_Cart_coords(grid->comm,grid->my_rank,2,coordinates);
   grid->my_row = coordinates[0];
   grid->my_col = coordinates[1];

  //Set up row and column communicators
  varying_coords[0] = 0; varying_coords[1] = 1;
  MPI_Cart_sub(grid->comm, varying_coords, &(grid->row_comm));
  varying_coords[0] = 1; varying_coords[1] = 0;
  MPI_Cart_sub(grid->comm,varying_coords, &(grid->col_comm));
  
  grid->first_row = grid->last_row = grid->first_column = grid->last_column= 0;

  if (grid->my_row == 0)
  {
     grid->first_row = 1;
  }
  else if (grid->my_row == (grid->q-1))
  {
     grid->last_row = 1;
  }


  if (grid->my_col == 0)
  {
     grid->first_column = 1;
  }
  else if (grid->my_col == (grid->q-1))
  {
     grid->last_column = 1;
  }
}

//*******************************************************************
// Name::calcIntTempDistribution()
// Parameters: 2 float pointers
// Calculates the temparture of interior point by find the avergae 
// of the four adjacent points 
//********************************************************************
void calcIntTempDistribution(float *h,float *g)
{
   int iteration = 0;
   
   int Continue;

   do 
   {
      //compute averages
      for (int i = 1; i < (N-1); i++)
      {
        for(int j = 1; j < (N-1); j++)
        {
            g[i* N + j] = 0.25 * (h[(i-1) * N + j] + h[(i+1) * N + j]
                                       +h[i* N + j-1]+h[i* N + j+1]);
        }
      }
      
      Continue = 0;
      
      // test convergence and update new value array
      for (int i = 1; i < (N-1); i++)
      {
        for (int j = 1; j<(N-1); j++)
        {
            if( converged(g[i*N + j],h[i* N + j]) == 0)
            {
                Continue = 1;
            }
            h[i* N + j] = g[i * N + j];
        }
      }
    //  printf("g: \n");
    //  print(g);

    //  printf("h: \n");
    //  print(h);
     iteration++;
   }while(Continue == 1 && iteration < limit);
//    printf("Blah");
//    printf("%d\n", iteration);
}

//*******************************************************************
// Name::converged()
// Parameters: 2 floats
// Tests for convergence of two points. Returns true if the error is 
// within error tolerance; false otherwise 
//********************************************************************
// bool converged (float newValue, float oldValue )
 int converged (float newValue, float oldValue )
{
    float er = (newValue-oldValue)/newValue;
    //printf("er %f\n", er);
    if (er < 0) er = -er;

    if (er <= eT) 
    {
        return 1;
    }
    else 
    {
        return 0;
    }
}


//*******************************************************************
// Name::initMetalPlate()
// Parameters: 1 2d float array, 1 float
// Initializes the metal sheet with the intial values of the edges
// and guess values for interior points 
//********************************************************************
void initMetalPlate(float *h, float * g, float edgeTemp, int N)
{
   //we reduce the temparture by this value with every 
   // outer loop iteration
   float reduceFactor = edgeTemp/N;

   int row = 0;
   int col;

   for( int i = 0; i < (N/2); i++)
   {
        
        row = i;
        for (col = i; col < N-i; col++ )
        {
            h[row * N + col] = edgeTemp;
            g[row * N + col] = edgeTemp;
        }

        col--;

        for (row = row+1; row < N-i; row++)
        {
            h[row * N + col] = edgeTemp;
            g[row * N + col] = edgeTemp;
        }

        row--;
        col--;

        for(col = col; col >=i; col--)
        {
           h[row * N + col] = edgeTemp;
           g[row * N + col] = edgeTemp;
        }
        
        row--;
        col++;

        for(row = row; row >i; row--)
        {
            h[row * N + col] = edgeTemp;
            g[row * N + col] = edgeTemp;
        }
      
      edgeTemp = edgeTemp - reduceFactor;
    }

    // print(g);
    // printf("\n\n");
    
}


