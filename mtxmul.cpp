#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <time.h>

int main(int argc, char *argv[])
{
	double **mtx1,**mtx2,**mtx;
	double *line;
	int i, j, k, mtxsize;
	clock_t alloc_timer, cal_timer;

	mtxsize = atoi(argv[1]);
	if (mtxsize <= 0)
	{
		std::cout<< "input argument is wrong"<<std::endl;
		return 1;
	}

	//Allocating matrix 1 and matrix 2
	alloc_timer = clock();
	mtx1 = new double * [mtxsize];
	mtx2 = new double * [mtxsize];
	mtx  = new double * [mtxsize];
	srand(time(NULL));
	for(i=0;i<mtxsize;i++)
	{
		try
		{
			*(mtx1+i) = new double [mtxsize];
        	*(mtx2+i) = new double [mtxsize];
			*(mtx+i)  = new double [mtxsize];
		}
		catch (std::exception& excpt)
		{
			std::cout<< "Standard exception: " << excpt.what() <<std::endl;
			return 1;
		}
        for(j=0;j<mtxsize;j++)
		{
			*(*(mtx1+i)+j) = double(rand()-RAND_MAX/2.0)/double(RAND_MAX);
			*(*(mtx2+i)+j) = double(rand()-RAND_MAX/2.0)/double(RAND_MAX);
			*(*(mtx +i)+j) = 0.0;
		}
	}
	alloc_timer = clock()-alloc_timer;

    // print out matrices mtx1 and mtx2
	/*
	for(i=0;i<mtxsize;i++)
	{
		std::cout<<"[ ";
		std::cout<<std::setw(10);
		for(j=0;j<mtxsize;j++)
			std::cout<<*(*(mtx1+i)+j)<<" ";
		std::cout<<"]    [";
		for(j=0;j<mtxsize;j++)
			std::cout<<*(*(mtx2+i)+j)<<" ";
		std::cout<<"]"<<std::endl;
	}
    std::cout<<"="<<std::endl;
	*/

	//calculating mtx multiplication
	cal_timer = clock();
	for(i=0;i<mtxsize;i++)
	{
		//std::cout<<"[";
		for(j=0;j<mtxsize;j++)
		{
			for(k=0;k<mtxsize;k++)
			{
				*(*(mtx+i)+j) = *(*(mtx+i)+j) + (*(*(mtx1+i)+k)) * (*(*(mtx2+k)+j));
			}
			//std::cout<< *(*(mtx+i)+j) << " ";
		}
		//std::cout<<"]"<<std::endl;
	}
	cal_timer = clock()-cal_timer;

	std::printf("Allocating time: %12.6f \n",(double)alloc_timer/(double)CLOCKS_PER_SEC);  
	std::printf("Calculating time: %12.6f \n",(double)cal_timer/(double)CLOCKS_PER_SEC);

    for(i=0;i<mtxsize;i++)
	{
		delete [] *(mtx+i) ;
		delete [] *(mtx1+i);
		delete [] *(mtx2+i);
	}
	delete [] mtx;
	delete [] mtx1;
	delete [] mtx2;

	return 0;
}
