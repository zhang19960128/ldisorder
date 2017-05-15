#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H
#include <iostream>
#include <vector>
typedef struct particle{
	int listaround[6]={-1,-1,-1,-1,-1,-1};// store the particles around it
    int tid;//use the tid to mark each particle.
    int tidx;//use the tidx to mark each particle.
    int tidy;//use the tidy to mark each particle.
	double x;//give the x-coordinate of each particle
	double y;//give the y-coordinate of each particle
}node;
typedef struct{
    double k;
    int number;//to store nx^2+ny^2;
    std::vector<std::vector<int> > pairs;// to store pairs of (nx,ny);
}k_vector;
#endif
