#include <iostream>
#include <fstream>
#include <list>
#include "datastructure.h" 
#include "functionprototype.h"
#include <stdio.h>
#include "mkl_lapacke.h"
#include <malloc.h>
#include <new>
#include <list>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "period.h"
#define pi 3.141592653
int main(int argc,char* argv[]){
	int n=4096;//n is the number of the particles in our system
	int m = sqrt(n);
    double ita=atof(argv[1]);
	node **particle_node_II = new node*[m];//store the point in two dimensional array
	node *particle_node_I = new node[m*m];//store the point in one dimensional array
	for (int i = 0; i < m; i++){
		particle_node_II[i] = &particle_node_I[i*m];
	}
    double *p=NULL;
    for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
            p=generateposition(ita);
            particle_node_II[i][j].y=i*sqrt(3)/2.0+p[1];
            if(i%2==1) {
                particle_node_II[i][j].x=0.5+j*1.0+p[0];
                particle_node_II[i][j].listaround[0]=perodic(i-1,m)*m+j;
                particle_node_II[i][j].listaround[1]=perodic(i-1,m)*m+perodic(j+1,m);
                particle_node_II[i][j].listaround[2]=(i)*m+perodic(j-1,m);
                particle_node_II[i][j].listaround[3]=(i)*m+perodic(j+1,m);
                particle_node_II[i][j].listaround[4]=perodic(i+1,m)*m+j;
                particle_node_II[i][j].listaround[5]=perodic(i+1,m)*m+perodic(j+1,m);
            }
            else {
                particle_node_II[i][j].x=j*1.0+p[0];
                particle_node_II[i][j].listaround[0]=perodic(i-1,m)*m+perodic(j-1,m);
                particle_node_II[i][j].listaround[1]=perodic(i-1,m)*m+j;
                particle_node_II[i][j].listaround[2]=(i)*m+perodic(j-1,m);
                particle_node_II[i][j].listaround[3]=(i)*m+perodic(j+1,m);
                particle_node_II[i][j].listaround[4]=perodic(i+1,m)*m+perodic(j-1,m);
                particle_node_II[i][j].listaround[5]=perodic(i+1,m)*m+j;

            }
            particle_node_II[i][j].tidx=j;
            particle_node_II[i][j].tidy=i;
            particle_node_II[i][j].tid=i*m+j;
         }
    }
    //generate a adjacency matrix.
    double** ad_matrix=new double* [n];
    double* ad_matrix_one=new double [n*n];
    for(int i=0;i<n;i++) ad_matrix[i]=&ad_matrix_one[i*n];
    //give the initial value of the adjacency matrix;
    for(int i=0;i<n;i++){
        for(int j=i;j<n;j++){
            ad_matrix[i][j]=0;
            ad_matrix[j][i]=0;
        }
    }
    ad_kconstant(particle_node_II,ad_matrix,n,ita);
    //calcualte the hessian matrix
    double** matrix=new double* [2*n];
    double* matrix_oneD=new double[2*n*2*n];
    double* lamda=new double[2*n];
    for(int i=0;i<2*n;i++) matrix[i]=&matrix_oneD[2*n*i];
    for(int i=0;i<2*n;i++)
        for(int j=0;j<2*n;j++)
            matrix[i][j]=0;
    int a[2];
    int a1[2];
    int b;
    double distance;
    double distance_x;
    double distance_y;
    double matrix_xx;//the elements in the matrix;
    double matrix_yy;//the elements in the matrix;
    double matrix_xy;//the elements in the matrix;
    for(int i=0;i<n;i++){
        dimension_change(a,i,m);
        matrix_xx=0;
        matrix_yy=0;
        matrix_xy=0;
        for(int k=0;k<6;k++){
           b=particle_node_II[a[0]][a[1]].listaround[k];
           dimension_change(a1,b,m);
           if(b==-1) continue;
           distance_x=disx(particle_node_II[a[0]][a[1]].x,particle_node_II[a1[0]][a1[1]].x,(m));
           distance_y=disy(particle_node_II[a[0]][a[1]].y,particle_node_II[a1[0]][a1[1]].y,(m)*sqrt(3)/2);
           distance=sqrt(distance_x*distance_x+distance_y*distance_y);
           matrix_xx=ad_matrix[b][i]*distance_x*distance_x/distance/distance+matrix_xx;
           matrix_yy=ad_matrix[b][i]*distance_y*distance_y/distance/distance+matrix_yy;
           matrix_xy=ad_matrix[b][i]*distance_x*distance_y/distance/distance+matrix_xy;
           matrix[2*i][2*b]=-1*ad_matrix[b][i]*distance_x*distance_x/distance/distance;
           matrix[2*i][2*b+1]=-1*ad_matrix[b][i]*distance_x*distance_y/distance/distance;
           matrix[2*i+1][2*b]=matrix[2*i][2*b+1];
           matrix[2*i+1][2*b+1]=-1*ad_matrix[b][i]*distance_y*distance_y/distance/distance;
         }
        matrix[2*i][2*i]=matrix_xx;
        matrix[2*i][2*i+1]=matrix_xy;
        matrix[2*i+1][2*i]=matrix_xy;
        matrix[2*i+1][2*i+1]=matrix_yy;
    }
     lapack_int info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',2*n,matrix_oneD,2*n,lamda);
     if(info==1){
        std::cout<<"LAPACKE_ERROR"<<std::endl;
        return 0;
     }
    for(int i=0;i<2*n;i++){
        lamda[i]=sqrt(lamda[i]);
    }
    std::fstream eigen;
    eigen.open(argv[2],std::fstream::out);
    for(int i=0;i<2*n;i++){
        eigen<<lamda[i]<<std::endl;
    }
    eigen.close();
}
