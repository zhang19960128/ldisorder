#include "datastructure.h"
#include <iostream>
#include "mkl_lapacke.h"
#include "myrand.h"
#include <cmath>
#include <complex>
#include <list>
extern lapack_int LAPACKE_dsyev(int matrix_layout,char jobz,char uplo,lapack_int n,double* a,lapack_int lda,double* w);
double k_constant(double inta){
    return genrand()/2.0*inta+1.000;
}
int perodic(int x,int m){
    if(x<0) return x+m;
    else if(x>=m) return x-m;
    else return x;
}
void dimension_change(int *a,int j,int Dm){
    a[1]=j%Dm;//i
    a[0]=(j-j%Dm)/Dm;//j;
}
void ad_kconstant(node **particle_node_II,double **ad_matrix,int n,double ita){
     int a;
     int m=sqrt(n);
     double k_con;
     for(int i=0;i<m;i++)
         for(int j=0;j<m;j++){
             a=i*m+j;
             for(int k=0;k<6;k++){
                 if(particle_node_II[i][j].listaround[k]==-1) continue;
                 if(ad_matrix[particle_node_II[i][j].listaround[k]][a]==0){
                    k_con=1+ita*(genrand()-0.5);
                    ad_matrix[particle_node_II[i][j].listaround[k]][a]=k_con;
                    ad_matrix[a][particle_node_II[i][j].listaround[k]]=k_con;
                 }
                 else{
                     ad_matrix[a][particle_node_II[i][j].listaround[k]]=ad_matrix[particle_node_II[i][j].listaround[k]][a];
                 }
             }
         }
};
void listinsert(std::list<k_vector> &k_list,k_vector in){
    for(auto a=k_list.begin();a!=k_list.end();a++){
        if(a->number==in.number){
            (a->pairs).push_back(*(in.pairs.begin()));
            return;
        }
        else if(in.number < a->number){
               k_list.insert(a,in);
            return;
               }
             
    }
    k_list.push_back(in);
}
std::list<k_vector> sort(std::vector<k_vector> k){
    std::list<k_vector> k_list;
     k_list.push_back(*(k.begin()));
     for(auto a=(k.begin()+1);a!=k.end();a++){
         listinsert(k_list,*a);
     }
     return k_list;
}
double sumFnl(double ** matrix,k_vector kv,node* particle_node_I,int m,int n){
//matrix is the eigenvector of the hessian matrix; particle_node_II constains the information of all the paticles;
//m is the length of each side.
//use n to represent the nth mode
    double result=0;
    double qx,qy,ex,ey;
    std::complex<double> result1(0,0);
    std::complex<double> I(0,1);
    std::complex<double> delay(0,0);
    for(auto a=kv.pairs.begin();a!=kv.pairs.end();a++){
         qx=(*a)[0]*2*3.141592653/m;
         qy=(*a)[1]*2*3.141592653/m;
         result1=0;
         for(int i=0;i<2*m*m;i=i+2){
             ex=matrix[i][n];
             ey=matrix[i+1][n];
             delay=I*(qx*particle_node_I[i/2].x+qy*particle_node_I[i/2].y);
             if(sqrt(qx*qx+qy*qy)==0){
                 continue;
             }
             result1=(ex*qx+ey*qy)/sqrt(qx*qx+qy*qy)*std::exp(delay)+result1;
            // std::cout<<result1<<std::endl;
         }
         result=result+std::abs(result1)*std::abs(result1);
    }
    return result;
}
double sumFnt(double ** matrix,k_vector kv,node* particle_node_I,int m,int n){
//matrix is the eigenvector of the hessian matrix; particle_node_II constains the information of all the paticles;
//m is the length of each side.
//use n to represent the nth mode
    double result=0;
    double qx,qy,ex,ey;
    std::complex<double> result1(0,0);
    std::complex<double> I(0,1);
    std::complex<double> delay(0,0);
    for(auto a=kv.pairs.begin();a!=kv.pairs.end();a++){
         qx=(*a)[0]*2*3.141592653/m;
         qy=(*a)[1]*2*3.141592653/m;
         result1=0;
         for(int i=0;i<2*m*m;i=i+2){
             ex=matrix[i][n];
             ey=matrix[i+1][n];
             delay=I*(qx*particle_node_I[i/2].x+qy*particle_node_I[i/2].y);
             if(sqrt(qx*qx+qy*qy)==0){
                 continue;
             }
             result1=(ex*qy-ey*qx)/sqrt(qx*qx+qy*qy)*std::exp(delay)+result1;
            // std::cout<<result1<<std::endl;
         }
         result=result+std::abs(result1)*std::abs(result1);
    }
    return result;
}
double* generateposition(double ita){
    double* a=new double[2];
    double chai=genrand()-0.5;
    double theta=(genrand()-0.5)*3.141592653;
    a[0]=ita*chai*cos(theta);
    a[1]=ita*chai*sin(theta);
    return a;
}
