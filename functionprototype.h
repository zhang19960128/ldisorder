#ifndef FUNCTIONPROTOTYPE_H
double k_constant(double); 
int perodic(int,int); 
void dimension_change(int *,int,int);
void ad_kconstant(node **,double **,int,double);
void listinsert(std::list<k_vector> &,k_vector);
std::list<k_vector> sort(std::vector<k_vector> );
double sumFnl(double **,k_vector,node *,int m,int n);//node contain all the information of the particles.double ** is the e_nj,k_vector is the vector;
double sumFnt(double **,k_vector,node *,int m,int n);
double* generateposition(double ita);
#endif
