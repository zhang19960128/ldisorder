#ifndef PERIOD_H
#define PERIOD_H
#include "math.h"
double disx(double x,double x1,double scale){
    //x is the center,x1-x under periodic condition;
    return (x1-x)-round((x1-x)/scale)*scale;
}
double disy(double y,double y1,double scale){
    // y is the center, y1-y under periodic condition;
    return (y1-y)-round((y1-y)/scale)*scale;
}
#endif
