#include "model.h"
#include <iostream>
#include <stdlib.h> 
#include <cmath>

using namespace std;


int main(int argc, char *argv[])
{
    int choice = atoi(argv[2]);
    int solver = atoi(argv[1]);
    Model myModel = Model(choice);
    double T = myModel.simulate(solver, false);
    int NumOfStep = std::floor(Model::default_NumOfStep*myModel.Tf/T);
    myModel = Model(choice,NumOfStep);
    T = myModel.simulate(solver, true);
    std::cout<<T<<" - \""<<myModel.filename<<"\""<<std::endl;
    return 0;
}

