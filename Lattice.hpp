#ifndef __Lattice__
#define __Lattice__

#include <random>
#include <set>
#include <iostream>
#include <fstream>

class Lattice{
public:
    Lattice(int L, double J, double T, int Seed, int Initial_config);
    ~Lattice();
    int MetropolisMove();
    int WolffClusterMove();
    double computeE();
    double computeM();
    void PrintS();
    void PrintS(std::string);
    int GetS(int pos);

private:
    int* _S;
    int** _nbr;
    int _N;
    double _J;
    double _T;
    int _L;

    double p;
    double beta;

    int n_it;

    std::mt19937 RND;

    std::set<int> Pocket;
    std::set<int> Cluster;
    std::set<int>::const_iterator it;

};

#endif
