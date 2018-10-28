#include "Lattice.hpp"

using namespace std;

Lattice::Lattice(int L, double J, double T, int Seed, int Initial_config){

    std::uniform_real_distribution<> uniform(0.0,1.0);

    _L=L;
    _N=L*L;
    _J=J;
    _T=T;

    RND.seed(Seed);

    beta=1.0/T;
    p  = 1.0 - std::exp(-2.0*J*beta);

    n_it=0;

    //create Lattice

    _S=new int[_N];
    if(Initial_config==0){
        cout<<"Initial configuration: T=0"<<endl<<endl;
        for(int i=0; i<_N; i++)
            _S[i]=1;
    }
    else{
        cout<<"Initial configuration: T=Infinite"<<endl<<endl;
        int values[2]={-1,1};
        for(int i=0; i<_N; i++)
            _S[i]=values[int(uniform(RND)*2)];
    }

    // create neighbour list
    //a%b
    //(b + (a%b)) % b

    _nbr=new int*[_N];
    for(int i = 0; i < _N; ++i)
        _nbr[i]=new int[4];

    for(int i=0; i<_N; i++){
        _nbr[i][0]=(i/L) * L + (L+((i+1)%L))%L;
        _nbr[i][1]= (_N + ((i+L) % _N) )%_N;
        _nbr[i][2]=(i/L) * L + (L + ((i-1)%L)  )%L;
        _nbr[i][3]= (_N+((i-L) % _N ))%_N;
    }

}

int Lattice::WolffClusterMove(){

    std::uniform_real_distribution<> uniform(0.0,1.0);

    int k=uniform(RND)*_N;
    Pocket.insert(k);
    Cluster.insert(k);
    while(Pocket.size()>0){
        it=Pocket.begin();
        advance(it,uniform(RND)*Pocket.size());
        int j=*it;
        for(int n=0; n<4; n++){
            const bool nbr_not_in_C = Cluster.find(_nbr[j][n]) == Cluster.end();
            if(_S[_nbr[j][n]]==_S[j] && nbr_not_in_C && uniform(RND)<p ){
                Pocket.insert(_nbr[j][n]);
                Cluster.insert(_nbr[j][n]);
            }
        }
        Pocket.erase(j);
    }

    for (it=Cluster.begin(); it!=Cluster.end(); ++it)
        _S[*it]*=-1;

    int cluster_dim=Cluster.size();

    Cluster.clear();

    return cluster_dim;
}

int Lattice::MetropolisMove(){

    int accepted_flips=0;

    std::uniform_real_distribution<> uniform(0.0,1.0);

    for(int i=0; i<_N; i++){
        int k=uniform(RND)*_N;
        double Snn=_S[_nbr[k][0]]+_S[_nbr[k][1]]+_S[_nbr[k][2]]+_S[_nbr[k][3]];      // sum(S[nn] for nn in nbr[k])
        double delta_E=2.0*_J*_S[k]*Snn;

        if( uniform(RND) < std::exp(-beta * delta_E) ){
            _S[k] *= -1;
            accepted_flips++;
        }

    }

    return accepted_flips;

}

double Lattice::computeE(){
    double _E=0;
    double Sn;
    for(int i=0; i<_N; i++){
        Sn=_S[_nbr[i][0]]+_S[_nbr[i][1]]+_S[_nbr[i][2]]+_S[_nbr[i][3]];
        _E+=-_J*_S[i]*Sn;
    }
    return 0.5*_E/_N;
}

double Lattice::computeM(){
    double _M=0;
    for(int i=0; i<_N; i++)
        _M+=_S[i];
    return _M/_N;
}

void Lattice::PrintS(string filename){
    ofstream myfile;
    myfile.open("IMAGES/"+filename+".pgm");
    myfile<<"P2"<<endl<<_L<<" "<<_L<<endl<<"255"<<endl;
    for(int i=0; i<_N; i++){
        myfile<<(_S[i]+1)*127<<endl;
    }
    myfile.close();
}

void Lattice::PrintS(){
    string filename="output"+to_string(n_it);
    this->PrintS(filename);
    n_it++;
}

int Lattice::GetS(int pos){
    return _S[pos];
}

Lattice::~Lattice(){
    delete[] _S;
    for(int i=0; i<_N; ++i)
        delete[] _nbr[i];
    delete[] _nbr;
}
