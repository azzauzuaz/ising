#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "Lattice.hpp"
#include "Histogram.hpp"
#include "ising.hpp"

using namespace std;

int main(){

    Input();

    cout<<"Equilibration..."<<endl;

    begin_time = clock();

    for(int eq_step=0; eq_step<N_EQUILIBRATION_WOLFF_STEPS; eq_step++){
        cluster_dim=S->WolffClusterMove();
        hist->add_x(cluster_dim);
    }
    for(int eq_step=0; eq_step<N_EQUILIBRATION_METROPOLIS_STEPS; eq_step++){
        S->MetropolisMove();
    }

    S->PrintS("equilibrated");

    end_time = clock();
    double elapsed_secs = double(end_time - begin_time) / CLOCKS_PER_SEC;

    cout<<"Equilibration time: "<<elapsed_secs<<endl<<endl;

    hist->print_hist("cluster_histogram.dat");

    cout<<"Simulation..."<<endl<<endl;
    for(int iblock=0; iblock<N_BLOCKS; iblock++){

        Reset();

        for(int step=0; step<N_WOLFF_STEPS_PER_BLOCK; step++){
            cluster_dim=S->WolffClusterMove();
            hist->add_x(cluster_dim);
            Measure();
            Accumulate();
        }
        for(int step=0; step<N_METROPOLIS_STEPS_PER_BLOCK; step++){
            accepted_flips+=S->MetropolisMove();
            Measure();
            Accumulate();
        }

        S->PrintS();
        Average();

    }

    PrintFinal();

    return 0;
}

//##############################################################################

void Input(){

    cout << "Ising Model " << endl;
    cout << "Monte Carlo simulation  " << endl<<endl;

    Read_From_File("input.conf");

    cout << "Simulation parameters: "<<endl;
    cout << "Lattice dimension: "<< L <<" X "<<L<<endl;
    cout << "Temperature: "<<T<<endl;
    cout << "Interaction constant J: "<<J<<endl;
    cout << "Equilibration steps:"<<endl;
    cout << N_EQUILIBRATION_WOLFF_STEPS << " WOLFF"<<endl;
    cout << N_EQUILIBRATION_METROPOLIS_STEPS << " METROPOLIS" <<endl;
    cout << "Number of Blocks: "<<N_BLOCKS<<endl;
    cout << "Steps per block: "<<endl;
    cout << N_WOLFF_STEPS_PER_BLOCK << " WOLFF" <<endl;
    cout << N_METROPOLIS_STEPS_PER_BLOCK << " METROPOLIS" <<endl;

    StepPerBlock=N_WOLFF_STEPS_PER_BLOCK+N_METROPOLIS_STEPS_PER_BLOCK;

    S=new Lattice(L, J, T, IDUM, INITIAL_CONFIG);

    N=L*L;

    hist=new Histogram(0., N, 100, true);

    sisj=new double[R_CORR];
    sum_sisj=new double[R_CORR];
    corr_func=new double[R_CORR];
    avg_corr_func=new double[R_CORR];
    for(int r=0; r<R_CORR; r++){
        avg_corr_func[r]=0;
    }

    m_history=new double[T_CORR];
    t_corr_func=new double[T_CORR];
    avg_t_corr_func=new double[T_CORR];
    for(int t=0; t<T_CORR; t++){
        avg_t_corr_func[t]=0;
    }

    S->PrintS("0start");
}

void Reset(){

    begin_time = clock();

    sum_E=0;
    sum_M=0;
    sum_E2=0;
    sum_M2=0;

    for(int r=0; r<R_CORR; r++){
        sum_sisj[r]=0;
    }
    for(int t=0; t<T_CORR; t++){
        t_corr_func[t]=0;
    }
}

void Measure(){

    _E=S->computeE();
    _M=S->computeM();
    _E2=_E*_E;
    _M2=_M*_M;

    for(int r=0; r<R_CORR; r++){
        sisj[r]=0;
    }

    for(int icorr=0; icorr<N; icorr++){
        int si=S->GetS(icorr);

        for(int r=0; r<R_CORR; r++){
            int _r=r+1;
            double rx=S->GetS((icorr/L) * L + (L+((icorr+_r)%L))%L);
            double down=S->GetS((N + ((icorr+_r*L) % N) )%N);
            double lx=S->GetS((icorr/L) * L + (L + ((icorr-_r)%L)  )%L);
            double up=S->GetS((N+((icorr-_r*L) % N ))%N);

            double sj=(up+down+lx+rx)/4.;
            sisj[r]+=si*sj/N;
        }
    }

    //spostare dentro Accumulate ?
    //nota che se T_CORR è maggiore del numero di steps per block ci sono problemi
    if(block_norm>0){
        double m_avg=global_M/block_norm;
        for(int t=0; t<T_CORR; t++)
            t_corr_func[t]+=(_M-m_avg)*(m_history[ (T_CORR+( (t_count-t-1)%T_CORR ))%T_CORR ]-m_avg);

    }

    m_history[t_count]=_M;
    t_count++;
    if(t_count==T_CORR)
        t_count=0;

}

void Accumulate(){
    sum_E+=_E;
    sum_M+=_M;
    sum_E2+=_E2;
    sum_M2+=_M2;

    for(int r=0; r<R_CORR; r++){
        sum_sisj[r]+=sisj[r];
    }
}

void Average(){

    int wd=12;

    double energy=sum_E/StepPerBlock;
    double magnetization=sum_M/StepPerBlock;
    double heat_capacity=(sum_E2/StepPerBlock-energy*energy)/(T*T);
    double magnetic_susceptibility=(sum_M2/StepPerBlock-magnetization*magnetization)/T;

    ofstream en_output, magn_output, susc_output, spec_heat_output;
    ofstream corr_output, corr_time_output;
    ofstream corr_avg_output, avg_corr_time_output;
    en_output.open("internal_energy.dat", ios::app);
    magn_output.open("magnetization.dat", ios::app);
    corr_output.open("correlation_function.dat", ios::app);
    corr_time_output.open("correlation_time.dat", ios::app);
    susc_output.open("magnetic_susceptibility.dat", ios::app);
    spec_heat_output.open("heat_capacity.dat", ios::app);
    corr_avg_output.open("average_correlation_function.dat", ios::trunc);
    avg_corr_time_output.open("average_correlation_time.dat", ios::trunc);

    block_norm++;

    for(int r=0; r<R_CORR; r++){
        corr_func[r]=sum_sisj[r]/StepPerBlock - magnetization*magnetization;
        avg_corr_func[r]+=corr_func[r];

        corr_output<<setw(wd)<<r+1<<setw(wd)<<corr_func[r]<<endl;
        corr_avg_output<<setw(wd)<<r+1<<setw(wd)<<avg_corr_func[r]/block_norm<<endl;
    }

    for(int t=0; t<T_CORR; t++){
        avg_t_corr_func[t]+=t_corr_func[t];
        t_corr_func[t]=t_corr_func[t]/StepPerBlock;

        corr_time_output<<setw(wd)<<t+1<<setw(wd)<<t_corr_func[t]<<endl;
        avg_corr_time_output<<setw(wd)<<t+1<<setw(wd)<<avg_t_corr_func[t]/(StepPerBlock*(block_norm-1))<<endl;
    }

    global_E+=energy;
    global_E2+=energy*energy;
    global_M+=magnetization;
    global_M2+=magnetization*magnetization;
    error_E=Error(global_E, global_E2, block_norm);
    error_M=Error(global_M, global_M2, block_norm);

    global_Cv+=heat_capacity;
    global_chi+=magnetic_susceptibility;
    global_Cv2+=heat_capacity*heat_capacity;
    global_chi2+=magnetic_susceptibility*magnetic_susceptibility;
    error_Cv=Error(global_Cv, global_Cv2, block_norm);
    error_chi=Error(global_chi, global_chi2, block_norm);

    en_output<<setw(wd)<<block_norm<<setw(wd)<<energy<<setw(wd)<<global_E/block_norm<<setw(wd)<<error_E<<endl;
    magn_output<<setw(wd)<<block_norm<<setw(wd)<<magnetization<<setw(wd)<<global_M/block_norm<<setw(wd)<<error_M<<endl;
    susc_output<<setw(wd)<<block_norm<<setw(wd)<<magnetic_susceptibility<<setw(wd)<<global_chi/block_norm<<setw(wd)<<error_chi<<endl;
    spec_heat_output<<setw(wd)<<block_norm<<setw(wd)<<heat_capacity<<setw(wd)<<global_Cv/block_norm<<setw(wd)<<error_Cv<<endl;

    en_output.close();
    magn_output.close();
    corr_output.close();
    corr_time_output.close();
    susc_output.close();
    spec_heat_output.close();
    corr_avg_output.close();
    avg_corr_time_output.close();

    end_time = clock();
    double elapsed_secs = double(end_time - begin_time) / CLOCKS_PER_SEC;

    cout << "Block number: " << block_norm << endl;
    cout << "Elapsed time: " << elapsed_secs << endl;
    if(N_METROPOLIS_STEPS_PER_BLOCK>0){
        cout<<"Accepted flips rate: "<<(float)accepted_flips/(N_METROPOLIS_STEPS_PER_BLOCK*N)<<endl;
        accepted_flips=0;
    }
    cout << "----------------------------" << endl << endl;

}

void PrintFinal(){

    int wd=12;

    cout<<"E:            "<<global_E/block_norm<<endl;
    cout<<"St. dev. E:   "<<error_E<<endl;
    cout<<"M:            "<<global_M/block_norm<<endl;
    cout<<"St. dev. M:   "<<error_M<<endl;
    cout<<"Cv:           "<<global_Cv/block_norm<<endl;
    cout<<"St. dev. Cv:  "<<error_Cv<<endl;
    cout<<"Chi:          "<<global_chi/block_norm<<endl;
    cout<<"St. dev. chi: "<<error_chi<<endl<<endl;

    hist->print_hist("cluster_histogram.dat");

    delete[] sisj;
    delete[] corr_func;
    delete[] avg_corr_func;
    delete[] sum_sisj;
    delete[] m_history;
    delete[] t_corr_func;
    delete[] avg_t_corr_func;
    delete S;
    delete hist;

}

//##############################################################################

double Error(double sum, double sum2, int iblk){
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Read_From_File(string filename){
    ifstream read;
    cout<<"Reading input file: "<<filename<<" ..."<<endl;
    read.open(filename, ios::in);
    if(read.fail()){
        cout<< "ERROR: input file not found! "<<  endl;
        exit(1);
    }
    string temp, word;

    while (!read.eof()){
        read>>word;
        if (word=="L") {
            read>> temp;
            L=stoi(temp);
        }
        if (word=="T") {
            read>> temp;
            T=stof(temp);
        }
        if (word=="J") {
            read>> temp;
            J=stof(temp);
        }
        if (word=="N_BLOCKS") {
            read>> temp;
            N_BLOCKS=stoi(temp);
        }
        if (word=="N_EQUILIBRATION_WOLFF_STEPS") {
            read>> temp;
            N_EQUILIBRATION_WOLFF_STEPS=stoi(temp);
        }
        if (word=="N_EQUILIBRATION_METROPOLIS_STEPS") {
            read>> temp;
            N_EQUILIBRATION_METROPOLIS_STEPS=stoi(temp);
        }
        if (word=="N_WOLFF_STEPS_PER_BLOCK") {
            read>> temp;
            N_WOLFF_STEPS_PER_BLOCK=stoi(temp);
        }
        if (word=="N_METROPOLIS_STEPS_PER_BLOCK") {
            read>> temp;
            N_METROPOLIS_STEPS_PER_BLOCK=stoi(temp);
        }
        if (word=="SEED") {
            read>> temp;
            IDUM=stoi(temp);
        }
        if (word=="INITIAL_CONFIG") {
            read>> temp;
            INITIAL_CONFIG=stoi(temp);
        }
        if (word=="R_CORR") {
            read>> temp;
            R_CORR=stoi(temp);
        }
        if (word=="T_CORR") {
            read>> temp;
            T_CORR=stoi(temp);
        }

    }

    read.close();

}
