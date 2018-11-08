#ifndef __Ising__
#define __Ising__

//input
Lattice* S;
int L;
double J;
int N_BLOCKS;
double T;
int N_WOLFF_STEPS_PER_BLOCK;
int N_METROPOLIS_STEPS_PER_BLOCK;
int IDUM;
int INITIAL_CONFIG;
int N_EQUILIBRATION_WOLFF_STEPS;
int N_EQUILIBRATION_METROPOLIS_STEPS;
int R_CORR;
int T_CORR;
int T_CORR_STEP;
bool ABS_MAGN;

int N;

Histogram* hist;

int StepPerBlock;

double* corr_func;
double* avg_corr_func;
double* avg_corr_func2;
double* sisj;
double* sum_sisj;

int t_count=0;
double* avg_t_corr_func;
double* avg_t_corr_func2;
double* t_corr_func;
double* m_history;

//Measure
double _E;
double _M;
double _M2;
double _E2;
double sum_E=0, sum_M=0, sum_E2=0, sum_M2=0;

//Averages
int block_norm=0;
int accepted_flips=0;
int cluster_dim;

double error_E=0, error_M=0;
double error_Cv=0, error_chi=0;
double global_E=0, global_M=0, global_M2=0, global_E2=0;
double global_Cv=0, global_Cv2=0, global_chi=0, global_chi2=0;

//time measurement

clock_t begin_time;
clock_t end_time;

//functions
void Input();
void Reset();
void Measure();
void Accumulate();
void Average();
double Error(double, double, int);
void PrintFinal();
void Read_From_File(std::string);

#endif
