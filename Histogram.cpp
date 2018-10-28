#include "Histogram.hpp"

using namespace std;

Histogram::Histogram(double x_min, double x_max, int n_bin, bool discard_out_of_range){
    _hist=new int[n_bin];
    for(int i=0; i<n_bin; i++)
        _hist[i]=0;
    _bin_width=(x_max-x_min)/n_bin;
    _x_max=x_max;
    _x_min=x_min;
    _n_bin=n_bin;
    _hist_norm=0;
    _discard_out_of_range=discard_out_of_range;
};
Histogram::~Histogram(){
    delete[] _hist;
}
void Histogram::add_x(double x){
    int pos=(x-_x_min)/(_bin_width);
    if(_discard_out_of_range){
        if(x>=_x_min && x<=_x_max){
            _hist[pos]++;
            _hist_norm++;
        }
    }
    else{
        if(x<_x_min) x=_x_min;
        if(x>_x_max) x=_x_max;
        _hist[pos]++;
        _hist_norm++;
    }
};
void Histogram::print_hist(string histname){
    ofstream output;
    output.open(histname);
    double x_pos=_x_min+_bin_width/2.;
    for(int i=0; i<_n_bin; i++){
        output<<x_pos<<"  "<<(double)_hist[i]/_hist_norm<<endl;
        x_pos+=_bin_width;
    }
}
