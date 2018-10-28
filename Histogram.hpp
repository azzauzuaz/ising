#include <fstream>

#ifndef __Histogram__
#define __Histogram__

class Histogram{
public:
    Histogram(double x_min, double x_max, int n_bin, bool discard_out_of_range);
    ~Histogram();
    void add_x(double x);
    void print_hist(std::string);

private:
    int* _hist;
    double _bin_width;
    double _x_max;
    double _x_min;
    int _n_bin;
    int _hist_norm;
    bool _discard_out_of_range;
};

#endif
