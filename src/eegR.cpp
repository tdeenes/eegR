#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector consectrueRcpp(LogicalMatrix x) {
    int nrow = x.nrow(), ncol = x.ncol();
    bool prev;
    int colmax;
    int runlength;
    IntegerVector lengths(ncol);
    for(int ci = 0; ci < ncol; ++ci) {
        colmax = 0;
        runlength = 0;
        if (x(1,ci)) {
            prev = true;
            runlength++;
            colmax++;
        } else {
            prev = false;
        }
        for(int ri = 1; ri < nrow; ++ri) {
            if (x(ri,ci)) {
                runlength++; 
                prev = true;
            } else if (prev) {
                if (runlength>colmax) {
                    colmax = runlength;
                    runlength = 0;
                }
                prev = false;
            }
        }
        if (runlength>colmax) {
            colmax = runlength;
        }
        lengths[ci] = colmax;
    }
    return lengths;
}

// [[Rcpp::export]]
List rleRcpp(NumericMatrix x) {
    int nrow = x.nrow(), ncol = x.ncol();
    std::vector<int> lengths;
    std::vector<double> values;
    std::vector<int> matcol;
    
    // Initialise first value
    int i = 0;
    double prev;
    
    for(int ci = 0; ci < ncol; ++ci) {
        prev = x(0,ci);
        values.push_back(prev);
        lengths.push_back(1);
        matcol.push_back(ci+1);
        for(int ri = 1; ri < nrow; ++ri) {
            if (prev == x(ri,ci)) {
                lengths[i]++;
            } else {
                prev = x(ri,ci);
                values.push_back(prev);
                lengths.push_back(1);
                matcol.push_back(ci+1);
                i++;
            }
        }
        i++;
    }
    return List::create(_["lengths"] = lengths, _["values"] = values, _["matrixcolumn"] = matcol);
}
