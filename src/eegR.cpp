#include <stdlib.h>
#include <Rcpp.h>
using namespace Rcpp;

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

int indexFn (const int& chdim, const int& dimc, const int& dimt, const int& ic, const int& it) {
    return (chdim == 0) ? ((it*dimc) + ic) : ((ic*dimt) + it);
}


// [[Rcpp::export]] 
NumericMatrix tfce(NumericMatrix input_data, int chan_dim, IntegerMatrix ChN, NumericVector EH, int num_steps, bool has_negative, bool has_positive) {
    double valToAdd;
    int ii, iT, iC, tiC, tiT;
    int maxt, mint, temp, growingInd, growingCur, ChCurr, idx;
    double thresh0, thresh, delta;
    int dims[2] = {input_data.nrow(), input_data.ncol()};
    int dimsChN[2] = {ChN.nrow(), ChN.ncol()};
    int numVoxels = dims[0] * dims[1];
    chan_dim = chan_dim - 1;
    int dimC = dims[chan_dim];
    int time_dim = 1 - chan_dim;
    int dimT = dims[time_dim];
    double E = EH[0];
    double H = EH[1];
    NumericMatrix outData(dims[0], dims[1]);
    // calculate steps (thresholds)
    double fmax = 0.0;
    for (ii = 0; ii < numVoxels; ++ii) {
        if (has_positive && input_data[ii] > fmax) fmax = input_data[ii];
        else if (has_negative && -input_data[ii] > fmax) fmax = -input_data[ii];
    }
    delta = fmax/num_steps;
    thresh0 = delta/2.0;
    
    std::vector<bool> flagUsed;
    std::vector<int> grow_t;
    std::vector<int> grow_c;
    
    // main part
    for (ii = 0; ii < num_steps; ii++) {
        thresh = thresh0 + (double)ii*delta;
		valToAdd = 0.0;
        flagUsed.resize(numVoxels);
       
        for (iT = 0; iT < dimT; ++iT) {
    		for (iC = 0; iC < dimC; ++iC) {
    			// temp is the current point in the grid
    			temp = indexFn(chan_dim, dimC, dimT, iC, iT);
    			// Check if this point has been seen (flagUsed) and whether its over the threshold
    			if (has_positive && !flagUsed[temp] && input_data[temp] >= thresh) {
    				// make the flagUsed so that algorithm doesn't visit this point again
    				flagUsed[temp] = true;
    				growingInd = 1;
    				growingCur = 0;

    				// Define the current coordinates to continue with
                    grow_c.reserve(numVoxels);
                    grow_t.reserve(numVoxels);
                    grow_c.push_back(iC);
                    grow_t.push_back(iT);
                    
    				while (growingCur < growingInd) {
                        //This just limits not to overrun borders <-- in our case the zero padding
    				    //And creates a 3x3 windows of scanning for valid points (one point around the current point)
    				    //E.g. In 2 dimensions... if the current point is 2,j then maxi = Min of dimension size or 2+2
    				    //thus maxi = 4... mini = max of 0 or 2-1... so 1...
    				    //so in the next set of loops ti will look at i=1,i=2,and i=3 (not 4 because its set to < maxi)
                        maxt = MIN(dimT, grow_t[growingCur] + 2);
    			        mint = MAX(0, grow_t[growingCur] - 1);
    					
    					// start of the smaller scanning window
    					for (tiT = mint; tiT < maxt; ++tiT) {
    						for (tiC = 0; tiC < dimsChN[1]; ++tiC) {
    							idx = (tiC*dimsChN[0]) + grow_c[growingCur];
    							ChCurr = ChN[idx];
    							if (ChCurr == 0) {
    								break;
    							}  							
    							ChCurr = ChCurr - 1;
    							
    							temp = indexFn(chan_dim, dimC, dimT, ChCurr, tiT);
    							if (!flagUsed[temp] && input_data[temp] >= thresh) {
    								flagUsed[temp] = true;
    								grow_c.push_back(ChCurr);
    								grow_t.push_back(tiT);
    								growingInd++;
    								//Here the growing index increases everytime we find a value above the threshhold
    								//This grows depending on the origin of our supporting weight which is the outer loop
    								//-#
    							}
    						}	
    					}
    				    // GrowingCur increases one and thus looks at the next point that was found in the previous loop
    				    growingCur++;
    				}
    				// Reset back to 0 so that when the next "while loop" runs it adds the value to each point found
                    growingCur = 0;
    				valToAdd = pow(growingInd, E) * pow(thresh, H) * delta;
                    
    				// Adds the valToAdd to the points that it found in that cluster
    				while (growingCur < growingInd) {
                        temp = indexFn(chan_dim, dimC, dimT, grow_c[growingCur], grow_t[growingCur]);
    				    outData[temp] += valToAdd;
    				    // GrowingCur adds one so that next coordinate is added until growingCurr is less than the number of points found
    				    growingCur++;
    				}
                    grow_c.clear();
                    grow_t.clear();
    			}
    			// Check the same for negative values
    			if (has_negative && !flagUsed[temp] && -input_data[temp] >= thresh) {
    			    flagUsed[temp] = true; 
    			    growingInd = 1; 
    			    growingCur = 0;
    			    grow_c.reserve(numVoxels);
    			    grow_t.reserve(numVoxels);
    			    grow_c.push_back(iC);
    			    grow_t.push_back(iT);
    			    while (growingCur < growingInd) {
    			        maxt = MIN(dimT, grow_t[growingCur] + 2);
    			        mint = MAX(0, grow_t[growingCur] - 1);
    			        for (tiT = mint; tiT < maxt; ++tiT) {
    			            for (tiC = 0; tiC < dimsChN[1]; ++tiC) {
    			                idx = (tiC*dimsChN[0]) + grow_c[growingCur];
    			                ChCurr = ChN[idx];
    			                if (ChCurr == 0) {
    			                    break;
    			                }  							
    			                ChCurr = ChCurr - 1;
    			                temp = indexFn(chan_dim, dimC, dimT, ChCurr, tiT);
    			                if (!flagUsed[temp] && -input_data[temp] >= thresh) {
    			                    flagUsed[temp] = true;
    			                    grow_c.push_back(ChCurr);
    			                    grow_t.push_back(tiT);
    			                    growingInd++;
    			                }
    			            }	
    			        }
    			        growingCur++;
    			    }
    			    growingCur = 0;
    			    valToAdd = pow(growingInd, E) * pow(thresh, H) * delta;
    			    while (growingCur < growingInd) {
    			        temp = indexFn(chan_dim, dimC, dimT, grow_c[growingCur], grow_t[growingCur]);
    			        outData[temp] -= valToAdd;
    			        growingCur++;
    			    }
    			    grow_c.clear();
    			    grow_t.clear();
    			}
    		}
    	}
        flagUsed.clear();
    }
    return outData;
}

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

// [[Rcpp::export]]
CharacterVector charmatCollapse(CharacterMatrix x, int along_dim) {
    int nrow = x.nrow(), ncol = x.ncol();
    if (along_dim == 1) {
        CharacterVector out = no_init(ncol);
        for( int j = 0; j < ncol; j++ ) {
            CharacterMatrix::Column tmp = x(_, j);
            out[j] = collapse( tmp );
        }
        return out;
    } else {
        CharacterVector out = no_init(nrow);
        for( int i = 0; i < nrow; i++ ) {
            CharacterMatrix::Row tmp = x(i, _);
            out[i] = collapse( tmp );
        }
        return out;
    }
}

// [[Rcpp::export]]
NumericMatrix groupsum(NumericMatrix x, IntegerVector g, int ug) {
    int nrow = x.nrow(), ncol = x.ncol();
    int gind;
    NumericMatrix out(ug, ncol);
    for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nrow; i++) {
            gind = g[i] - 1;
            out(gind, j) = out(gind, j) + x(i, j);
        }
    }
    return out;
}

//
// fast sweep ------
//

// operators
template <typename opType>
opType f_mtype_operator(opType x, opType y, std::string f) {
    if (f == "*") {
        return x * y;
    } else if (f == "+") {
        return x + y;
    } else if (f == "-") {
        return x - y;
    } else {
        throw std::invalid_argument( "operator not implemented" );
    }
}

double f_double_operator(double x, double y, std::string f) {
    if (f == "/") {
        return x / y;
    } else if (f == "^") {
        return pow(x, y);
    } else {
        throw std::invalid_argument( "operator not implemented" );
    }
}

template <typename logType>
bool f_logical_operator(logType x, logType y, std::string f) {
    if (f == "<") {
        return x < y;
    } else if (f == ">") {
        return x > y;
    } else if (f == "==") {
        return x == y;
    } else if (f == ">=") {
        return x >= y;
    } else if (f == "<=") {
        return x <= y;
    } else if (f == "!=") {
        return x != y;
    } else {
        throw std::invalid_argument( "operator not implemented" );
    }
}

// [[Rcpp::export]]
SEXP sweepcol_multitype_cpp(SEXP x, SEXP y, std::string fun){
    switch( TYPEOF(x) ){
        case REALSXP: {
            NumericMatrix X(x);
            NumericVector Y(y);
            int nrow = X.nrow();
            int ncol = X.ncol();
            NumericMatrix output(nrow, ncol);
            double yy;
            for(int i = 0; i < ncol; i++) {
                yy = Y[i];
                for(int j = 0; j < nrow; j++) {
                    output(j, i) = f_mtype_operator(X(j, i), yy, fun);
                }
            }
            return wrap( output );
        }
        case INTSXP: {
            IntegerMatrix X(x);
            IntegerVector Y(y);
            int nrow = X.nrow();
            int ncol = X.ncol();
            IntegerMatrix output(nrow, ncol);
            int yy;
            for(int i = 0; i < ncol; i++) {
                yy = Y[i];
                for(int j = 0; j < nrow; j++) {
                    output(j, i) = f_mtype_operator(X(j, i), yy, fun);
                }
            }
            return wrap( output );
        }
        default: {
            return R_NilValue;
        }
    }
}

// [[Rcpp::export]]
NumericMatrix sweepcol_double_cpp(NumericMatrix x, NumericVector y, std::string fun){
    NumericMatrix output(x.nrow(), x.ncol());
    double yy;
    for(int i = 0; i < x.ncol(); i++) {
        yy = y[i];
        for(int j = 0; j < x.nrow(); j++) {
            output(j, i) = f_double_operator(x(j, i), yy, fun);
        }
    }
    return output ;
} 

// [[Rcpp::export]]
LogicalMatrix sweepcol_logical_cpp(SEXP x, SEXP y, std::string fun){
    switch( TYPEOF(x) ){
        case REALSXP: {
            NumericMatrix X(x);
            NumericVector Y(y);
            int nrow = X.nrow();
            int ncol = X.ncol();
            LogicalMatrix output(nrow, ncol);
            double yy;
            for(int i = 0; i < ncol; i++) {
                yy = Y[i];
                for(int j = 0; j < nrow; j++) {
                    output(j, i) = f_logical_operator(X(j, i), yy, fun);
                }
            }
            return output;
        }
        case INTSXP: {
            IntegerMatrix X(x);
            IntegerVector Y(y);
            int nrow = X.nrow();
            int ncol = X.ncol();
            LogicalMatrix output(nrow, ncol);
            int yy;
            for(int i = 0; i < ncol; i++) {
                yy = Y[i];
                for(int j = 0; j < nrow; j++) {
                    output(j, i) = f_logical_operator(X(j, i), yy, fun);
                }
            }
            return output;
        }
        default: {
            return R_NilValue;
        }
    }
} 

