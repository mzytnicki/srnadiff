#include <Rcpp.h>
#include <vector>
#include <valarray>
#include <algorithm>
#include <math.h>
#include "rcpp_utils.cpp"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

double computeMedian(std::vector < std::pair < double, int > > &table){
    sort(table.begin(), table.end());
    int i, s = 0, size = 0;
    for (auto &i: table) {
        size += i.second;
    }
    for (i = 0; s < size/2; ++i) {
        s += table[i].second;
    }
    return table[i].first;
}


// [[Rcpp::export]]
NumericVector rcpp_normalization(List &coverages, NumericVector &librarySizes) {
    GenomeIterator iterator (coverages);
    std::vector < std::pair < std::valarray < double >, int > > normValues;
    std::vector < std::pair < double, int > > normValuesPerSample;
    std::vector < std::pair < double, int > > avgValues;
    std::vector < std::pair < double, int > > avgNormValues;
    int nSamples     = coverages.size();
    std::valarray < double > theseValues (nSamples);
    NumericVector factors (nSamples);
    for (; ! iterator.isOver(); iterator.getNext()) {
        theseValues = iterator.getRawValuesDouble();
        if (theseValues.min() >= 1) {
            theseValues /= exp(log(theseValues).sum() / nSamples);
            normValues.emplace_back(theseValues, iterator.getStep());
        }
    }
    for (int sample = 0; sample < nSamples; ++sample) {
        normValuesPerSample.clear();
        normValuesPerSample.reserve(normValues.size());
        for (auto &normValue: normValues) {
            normValuesPerSample.emplace_back((normValue.first)[sample],
                                             normValue.second);
        }
        factors[sample] = computeMedian(normValuesPerSample) /
            librarySizes[sample];
    }
    factors = factors / (exp(sum(log(factors)) / nSamples));
    return factors;
}
