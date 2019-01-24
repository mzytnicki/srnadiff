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

//' Normalize counts
//'
//' @param lengths          the sizes of the RLEs (one list per chromosome)
//' @param values           the values of the RLEs (one list per chromosome)
//' @param chromosomeSizes  the sizes of the chromosomes
//' @param librarySizes     number of elements per sample
//' @return                 the normalization factors
// [[Rcpp::export]]
NumericVector rcpp_normalization(ListOf < ListOf < IntegerVector > > &lengths,
                        ListOf < ListOf < IntegerVector > > &values,
                        IntegerVector &chromosomeSizes,
                        IntegerVector &librarySizes) {
    GenomeIterator iterator (lengths, values, chromosomeSizes);
    std::vector < std::pair < std::valarray < double >, int > > normValues;
    std::vector < std::pair < double, int > > normValuesPerSample;
    std::vector < std::pair < double, int > > avgValues;
    std::vector < std::pair < double, int > > avgNormValues;
    int nSamples     = lengths[0].size();
    //int nChromosomes = chromosomeSizes.size();
    std::valarray < double > theseValues (nSamples);
    //std::valarray < double > sums (nSamples);
    //std::valarray < double > normalizedSums (nSamples);
    //std::valarray < double > librarySizesArray (nSamples);
    NumericVector factors (nSamples);
    //std::valarray < double > factors (nSamples);
    /*
    for (int i = 0; i < nSamples; ++i) {
        librarySizesArray[i] = librarySizes[i];
    }
    */
    for (; ! iterator.isOver(); iterator.getNext()) {
        theseValues = iterator.getValuesDouble();
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
    /*
    librarySizesArray *= factors;
    for (iterator.reset(); ! iterator.isOver(); iterator.getNext()) {
        sums           += iterator.getValuesDouble();
        normalizedSums += (iterator.getValuesDouble() / librarySizesArray);
    }
    std::nth_element(std::begin(sums), std::begin(sums) + nSamples / 2,
                     std::end(sums));
    std::nth_element(std::begin(normalizedSums),
                     std::begin(normalizedSums) + nSamples / 2,
                     std::end(normalizedSums));
    double fs = sums[nSamples / 2] / normalizedSums[nSamples / 2];
    for (int chromosomeId = 0; chromosomeId < nChromosomes; ++chromosomeId) {
        for (int sampleId = 0; sampleId < nSamples; ++sampleId) {
            for (int i = 0; i < values[chromosomeId][sampleId].size(); ++i) {
                values[chromosomeId][sampleId][i] =
                    static_cast<int>(round(static_cast<double>(
                            values[chromosomeId][sampleId][i]) * fs));
            }
        }
    }
    */
}
