#include <Rcpp.h>
#include <vector>
#include <map>
#include <algorithm>
#include <valarray>
#include "rcpp_utils.cpp"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' Compute unique counts.
//'
//' @param lengths          the sizes of the RLEs (one list per chromosome)
//' @param values           the values of the RLEs (one list per chromosome)
//' @param chromosomeSizes  the sizes of the chromosomes
//' @param minDepth         the minimum read coverage
//' @return                 the unique counts
// [[Rcpp::export]]
IntegerMatrix rcpp_buildHmm(ListOf < ListOf < IntegerVector > > &lengths,
                   ListOf < ListOf < IntegerVector > > &values,
                   IntegerVector &chromosomeSizes, int minDepth) {
    std::vector < std::vector < unsigned int > > outputValues;
    std::vector < std::vector < unsigned int > >::iterator it;
    for (GenomeIterator iterator (lengths, values, chromosomeSizes); ; iterator.getNext()) {
                if (iterator.hasChangedChromosome() || iterator.isOver()) {
                        sort(outputValues.begin(), outputValues.end());
                        auto it = std::unique(outputValues.begin(), outputValues.end());
                        outputValues.resize(std::distance(outputValues.begin(), it));
                        if (iterator.isOver()) {
                                IntegerMatrix matrix(outputValues.size(), lengths[0].size());
                                for (size_t i = 0; i < outputValues.size(); ++i) {
                                        matrix.row(i) = IntegerVector(outputValues[i].begin(), outputValues[i].end());
                                }
                                return matrix;
                        }
                }
                if (iterator.getValues().max() >= minDepth) {
                        outputValues.push_back(iterator.getValuesVector());
                }
    }
}

//' Run the Viterbi algorithm on the HMM.
//'
//' @param chromosomeSizes   the sizes of the chromosomes
//' @param transitions       the transition log-probabilities
//' @param emissions         the emission log-probabilities
//' @param emissionThreshold the emission threshold
//' @param starts            the start log-probabilities
//' @param counts            the unique counts
//' @param pvalues           the p-values of the counts
//' @param lengths           the sizes of the RLEs (one list per chromosome)
//' @param values            the values of the RLEs (one list per chromosome)
//' @param minDepth          the minimum read coverage
//' @param minSize           the minimum size region
//' @param maxSize           the maximum size region
//' @return                  a segmentation of the chromosomes
// [[Rcpp::export]]
DataFrame rcpp_viterbi(IntegerVector &chromosomeSizes,
                       NumericMatrix &transitions,
                       NumericMatrix &emissions,
                       float         emissionThreshold,
                       NumericVector &starts,
                       IntegerVector &counts,
                       NumericVector &pvalues,
                       ListOf < ListOf < IntegerVector > > &lengths,
                       ListOf < ListOf < IntegerVector > > &values, int minDepth, int minSize, int maxSize) {
    static const unsigned int NO_DIFF_CLASS = 0;
    static const unsigned int DIFF_CLASS    = 1;
    static const unsigned int N_CLASSES     = 2;
    size_t nSamples     = lengths[0].size();
    size_t nChromosomes = chromosomeSizes.size();
    std::vector < double > startsArray (starts.begin(), starts.end());
    std::vector < double > previousP = startsArray;
    std::vector < double > currentP(2);
    std::vector < std::vector < unsigned int > > previousStates;
    std::vector < unsigned int > statePos;
    std::map < std::vector < unsigned int >, double > pvalueMap;
    std::vector < unsigned int > allDiffStarts, allDiffEnds;
    std::vector < std::string > allDiffChromosomes;
    bool valueChange = true;
    unsigned int step = 0;
    double pvalue;
    for (size_t i = 0; i < pvalues.size(); ++i) {
        std::vector < unsigned int > row(nSamples);
        for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
            row[sampleId] = counts(i, sampleId);
        }
        pvalueMap[row] = pvalues[i];
    }
    for (GenomeIterator iterator (lengths, values, chromosomeSizes); ; iterator.getNext(step)) {
        if (iterator.hasChangedChromosome() || iterator.isOver()) {
            unsigned int chromosomeId   = iterator.getChromosomeId()-1;
            unsigned int chromosomeSize = chromosomeSizes[chromosomeId];
            unsigned int currentState = (currentP[NO_DIFF_CLASS] <= currentP[DIFF_CLASS])? NO_DIFF_CLASS: DIFF_CLASS;
            unsigned int pos, previousPos = chromosomeSize+1;
            bool inDiff = (currentState == DIFF_CLASS);
            std::vector < unsigned int > diffStarts, diffEnds;
            std::string chromosome = as < std::string >(as< CharacterVector >(chromosomeSizes.names()) [chromosomeId]);
            if (inDiff) {
                diffEnds.push_back(chromosomeSize);
            }
            for (size_t i = previousStates.size(); i > 0; --i) {
                pos = statePos[i-1];
                if (previousPos != pos+1) {
                    if (inDiff) {
                        diffStarts.push_back(previousPos);
                    }
                    inDiff       = false;
                    currentState = NO_DIFF_CLASS;
                }
                currentState = previousStates[i-1][currentState];
                if ((currentState == DIFF_CLASS) && (! inDiff)) {
                    diffEnds.push_back(pos);
                    inDiff = true;
                }
                else if ((currentState == NO_DIFF_CLASS) && (inDiff)) {
                    diffStarts.push_back(previousPos);
                    inDiff = false;
                }
                previousPos = pos;
            }
            if (inDiff) {
                diffStarts.push_back(previousPos);
            }
            std::reverse(diffStarts.begin(), diffStarts.end());
            std::reverse(diffEnds.begin(), diffEnds.end());
						unsigned int cpt = 0;
            for (size_t i = 0; i < diffStarts.size(); ++i) {
                unsigned int size = diffEnds[i] - diffStarts[i] + 1;
                if ((minSize <= size) && (size <= maxSize)) {
                    allDiffStarts.push_back(diffStarts[i]);
                    allDiffEnds.push_back(diffEnds[i]);
                    allDiffChromosomes.push_back(chromosome);
										++cpt;
                }
            }
            if (iterator.isOver()) {
                return DataFrame::create(_["seqnames"] = allDiffChromosomes,
                                                                     _["start"] = allDiffStarts ,
                                                                     _["end"] = allDiffEnds);
            }
            valueChange = true;
            previousP   = startsArray;
            previousStates.clear();
            statePos.clear();
        }
        std::vector < unsigned int > previousState(2);
        if (valueChange) {
            if (iterator.getValues().max() >= minDepth) {
                pvalue = pvalueMap[iterator.getValuesVector()];
            }
            else {
                pvalue = 1.0;
            }
        }
        currentP[NO_DIFF_CLASS] = currentP[DIFF_CLASS] = INFINITY;
        for (unsigned int pc = 0; pc < N_CLASSES; ++pc) {
            for (unsigned int nc = 0; nc < N_CLASSES; ++nc) {
                unsigned int pvalueIndex = (nc == NO_DIFF_CLASS)? ((pvalue <= emissionThreshold)? 0: 1): ((pvalue <= emissionThreshold)? 0: 1);
                double                 p = emissions(nc, pvalueIndex) + transitions(pc, nc) + previousP[pc];
                if (p < currentP[nc]) {
                    previousState[nc] = pc;
                    currentP[nc]      = p;
                }
            }
        }
        previousP = currentP;
        step = 1;
        if (previousState[NO_DIFF_CLASS] != NO_DIFF_CLASS ||
            previousState[DIFF_CLASS] != NO_DIFF_CLASS) {
            statePos.push_back(iterator.getPosition());
            previousStates.push_back(previousState);
        }
        else if (previousState[NO_DIFF_CLASS] == NO_DIFF_CLASS && previousState[DIFF_CLASS] == NO_DIFF_CLASS && pvalue == 1.0 && currentP[NO_DIFF_CLASS] < currentP[DIFF_CLASS]) {
            step = 0;
        }
    }
}
