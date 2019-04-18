#include <Rcpp.h>
#include <vector>
#include <map>
#include <algorithm>
#include <valarray>
#include "rcpp_utils.cpp"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix rcpp_buildHmm(List &coverages, int minDepth) {
    std::vector < std::vector < int > > outputValues;
    std::vector < std::vector < int > >::iterator it;
    for (GenomeIterator iterator (coverages); ;
         iterator.getNext()) {
        if (iterator.hasChangedChromosome() || iterator.isOver()) {
            sort(outputValues.begin(), outputValues.end());
            auto it = std::unique(outputValues.begin(), outputValues.end());
            outputValues.resize(std::distance(outputValues.begin(), it));
            if (iterator.isOver()) {
                IntegerMatrix matrix(outputValues.size(), coverages.size());
                for (size_t i = 0; i < outputValues.size(); ++i) {
                    matrix.row(i) = IntegerVector(outputValues[i].begin(),
                               outputValues[i].end());
                }
                return matrix;
            }
        }
        if (iterator.getValues().max() >= minDepth) {
            outputValues.push_back(iterator.getRawValuesVector());
        }
    }
}


// [[Rcpp::export]]
DataFrame rcpp_viterbi(List          &coverages,
                       NumericMatrix &transitions,
                       NumericMatrix &emissions,
                       float         emissionThreshold,
                       NumericVector &starts,
                       IntegerVector &counts,
                       NumericVector &pvalues,
                       int minDepth, int minSize, int maxSize) {
    static const int NO_DIFF_CLASS = 0;
    static const int DIFF_CLASS    = 1;
    static const int N_CLASSES     = 2;
    int nSamples                   = coverages.size();
    std::vector < double > startsArray (starts.begin(), starts.end());
    std::vector < double > previousP = startsArray;
    std::vector < double > currentP(2);
    std::vector < std::vector < int > > previousStates;
    std::vector < int > statePos;
    std::map < std::vector < int >, double > pvalueMap;
    std::vector < int > allDiffStarts, allDiffEnds;
    std::vector < std::string > allDiffChromosomes;
    bool valueChange = true;
    int step = 0;
    double pvalue;
    for (int i = 0; i < pvalues.size(); ++i) {
        std::vector < int > row(nSamples);
        for (int sampleId = 0; sampleId < nSamples; ++sampleId) {
            row[sampleId] = counts(i, sampleId);
        }
        if (! isnan(pvalues[i])) {
            pvalueMap[row] = pvalues[i];
        }
    }
    for (GenomeIterator iterator (coverages); ; iterator.getNext(step)) {
        if (iterator.hasChangedChromosome() || iterator.isOver()) {
            int chromosomeId     = iterator.getChromosomeId()-1;
            int currentState     = (currentP[NO_DIFF_CLASS] <=
                                             currentP[DIFF_CLASS])?
                                             NO_DIFF_CLASS: DIFF_CLASS;
            int pos, previousPos = iterator.getChromosomeSize(chromosomeId)+1;
            bool inDiff = (currentState == DIFF_CLASS);
            std::vector < int > diffStarts, diffEnds;
            std::string chromosome = iterator.getChromosome(chromosomeId);
            if (inDiff) {
                diffEnds.push_back(iterator.getChromosomeSize(chromosomeId));
            }
            for (int i = previousStates.size(); i > 0; --i) {
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
            int cpt = 0;
            for (size_t i = 0; i < diffStarts.size(); ++i) {
                int size = diffEnds[i] - diffStarts[i] + 1;
                if ((minSize <= size) && (size <= maxSize)) {
                    allDiffStarts.push_back(diffStarts[i]-1);
                    allDiffEnds.push_back(diffEnds[i]-1);
                    allDiffChromosomes.push_back(chromosome);
                    ++cpt;
                }
            }
            if (iterator.isOver()) {
                return DataFrame::create(_["seqnames"] = allDiffChromosomes,
                                         _["start"] = allDiffStarts,
                                         _["end"] = allDiffEnds);
            }
            valueChange = true;
            previousP   = startsArray;
            previousStates.clear();
            statePos.clear();
        }
        std::vector < int > previousState(2);
        if (valueChange) {
            pvalue = 1.0;
            if (iterator.getValues().max() >= minDepth) {
                auto pvalueIt = pvalueMap.find(iterator.getRawValuesVector());
                if (pvalueIt != pvalueMap.end()) {
                    pvalue = pvalueMap[iterator.getRawValuesVector()];
                }
            }
        }
        currentP[NO_DIFF_CLASS] = currentP[DIFF_CLASS] = INFINITY;
        for (int pc = 0; pc < N_CLASSES; ++pc) {
            for (int nc = 0; nc < N_CLASSES; ++nc) {
                int pvalueIndex = (nc == NO_DIFF_CLASS)?
                    ((pvalue <= emissionThreshold)? 0: 1):
                     ((pvalue <= emissionThreshold)? 0: 1);
                double p = emissions(nc, pvalueIndex) + transitions(pc, nc) +
                    previousP[pc];
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
        else if (previousState[NO_DIFF_CLASS] == NO_DIFF_CLASS &&
                 previousState[DIFF_CLASS] == NO_DIFF_CLASS &&
                 pvalue == 1.0 &&
                 currentP[NO_DIFF_CLASS] < currentP[DIFF_CLASS]) {
            step = 0;
        }
    }
}
