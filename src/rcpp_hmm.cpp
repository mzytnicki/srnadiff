#include <Rcpp.h>
#include <vector>
#include <map>
#include <algorithm>
using namespace Rcpp;

//' Compute unique counts.
//'
//' @param lengths          the sizes of the RLEs (one list per chromosome)
//' @param values           the values of the RLEs (one list per chromosome)
//' @param chromosomeSizes  the sizes of the chromosomes
//' @return                 the unique counts
// [[Rcpp::export]]
List rcpp_buildHmm(ListOf < ListOf < IntegerVector > > &lengths,
                   ListOf < ListOf < IntegerVector > > &values,
                   IntegerVector &chromosomeSizes) {
    std::vector < std::vector < unsigned int > > outputValues;
    size_t nSamples     = lengths[0].size();
    size_t nChromosomes = chromosomeSizes.size();
    for (size_t chromosomeId = 0; chromosomeId < nChromosomes; ++chromosomeId) {
        std::vector < size_t > indices(nSamples, 0);
        std::vector < unsigned int > remainings(nSamples);
        for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
            remainings[sampleId] = lengths[chromosomeId][sampleId][0];
        }
        for (unsigned int pos = 0; pos < chromosomeSizes[chromosomeId];) {
            std::vector < unsigned int > theseValues (nSamples);
            for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
                theseValues[sampleId] =
                    values[chromosomeId][sampleId][indices[sampleId]];
            }
            unsigned int sumValue = 0, maxValue;
            for (unsigned int sampleId = 0; sampleId < nSamples; ++sampleId) {
                sumValue += theseValues[sampleId];
            }
            maxValue = *std::max_element(theseValues.begin(), theseValues.end());
            if (maxValue >= 5 && sumValue >= 10) {
                outputValues.push_back(theseValues);
            }
            unsigned int step = *std::min_element(remainings.begin(),
                                                  remainings.end());
            pos              += step;
            for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
                remainings[sampleId] -= step;
                if (remainings[sampleId] == 0) {
                    ++indices[sampleId];
                    remainings[sampleId] =
                        lengths[chromosomeId][sampleId][indices[sampleId]];
                }
            }
        }
        std::vector < std::vector < unsigned int> >::iterator it =
            std::unique(outputValues.begin(), outputValues.end());
        outputValues.resize(std::distance(outputValues.begin(),it));
    }
    return wrap(outputValues);
}

//' Run the Viterbi algorithm on the HMM.
//'
//' @param chromosomeSizes  the sizes of the chromosomes
//' @param transitions      the transition log-probabilities
//' @param emissions        the emission log-probabilities
//' @param starts           the start log-probabilities
//' @param counts           the unique counts
//' @param pvalues          the p-values of the counts
//' @param lengths          the sizes of the RLEs (one list per chromosome)
//' @param values           the values of the RLEs (one list per chromosome)
//' @return                 a segmentation of the chromosomes
// [[Rcpp::export]]
DataFrame rcpp_viterbi(IntegerVector &chromosomeSizes,
                       NumericMatrix &transitions,
                       NumericMatrix &emissions,
                       NumericVector &starts,
                       IntegerVector &counts,
                       NumericVector &pvalues,
                       ListOf < ListOf < IntegerVector > > &lengths,
                       ListOf < ListOf < IntegerVector > > &values) {
    static const unsigned int NO_DIFF_CLASS = 0;
    static const unsigned int DIFF_CLASS    = 1;
    static const unsigned int N_CLASSES     = 2;
    size_t nSamples     = lengths[0].size();
    size_t nChromosomes = chromosomeSizes.size();
    std::map < std::vector < unsigned int >, double > pvalueMap;
    std::vector < unsigned int > allDiffStarts, allDiffEnds;
    std::vector < std::string > allDiffChromosomes;
    for (size_t i = 0; i < pvalues.size(); ++i) {
        std::vector < unsigned int > row(nSamples);
        for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
            row[sampleId] = counts(i, sampleId);
        }
        pvalueMap[row] = pvalues[i];
    }
    for (size_t chromosomeId = 0; chromosomeId < nChromosomes; ++chromosomeId) {
        //Rcout << "chr " << as<CharacterVector>(chromosomeSizes.names())[chromosomeId] << std::endl;
        std::vector < size_t > indices(nSamples, 0);
        std::vector < unsigned int > remainings(nSamples),
                currentCounts(nSamples);
        bool valueChange = true;
        unsigned int chromosomeSize = chromosomeSizes[chromosomeId];
        double pvalue;
        std::vector < double > previousP (starts.begin(), starts.end()),
            currentP(2);
        std::vector < std::vector < unsigned int > > previousStates;
        std::vector < unsigned int > statePos;
        for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
            remainings[sampleId] = lengths[chromosomeId][sampleId][0];
        }
        for (unsigned int pos = 0; pos < chromosomeSize;) {
            std::vector < unsigned int > previousState(2);
            if (valueChange) {
                unsigned int maxCount = 0, sumCount = 0;
                for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
                    unsigned int tmp        =
                        values[chromosomeId][sampleId][indices[sampleId]];
                    currentCounts[sampleId] = tmp;
                    sumCount               += tmp;
                    maxCount                =
                        std::max<unsigned int>(maxCount, tmp);
                    if (maxCount >= 5 && sumCount >= 10) {
                        pvalue = pvalueMap[currentCounts];
                    }
                    else {
                        pvalue = 1.0;
                    }
                }
            }
            currentP[NO_DIFF_CLASS] = currentP[DIFF_CLASS] = INFINITY;
            for (unsigned int pc = 0; pc < N_CLASSES; ++pc) {
                for (unsigned int nc = 0; nc < N_CLASSES; ++nc) {
                    unsigned int pvalueIndex = (nc == NO_DIFF_CLASS)?
                        ((pvalue <= 0.85)? 0: 1): ((pvalue <= 0.15)? 0: 1);
                    double                 p = emissions(nc, pvalueIndex) +
                        transitions(pc, nc) + previousP[pc];
                    if (p < currentP[nc]) {
                        previousState[nc] = pc;
                        currentP[nc]      = p;
                    }
                }
            }
            previousP = currentP;
            unsigned int step = 1;
            if (previousState[NO_DIFF_CLASS] != NO_DIFF_CLASS ||
                previousState[DIFF_CLASS] != NO_DIFF_CLASS) {
                statePos.push_back(pos+1);
                previousStates.push_back(previousState);
            }
            else if (previousState[NO_DIFF_CLASS] == NO_DIFF_CLASS &&
                     previousState[DIFF_CLASS] == NO_DIFF_CLASS &&
                     pvalue == 1.0 &&
                     currentP[NO_DIFF_CLASS] < currentP[DIFF_CLASS]) {
                step = *std::min_element(remainings.begin(), remainings.end());
            }
            for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
                remainings[sampleId] -= step;
                if (remainings[sampleId] == 0) {
                    ++indices[sampleId];
                    remainings[sampleId] =
                        lengths[chromosomeId][sampleId][indices[sampleId]];
                }
            }
            pos += step;
        }
        unsigned int currentState =
            (currentP[NO_DIFF_CLASS] <= currentP[DIFF_CLASS])?
            NO_DIFF_CLASS: DIFF_CLASS;
        unsigned int pos, previousPos = chromosomeSize+1;
        bool inDiff = (currentState == DIFF_CLASS);
        std::vector < unsigned int > diffStarts, diffEnds;
        if (inDiff) {
            diffEnds.push_back(chromosomeSize);
            //Rcout << "1 end" << std::endl;
        }
        for (size_t i = previousStates.size(); i > 0; --i) {
            pos = statePos[i-1];
            if (previousPos != pos+1) {
                if (inDiff) {
                    diffStarts.push_back(previousPos);
                    //Rcout << "start " << previousPos << std::endl;
                }
                inDiff       = false;
                currentState = NO_DIFF_CLASS;
            }
            currentState = previousStates[i-1][currentState];
            if ((currentState == DIFF_CLASS) && (! inDiff)) {
                //Rcout << "end " << previousPos << std::endl;
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
        allDiffStarts.insert(allDiffStarts.end(), diffStarts.begin(),
                             diffStarts.end());
        allDiffEnds.insert(allDiffEnds.end(), diffEnds.begin(), diffEnds.end());
        allDiffChromosomes.insert(allDiffChromosomes.end(), diffStarts.size(),
          as < std::string >(as< CharacterVector >(chromosomeSizes.names())
                                      [chromosomeId]));
    }
    return DataFrame::create(_["seqnames"] = allDiffChromosomes,
                             _["start"] = allDiffStarts ,
                             _["end"] = allDiffEnds);
}
