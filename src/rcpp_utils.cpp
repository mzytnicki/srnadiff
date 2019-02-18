#include <Rcpp.h>
#include <vector>
#include <valarray>
#include <tuple>
#include <string>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

class GenomeIterator {
    private:
        List         &coverages;
        StringVector  chromosomes;
        NumericVector normalizationFactors;
        int nSamples;
        std::vector < IntegerVector > lengths;
        std::vector < IntegerVector > values;
        std::vector < bool >     chromosomesOver;
        std::valarray < int >    indices;
        std::valarray < int >    remainings;
        std::valarray < int >    theseValues;
        std::vector < int >      theseValuesVector;
        std::valarray < double > theseValuesDouble;
        std::valarray < int >    theseRawValues;
        std::vector < int >      theseRawValuesVector;
        std::valarray < double > theseRawValuesDouble;
        std::vector   < int >    chromosomeSizes;
        int step;
        long pos;
        int chromosomeId;
        bool chromosomeChange;
        bool over;

    public:
        GenomeIterator (List &c, NumericVector f):
                    coverages(c),
                    normalizationFactors(f),
                    nSamples(coverages.size()),
                    lengths(nSamples),
                    values(nSamples),
                    chromosomesOver(nSamples, false),
                    indices(nSamples),
                    remainings(nSamples),
                    theseValues(nSamples),
                    theseValuesVector(nSamples),
                    theseValuesDouble(nSamples),
                    theseRawValues(nSamples),
                    theseRawValuesVector(nSamples),
                    theseRawValuesDouble(nSamples),
                    chromosomeSizes(nSamples, 0),
                    pos(0),
                    chromosomeId(0),
                    chromosomeChange(false),
                    over(false) {
            if (coverages.size() != 0) {
                if (as<List>(coverages[0]).size() != 0) {
                    chromosomes = as<List>(as<S4>(coverages[0]).slot("listData")).names();
                }
                for (int sampleId = 1; sampleId < nSamples; ++sampleId) {
                    if (as<StringVector>(as<List>(as<S4>(coverages[0]).slot("listData")).names()).size() != chromosomes.size()) {
                        stop("Number of chromosomes differ between samples!");
                    }
                    for (int chromosomeId = 0; chromosomeId < chromosomes.size(); ++ chromosomeId) {
                        if (as<StringVector>(as<List>(as<S4>(coverages[0]).slot("listData")).names())[chromosomeId] != chromosomes[chromosomeId]) {
                            stop("Chromosomes differ between samples!");
                        }
                    }
                }
            }
            reset();
        }

        GenomeIterator (List &c):
                    GenomeIterator(c, NumericVector(as<List>(c[0]).size(), 1.0)) {}

        void reset (bool nextChromosome = false) {
            std::fill(chromosomesOver.begin(), chromosomesOver.end(), false);
            if (nextChromosome) {
                ++chromosomeId;
                chromosomeChange = true;
                if (chromosomeId == chromosomes.size()) {
                    over = true;
                    return;
                }
            }
            else {
                chromosomeId = 0;
            }
            for (int sampleId = 0; sampleId < nSamples; ++sampleId) {
                S4 rle = as<S4>(as<List>(as<S4>(coverages[sampleId]).slot("listData"))[chromosomeId]);
                lengths[sampleId]    = rle.slot("lengths");
                values[sampleId]     = rle.slot("values");
                indices[sampleId]    = 0;
                remainings[sampleId] = lengths[sampleId][0];
                theseRawValuesDouble[sampleId] =
                    theseRawValues[sampleId] =
                    theseRawValuesVector[sampleId] =
                        values[sampleId][0];
                theseValuesDouble[sampleId] = theseRawValues[sampleId] *
                                                normalizationFactors[sampleId];
                theseValues[sampleId] =
                    theseValuesVector[sampleId] =
                        round(theseValuesDouble[sampleId]);
            }
            step = remainings.min();
            pos  = 0;
            over = false;
        }

        std::valarray < int > &getRawValues () {
            return theseRawValues;
        }

        std::valarray < double > &getRawValuesDouble () {
            return theseRawValuesDouble;
        }

        std::vector < int > &getRawValuesVector () {
            return theseRawValuesVector;
        }

        std::valarray < int > &getValues () {
            return theseValues;
        }

        std::valarray < double > &getValuesDouble () {
            return theseValuesDouble;
        }

        std::vector < int > &getValuesVector () {
            return theseValuesVector;
        }

        int getStep () {
            return step;
        }

        bool isOver () {
            return over;
        }

        bool hasChangedChromosome () {
            return chromosomeChange;
        }

        int getChromosomeId () {
            return chromosomeId;
        }

        std::string getChromosome () {
            return as < std::string >(chromosomes[chromosomeId]);
        }

        std::string getChromosome (int i) {
            return as < std::string >(chromosomes[i]);
        }

        int getChromosomeSize () {
            return chromosomeSizes[chromosomeId];
        }

        int getChromosomeSize (int i) {
            return chromosomeSizes[i];
        }

        long getPosition () {
            return pos+1;
        }

        void getNext (int s = 0) {
            chromosomeChange = false;
            if (s == 0) {
                s = step;
            }
            pos += s;
            for (int sampleId = 0; sampleId < nSamples; ++sampleId) {
                if (! chromosomesOver[sampleId]) {
                    remainings[sampleId] -= s;
                    if (remainings[sampleId] == 0) {
                        int value = 0;
                        ++indices[sampleId];
                        if (indices[sampleId] == lengths[sampleId].size()) {
                            chromosomesOver[sampleId] = true;
                            if (std::all_of(chromosomesOver.begin(), chromosomesOver.end(), [] (bool b) { return b; })) {
                                reset(true);
                                return;
                            }
                            value = 0;
                            remainings[sampleId] = std::numeric_limits<int>::max();
                        }
                        else {
                            value = values[sampleId][indices[sampleId]];
                            remainings[sampleId] =
                                lengths[sampleId][indices[sampleId]];
                        }
                        theseRawValuesDouble[sampleId] =
                            theseRawValues[sampleId] =
                            theseRawValuesVector[sampleId] = value;
                        theseValuesDouble[sampleId] =
                            value * normalizationFactors[sampleId];
                        theseValues[sampleId] =
                            theseValuesVector[sampleId] =
                                round(theseValuesDouble[sampleId]);
                    }
                }
            }
            step = remainings.min();
            chromosomeSizes[chromosomeId] = std::max<int>(chromosomeSizes[chromosomeId], pos);
        }
};
