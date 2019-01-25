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
        ListOf < ListOf < IntegerVector > > &lengths;
        ListOf < ListOf < IntegerVector > > &values;
        IntegerVector &chromosomeSizes;
        NumericVector normalizationFactors;
        int nSamples;
        int nChromosomes;
        std::valarray < int >    indices;
        std::valarray < int >    remainings;
        std::valarray < int >    theseValues;
        std::vector < int >      theseValuesVector;
        std::valarray < double > theseValuesDouble;
        std::valarray < int >    theseRawValues;
        std::vector < int >      theseRawValuesVector;
        std::valarray < double > theseRawValuesDouble;
        int step;
        long pos;
        int chromosomeId;
        bool chromosomeChange;
        bool over;

    public:
        GenomeIterator (ListOf < ListOf < IntegerVector > > &l,
                ListOf < ListOf < IntegerVector > > &v,
                IntegerVector &cs, NumericVector f):
                    lengths(l),
                    values(v),
                    chromosomeSizes(cs),
                    normalizationFactors(f),
                    nSamples(lengths[0].size()),
                    nChromosomes(chromosomeSizes.size()),
                    indices(nSamples),
                    remainings(nSamples),
                    theseValues(nSamples),
                    theseValuesVector(nSamples),
                    theseValuesDouble(nSamples),
                    theseRawValues(nSamples),
                    theseRawValuesVector(nSamples),
                    theseRawValuesDouble(nSamples),
                    pos(0),
                    chromosomeId(0),
                    chromosomeChange(false),
                    over(false) {
                reset();
        }

        GenomeIterator (ListOf < ListOf < IntegerVector > > &l,
                ListOf < ListOf < IntegerVector > > &v, IntegerVector &cs):
                    GenomeIterator(l, v, cs, NumericVector(l[0].size(), 1.0)) {}

        void reset (bool nextChromosome = false) {
            if (nextChromosome) {
                ++chromosomeId;
                chromosomeChange = true;
                if (chromosomeId == nChromosomes) {
                    over = true;
                    return;
                }
            }
            else {
                chromosomeId = 0;
            }
            for (int sampleId = 0; sampleId < nSamples; ++sampleId) {
                indices[sampleId]    = 0;
                remainings[sampleId] = lengths[chromosomeId][sampleId][0];
                theseRawValuesDouble[sampleId] =
                    theseRawValues[sampleId] =
                    theseRawValuesVector[sampleId] =
                        values[chromosomeId][sampleId][0];
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
            return as < std::string >(
                    as< CharacterVector >(chromosomeSizes.names())
                        [chromosomeId]);
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
            if (pos >= chromosomeSizes[chromosomeId]) {
                reset(true);
                return;
            }
            for (int sampleId = 0; sampleId < nSamples; ++sampleId) {
                remainings[sampleId] -= s;
                if (remainings[sampleId] == 0) {
                    ++indices[sampleId];
                    remainings[sampleId] =
                        lengths[chromosomeId][sampleId][indices[sampleId]];
                    theseRawValuesDouble[sampleId] =
                        theseRawValues[sampleId] =
                        theseRawValuesVector[sampleId] =
                            values[chromosomeId][sampleId][indices[sampleId]];
                    theseValuesDouble[sampleId] =
                        theseRawValues[sampleId] *
                        normalizationFactors[sampleId];
                    theseValues[sampleId] =
                        theseValuesVector[sampleId] =
                            round(theseValuesDouble[sampleId]);
                }
            }
            step = remainings.min();
        }
};
