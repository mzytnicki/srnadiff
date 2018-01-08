#include <Rcpp.h>
#include <vector>
#include <valarray>
#include <tuple>
#include <string>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

class GenomeIterator {
	private:
		ListOf < ListOf < IntegerVector > > &lengths;
		ListOf < ListOf < IntegerVector > > &values;
		IntegerVector &chromosomeSizes;
		size_t nSamples;
		size_t nChromosomes;
		std::valarray < size_t > indices;
		std::valarray < unsigned int > remainings;
		std::valarray < unsigned int > theseValues;
		std::vector < unsigned int > theseValuesVector;
		std::valarray < double > theseValuesDouble;
		unsigned int step;
		unsigned long pos;
		size_t chromosomeId;
		bool chromosomeChange;
		bool over;

	public:
		GenomeIterator (ListOf < ListOf < IntegerVector > > &l,
				ListOf < ListOf < IntegerVector > > &v,
				IntegerVector &cs):
					lengths(l),
					values(v),
					chromosomeSizes(cs),
					nSamples(lengths[0].size()),
					nChromosomes(chromosomeSizes.size()),
					indices(nSamples),
					remainings(nSamples),
					theseValues(nSamples),
					theseValuesVector(nSamples),
					theseValuesDouble(nSamples),
					pos(0),
					chromosomeChange(false),
					over(false),
					chromosomeId(0) {
				reset();
		}

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
			for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
				indices[sampleId]           =0;
				remainings[sampleId]        =lengths[chromosomeId][sampleId][0];
				theseValues[sampleId]       =values[chromosomeId][sampleId][0];
				theseValuesVector[sampleId] =theseValues[sampleId];
				theseValuesDouble[sampleId] =theseValues[sampleId];
			}
			step = remainings.min();
			pos  = 0;
			over = false;
		}

		std::valarray < unsigned int > &getValues () {
			return theseValues;
		}

		std::valarray < double > &getValuesDouble () {
			return theseValuesDouble;
		}

		std::vector < unsigned int > &getValuesVector () {
			return theseValuesVector;
		}

		unsigned int getStep () {
			return step;
		}

		bool isOver () {
			return over;
		}

		bool hasChangedChromosome () {
			return chromosomeChange;
		}

		size_t getChromosomeId () {
			return chromosomeId;
		}

		std::string getChromosome () {
			return as < std::string >(
			        as< CharacterVector >(chromosomeSizes.names())
            		    [chromosomeId]);
		}

		unsigned long getPosition () {
			return pos+1;
		}

		void getNext (unsigned int s = 0) {
			chromosomeChange = false;
			if (s == 0) {
					s = step;
			}
			pos += s;
			if (pos >= chromosomeSizes[chromosomeId]) {
				reset(true);
				return;
			}
			for (size_t sampleId = 0; sampleId < nSamples; ++sampleId) {
				remainings[sampleId] -= s;
				if (remainings[sampleId] == 0) {
					++indices[sampleId];
					remainings[sampleId] =
					    lengths[chromosomeId][sampleId][indices[sampleId]];
					theseValues[sampleId] =
					    values[chromosomeId][sampleId][indices[sampleId]];
					theseValuesDouble[sampleId] = theseValues[sampleId];
					theseValuesVector[sampleId] = theseValues[sampleId];
				}
			}
			step = remainings.min();
		}
};
