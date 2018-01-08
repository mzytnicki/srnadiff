#include <Rcpp.h>
#include <vector>
#include <map>
#include <algorithm>
#include <valarray>
#include "rcpp_utils.cpp"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

bool addElement(unsigned int start, unsigned int end, std::string &chromosome,
                std::vector < unsigned int > &starts,
                std::vector < unsigned int > &ends,
                std::vector < std::string > &chromosomes,
                bool first, int distance, int size) {
		if (end - start + 1 < size) {
				return false;
		}
		if ((! first) && (start - ends.back() <= distance)) {
				ends.back() = end;
		}
		else {
				starts.push_back(start);
				ends.push_back(end);
				chromosomes.push_back(chromosome);
		}
		return true;
}

//' Compute naive method.
//'
//' @param lengths          the sizes of the RLEs (one list per chromosome)
//' @param values           the values of the RLEs (one list per chromosome)
//' @param chromosomeSizes  the sizes of the chromosomes
//' @param depth            minimum number of reads per position
//' @param distance         threshold to merge consecutive regions
//' @param size             minimum region size
//' @return                 the unique counts
// [[Rcpp::export]]
DataFrame rcpp_naive(ListOf < ListOf < IntegerVector > > &lengths,
                     ListOf < ListOf < IntegerVector > > &values,
                     IntegerVector &chromosomeSizes, int depth,
                     int distance, int size) {
    std::vector < std::vector < unsigned int > >::iterator it;
    std::vector < unsigned int > starts, ends;
    std::vector < std::string > chromosomes;
    unsigned int nSamples = lengths[0].size();
    bool         inRegion = false;
    bool         first    = true;
    unsigned int start;
    depth *= nSamples;
    for (GenomeIterator iterator (lengths, values, chromosomeSizes); ;
         iterator.getNext()) {
        if (iterator.hasChangedChromosome() || iterator.isOver()) {
            if (inRegion) {
                unsigned int previousChromosomeId =
                    iterator.getChromosomeId() - 1;
                std::string chromosome = as < std::string >(
                    as< CharacterVector >(chromosomeSizes.names())
                    [previousChromosomeId]);
                addElement(start, chromosomeSizes[previousChromosomeId],
                           chromosome, starts, ends, chromosomes, first,
                           distance, size);
            }
            inRegion = false;
            first    = true;
            if (iterator.isOver()) {
                return DataFrame::create(_["seqnames"] = chromosomes,
                                         _["start"]    = starts,
                                         _["end"]      = ends);
            }
        }
        if (iterator.getValues().sum() >= depth) {
            if (not inRegion) {
                inRegion = true;
                start    = iterator.getPosition();
            }
        }
        else if (inRegion) {
            std::string chromosome = iterator.getChromosome();
            if (addElement(start, iterator.getPosition(), chromosome, starts,
                           ends, chromosomes, first, distance, size)) {
                first = false;
            }
            inRegion = false;
        }
    }
}
