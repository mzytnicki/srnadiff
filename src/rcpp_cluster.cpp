#include <Rcpp.h>
#include <vector>
#include <tuple>
#include <functional>
#include <algorithm>
#include <limits>
using namespace Rcpp;

struct Region {
    unsigned long start, end, value;
    Region (): start(0), end(0), value(0) {}
    Region (unsigned long s, unsigned long e, unsigned long v):
        start(s), end(e), value(v) {}
    bool empty () {
        return (start == 0);
    }
    unsigned long getSize () {
        return end - start + 1;
    }
    unsigned long getDifference (Region &r) {
        if ((empty()) || (r.empty()))
            return std::numeric_limits<unsigned long>::max();
        return ((start >= r.start)? start - r.start: r.start - start) +
               ((end   >= r.end  )? end   - r.end  : r.end   - end);
    }
    bool operator< (const Region& r) const {
        if ((start < r.start) ||
            (start == r.start) && (end < r.end) ||
            (start == r.start) && (end < r.end) && (value < r.value))
            return true;
    }
};

//' Compute unique counts.
//'
//' @param lengths          the sizes of the RLEs (one list per chromosome)
//' @param values           the values of the RLEs (one list per chromosome)
//' @param chromosomeSizes  the sizes of the chromosomes
//' @param minDepth         minimum coverage
//' @param minSize          minimum region size
//' @param maxSize          maximum region size
//' @param minDifference    minimum difference between 2 regions
//' @return                 selected regions
// [[Rcpp::export]]
List rcpp_clustering(ListOf < ListOf < IntegerVector > > &lengths,
                   ListOf < ListOf < IntegerVector > > &values,
                   IntegerVector &chromosomeSizes, int minDepth,
                   int minSize, int maxSize, int minDifference) {
    static const short nConditions = 2;
    std::vector < unsigned int > starts, ends;
    std::vector < std::string > chromosomes;
    size_t nChromosomes = chromosomeSizes.size();
    std::vector < Region > selectedRegions[nConditions];
    for (size_t chrId = 0; chrId < nChromosomes; ++chrId) {
        for (int condition = 0; condition < nConditions; ++condition) {
            std::vector < Region > theseRegions;
            size_t        nValues =  lengths[condition][chrId].size();
            unsigned long start   = 1;
            for (size_t index = 0; index < nValues; ++index) {
                unsigned long length = lengths[condition][chrId][index];
                unsigned long value  = values[condition][chrId][index];
                unsigned long end    = start + length - 1;
                if (value >= minDepth) {
                    for (std::vector< Region >::reverse_iterator rit =
                         theseRegions.rbegin(); (rit != theseRegions.rend()) &&
                         (rit->end == start - 1) && (rit->value < value);
                         ++rit) {
                        rit->end = end;
                    }
                    theseRegions.emplace_back(start, end, value);
                }
                start = end + 1;
            }
            bool first = true;
            for (Region &region: theseRegions) {
                unsigned long size = region.getSize();
                if ((minSize <= size) && (size <= maxSize)) {
                    unsigned long difference = (first) ? minDifference + 1 :
                        region.getDifference(selectedRegions[condition].back());
                    if (difference >= minDifference) {
                        selectedRegions[condition].push_back(region);
                    }
                }
                first = false;
            }
        }
        std::vector < Region > mergedRegions;
        mergedRegions.reserve(selectedRegions[0].size() +
            selectedRegions[1].size());
        std::vector < Region >::iterator its[2] = {selectedRegions[0].begin(),
                                                   selectedRegions[1].begin()};
        while ((its[0] != selectedRegions[0].end()) &&
               (its[1] != selectedRegions[1].end())) {
            if (*its[0] < *its[1]) {
                mergedRegions.push_back(*its[0]);
                ++its[0];
            }
            else {
                mergedRegions.push_back(*its[1]);
                ++its[1];
            }
        }
        std::string chromosome = as < std::string >(as < CharacterVector> (
            chromosomeSizes.names())[chrId]);
        Region previousRegion;
        for (Region &region: mergedRegions) {
            if (region.getDifference(previousRegion)) {
                starts.push_back(region.start);
                ends.push_back(region.end);
                chromosomes.push_back(chromosome);
                previousRegion = region;
            }
        }
    }
    return DataFrame::create(_["seqnames"] = chromosomes,
                             _["start"]    = starts ,
                             _["end"]      = ends);
}
