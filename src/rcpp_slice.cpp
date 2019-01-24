#include <Rcpp.h>
#include <vector>
#include "rcpp_utils.cpp"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

struct Region {
    long start, end, value;
    Region (): start(0), end(0), value(0) {}
    Region (long s, long e, long v):
        start(s), end(e), value(v) {}
    bool empty () {
        return (start == 0);
    }
    long getSize () {
        return end - start + 1;
    }
    long getDifference (Region &r) {
        if ((empty()) || (r.empty()))
            return std::numeric_limits<long>::max();
        return ((start >= r.start)? start - r.start: r.start - start) +
               ((end   >= r.end  )? end   - r.end  : r.end   - end);
    }
    bool operator< (const Region& r) const {
        return (((start < r.start) ||
            ((start == r.start) && (end < r.end)) ||
            ((start == r.start) && (end < r.end) && (value < r.value))));
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
List rcpp_slice(ListOf < ListOf < IntegerVector > > &lengths,
                   ListOf < ListOf < IntegerVector > > &values,
                   IntegerVector &chromosomeSizes, int minDepth,
                   int minSize, int maxSize, int minDifference) {
    static const short nConditions = 2;
    std::vector < int > starts, ends;
    std::vector < std::string > chromosomes;
    int nChromosomes = chromosomeSizes.size();
    std::vector < Region > selectedRegions[nConditions];
    std::string chromosome;
    for (int chrId = 0; chrId < nChromosomes; ++chrId) {
        chromosome = as < std::string >(as < CharacterVector> (
            chromosomeSizes.names())[chrId]);
        for (int condition = 0; condition < nConditions; ++condition) {
            std::vector < Region > theseRegions;
            int  nValues = lengths[condition][chrId].size();
            long start   = 1;
            int  rLast   = 0;
            selectedRegions[condition].clear();
            for (int index = 0; index < nValues; ++index) {
                long length    = lengths[condition][chrId][index];
                long value     = values[condition][chrId][index];
                long end       = start + length - 1;
                int  nextRLast = theseRegions.size();
                for (int rit = theseRegions.size()-1; rit >= rLast; --rit) {
                    Region &region = theseRegions[rit];
                    if ((region.end == start - 1) && (region.value <= value)) {
                        region.end = end;
                        nextRLast = rit;
                    }
                }
                rLast = nextRLast;
                if (value >= minDepth) {
                    theseRegions.emplace_back(start, end, value);
                }
                start = end + 1;
            }
            for (Region &region: theseRegions) {
                long size = region.getSize();
                if ((minSize <= size) && (size <= maxSize)) {
                    long difference = minDifference + 1;
                    for (auto rit = selectedRegions[condition].rbegin();
                         rit != selectedRegions[condition].rend(); ++rit) {
                        long d = region.getDifference(*rit);
                        difference = std::min<long>(difference, d);
                        if (d > 10 * minDifference) {
                            break;
                        }
                    }
                    if (difference >= minDifference) {
                        selectedRegions[condition].push_back(region);
                    }
                }
            }
        }
        std::vector < Region > mergedRegions;
        mergedRegions.reserve(selectedRegions[0].size() +
            selectedRegions[1].size());
        std::vector < Region >::iterator its[2] =
            {selectedRegions[0].begin(), selectedRegions[1].begin()};
        while ((its[0] != selectedRegions[0].end()) &&
               (its[1] != selectedRegions[1].end())) {
            int condition = (*its[0] < *its[1])? 0: 1;
            mergedRegions.push_back(*its[condition]);
            ++its[condition];
        }
        for (int condition = 0; condition < nConditions; ++condition) {
            while (its[condition] != selectedRegions[condition].end()) {
                mergedRegions.push_back(*its[condition]);
                ++its[condition];
            }
        }
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
