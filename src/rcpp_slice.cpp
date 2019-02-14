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

struct ValuedRegion {
    int start, end;
    double meanLFC;
};

void getPositions (int start, int end,
                   NumericVector &values, IntegerVector &lengths,
                   int &startId, int &endId) {
    int pos = 0;
    int id;
    for (id = 0; id < lengths.size(); ++id) {
        pos += lengths[id];
        if (pos > start) {
            pos -= lengths[id];
            --id;
            break;
        }
    }
    startId  = id;
    for (++id; id < lengths.size(); ++id) {
        pos += lengths[id];
        if (pos > end) {
            pos -= lengths[id];
            --id;
            break;
        }
    }
    endId  = id;
}

double getMeanLFC (int start, int end,
                   NumericVector &values, IntegerVector &lengths, int startId, int startPos, int endId, int endPos) {
    double nValues = 0.0;
    double sum     = 0.0;
    int    id, nValuesHere;
    nValuesHere = lengths[startId] - (start - startPos);
    nValues     = nValuesHere;
    sum         = nValuesHere * abs(values[startId]);
    for (id = startId+1; id < endId; ++id) {
        nValues += lengths[id];
        sum     += lengths[id] * abs(values[id]);
    }
    nValuesHere = lengths[endId] - (endPos - end);
    nValues    += nValuesHere;
    sum        += nValuesHere * abs(values[endId]);
    return sum/nValues;
}

// [[Rcpp::export]]
List rcpp_slice2(S4 &logFoldChanges, S4 &regions, int minLength, int maxLength, double minLFC) {
    IntegerVector outputStarts, outputEnds;
    StringVector outputChromosomes;
    S4 regionsSeqnames = regions.slot("seqnames");
    IntegerVector regionsSeqnamesValues = regionsSeqnames.slot("values");
    StringVector regionsSeqnamesLevels = regionsSeqnamesValues.attr("levels");
    IntegerVector regionsSeqnamesLengths = regionsSeqnames.slot("lengths");
    S4 regionsRanges = regions.slot("ranges");
    IntegerVector regionStarts = regionsRanges.slot("start");
    IntegerVector regionWidths = regionsRanges.slot("width");
    List rleList = logFoldChanges.slot("listData");
    int seqnamesId = 0;
    int regionStart, regionEnd;
    //Rcout << "regionsSeqnamesValues: " << regionsSeqnamesValues << std::endl;
    //Rcout << "regionsSeqnamesLevels: " << regionsSeqnamesLevels << std::endl;
    //Rcout << "regionsSeqnamesLengths: " << regionsSeqnamesLengths << std::endl;
    String chromosome = regionsSeqnamesLevels[regionsSeqnamesValues[0]-1];
    int seqnamesRemaining = regionsSeqnamesLengths[0];
    for (int regionId = 0; regionId < regionStarts.size(); ++regionId) {
        //Rcout << "region id: " << regionId << std::endl;
        if (seqnamesRemaining == 0) {
            ++seqnamesId;
            //Rcout << "seqnames id to: " << seqnamesId << std::endl;
            seqnamesRemaining = regionsSeqnamesLengths[seqnamesId];
            chromosome = regionsSeqnamesLevels[regionsSeqnamesValues[seqnamesId]];
        }
        regionStart = regionStarts[regionId];
        regionEnd   = regionStarts[regionId] + regionWidths[regionId] - 1;
        --seqnamesRemaining;
        S4     rle        = rleList[chromosome];
        NumericVector values  = rle.slot("values");
        IntegerVector lengths = rle.slot("lengths");
        int startId, endId, currentId, currentStart, currentPos;
        double currentLength, currentSum, nextLength, nextSum;
        getPositions(regionStart, regionEnd, values, lengths, startId, endId);
        currentPos = regionStart;
        currentId  = startId;
        while (currentId <= endId) {
            //Rcout << "  Starting with " << values[currentId] << " @ " << currentId << std::endl;
            while ((currentId <= endId) && (abs(values[currentId]) < minLFC)) {
                //Rcout << "    Moving to next: " << values[currentId] << " @ " << currentId << std::endl;
                currentPos += lengths[currentId];
                ++currentId;
            }
            if (currentId > endId) {
                //Rcout << "    Exiting case 1" << std::endl;
                break;
            }
            currentStart = currentPos;
            currentLength = currentSum = 0.0;
            //Rcout << "    Over step 1 with " << currentStart << std::endl;
            while (currentId <= endId) {
                //Rcout << "    Next step: " << values[currentId] << std::endl;
                nextSum    = currentSum + abs(values[currentId]) * lengths[currentId];
                nextLength = currentLength + lengths[currentId];
                if (nextSum / nextLength < minLFC) {
                    //Rcout << "    Exiting case 2" << std::endl;
                    break;
                }
                currentSum    = nextSum;
                currentLength = nextLength;
                currentPos   += lengths[currentId];
                ++currentId;
            }
            if ((minLength <= currentPos - currentStart - 1) && (currentPos - currentStart - 1 <= maxLength)) {
                //Rcout << "    Storing " << currentStart << "-" << currentPos << std::endl;
                outputChromosomes.push_back(chromosome);
                outputStarts.push_back(currentStart);
                outputEnds.push_back(currentPos);
            }
        }
        /*
        double meanLFC = getMeanLFC(regionStart, regionEnd, values, lengths, startId, startPos, endId, endPos);
        ValuedRegion firstRegion = { startId, endId, meanLFC };
        std::vector < ValuedRegion > selectedRegions = { firstRegion };
        double previousMinLFC = -1.0;
        bool toBeSplit = false;
        Rcout << "Now region " << regionStart << " (" << startId << ") - " <<  regionEnd << " (" << endId << "), mean LFC: " << meanLFC << std::endl;
        Rcout << "Profile:" << std::endl;
        for (int i = startId; i <= endId; ++i) Rcout << "  " << values[i] << " " << lengths[i] << std::endl;
        while (! selectedRegions.empty()) {
            std::vector < ValuedRegion > splittedRegions;
            bool inRegion = false;
            int  endRegionId = 0;
            double currentSum = 0.0;
            double currentSize = 0.0;
            double currentLFC;
            double currentMinLFC = std::numeric_limits<double>::max();
            for (int i = selectedRegions.back().start; i <= selectedRegions.back().end; ++i) {
                if ((abs(values[i]) > previousMinLFC) && (abs(values[i]) < currentMinLFC)) {
                    currentMinLFC = abs(values[i]);
                }
            }
            Rcout << "  min LFC: current = " << currentMinLFC << ", previous = " << previousMinLFC << std::endl;
            if (currentMinLFC == std::numeric_limits<double>::max()) {
                for (ssize_t i = selectedRegions.size()-1; i >= 0; --i) {
                    Rcout << "  Writing element " << i << ": " << selectedRegions[i].start << "-" << selectedRegions[i].end << " with ref " << startPos << " <-> " << startId << std::endl;
                    int start = 0, end = 0, pos = startPos;
                    for (int i = startId; i < selectedRegions[i].end; ++i) {
                        pos += lengths[i];
                        if (i == selectedRegions[i].start) {
                            start = pos;
                        }
                        end = pos;
                    }
                    outputStarts.push_back(start);
                    outputEnds.push_back(end);
                    outputChromosomes.push_back(chromosome);
                }
                Rcout << "  Written " << selectedRegions.size() << " elements" << std::endl;
                break;
            }
            double minLFC = currentMinLFC;
            for (int i = selectedRegions.back().end; i >= selectedRegions.back().start; --i) {
                Rcout << "    " << i << " -> " << values[i] << std::endl;
                if ((abs(values[i]) > currentMinLFC) && (! inRegion)) {
                    Rcout << "      entering region" << std::endl;
                    inRegion = true;
                    endRegionId = i;
                    currentSum = currentSize = 0.0;
                }
                else if ((abs(values[i]) <= currentMinLFC) && (inRegion)) {
                    inRegion = false;
                    currentLFC = currentSum / currentSize;
                    Rcout << "      exiting region: " << currentLFC << std::endl;
                    splittedRegions.push_back( { i+1, endRegionId, currentLFC } );
                    if (currentLFC < currentMinLFC) {
                        toBeSplit = true;
                    }
                }
                if (inRegion) {
                    Rcout << "      in region" << std::endl;
                    currentSize += lengths[i];
                    currentSum  += lengths[i] * abs(values[i]);
                }
            }
            if (inRegion) {
                currentLFC = currentSum / currentSize;
                Rcout << "      exiting final region: " << currentLFC << std::endl;
                splittedRegions.push_back( { selectedRegions.back().start, endRegionId, currentLFC } );
                if (currentLFC < currentMinLFC) {
                    toBeSplit = true;
                }
            }
            Rcout << "  Found " << selectedRegions.size() << " sub regions" << std::endl;
            for (auto &r: selectedRegions) Rcout << "    " << r.start << "-" << r.end << " (" << r.meanLFC << ")" << std::endl;
            if (toBeSplit) {
                Rcout << "  To be split" << std::endl;
                selectedRegions.pop_back();
                selectedRegions.reserve(selectedRegions.size() + splittedRegions.size());
                selectedRegions.insert(selectedRegions.end(), splittedRegions.begin(), splittedRegions.end());
            }
            previousMinLFC = minLFC;
        }
     */
    }
    //Rcout << "Output:\n\t" << outputChromosomes << "\n\t" << outputStarts << "\n\t" << outputEnds << std::endl;
    return DataFrame::create(_["seqnames"] = outputChromosomes,
                             _["start"]    = outputStarts,
                             _["end"]      = outputEnds);
}

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
