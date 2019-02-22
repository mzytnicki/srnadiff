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
                   NumericVector &values, IntegerVector &lengths,
                   int startId, int startPos, int endId, int endPos) {
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
List rcpp_ir(S4 &logFoldChanges, int minLength, int maxLength, double minLFC) {
    List          rleList      = logFoldChanges.slot("listData");
    StringVector  chromosomes  = rleList.names();
    int           nChromosomes = chromosomes.size();
    IntegerVector outputStarts, outputEnds;
    StringVector outputChromosomes;
    for (int chromosomeId = 0; chromosomeId < nChromosomes; ++chromosomeId) {
        S4            rle           = rleList[chromosomeId];
        NumericVector values        = rle.slot("values");
        IntegerVector lengths       = rle.slot("lengths");
        int           regionStart   = -1;
        int           pos           = 1;
        double        currentSum    = 0.0, nextSum;
        int           currentLength = 0,   nextLength;
        double        currentValue  = 0;
        int           currentSign   = 0;
        for (int rleId = 0; rleId < values.size(); ++rleId) {
            if ((regionStart == -1) || (currentSign * values[rleId] >= 0)) {
                currentValue = abs(values[rleId]);
            }
            else {
                currentValue = 0;
            }
            nextLength = currentLength + lengths[rleId];
            nextSum    = currentSum    + currentValue * lengths[rleId];
            if ((currentValue >= minLFC) && (regionStart == -1)) {
                regionStart   = pos;
                nextSum       = currentValue;
                nextLength    = lengths[rleId];
                currentSign   = (values[rleId] >= 0)? 1: -1;
            }
            else if ((currentValue < minLFC) && (regionStart != -1) &&
                        (nextSum / nextLength < minLFC)) {
                if ((minLength <= pos - regionStart - 1) &&
                        (pos - regionStart - 1 <= maxLength)) {
                    outputChromosomes.push_back(chromosomes[chromosomeId]);
                    outputStarts.push_back(regionStart);
                    outputEnds.push_back(pos-1);
                }
                regionStart = -1;
            }
            currentLength = nextLength;
            currentSum    = nextSum;
            pos           += lengths[rleId];
        }
        if ((regionStart != -1) && (minLength <= pos - regionStart - 1) &&
                (pos - regionStart - 1 <= maxLength)) {
            outputChromosomes.push_back(chromosomes[chromosomeId]);
            outputStarts.push_back(regionStart);
            outputEnds.push_back(pos-1);
        }
    }
    return DataFrame::create(_["seqnames"] = outputChromosomes,
                             _["start"]    = outputStarts,
                             _["end"]      = outputEnds);
}
