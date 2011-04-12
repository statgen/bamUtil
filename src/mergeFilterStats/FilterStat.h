#ifndef __FILTER_STAT_H
#define __FILTER_STAT_H

#define FILTER_STAT_COUNTS 6

#include <vector>
//#include <boost/thread/mutex.hpp>

class FilterStat {
 private:
  String sAnchorVcf;
  String sChrom;
  std::vector<int> vPos;
  std::vector<std::string> vAl1;
  std::vector<std::string> vAl2;
  std::vector<double> vAFs;
  std::vector<int> vCounts;  // RF,RR,AF,AR,OF,OR
  //boost::mutex mutex;

 public:
  bool loadAnchorVcf(const char* file);
  bool appendStatVcf(const char* file);
  bool writeMergedVcf(const char* outFile);
};

#endif


