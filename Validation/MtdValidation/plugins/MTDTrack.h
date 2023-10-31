#ifndef Validation_MtdValidation_MTDTrack_h
#define Validation_MtdValidation_MTDTrack_h

using namespace std;

struct MTDTrack {
  float pv_x;
  float pv_y;
  float pv_z;
  float pv_t;
  int pv_isFake;

  float genpv_x;
  float genpv_y;
  float genpv_z;
  float genpv_t;

  
  vector<float>  p;
  vector<float>  pt;
  vector<float>  eta;
  vector<float>  phi;
  vector<float>  x;
  vector<float>  y;
  vector<float>  z;
  vector<float>  dz;
  vector<float>  dxy;
  vector<float>  dzErr;
  vector<float>  dxyErr;
  vector<float>  chi2;
  vector<int>    ndof;
  vector<int>    nTrackerHits;
  vector<int>    nPixelBarrelHits;
  vector<int>    nPixelEndcapHits;
  vector<int>    nStripHits;
  vector<int>    nTimingBTLHits;
  vector<int>    nTimingETLHits;
  
  vector<float>  btlMatchChi2;
  vector<float>  btlMatchTimeChi2;
  vector<float>  etlMatchChi2;
  vector<float>  etlMatchTimeChi2;

  vector<float>  tmtd;
  vector<float>  t0;
  vector<float>  sigmat0;
  vector<float>  t0Pid;
  vector<float>  sigmat0Pid;
  vector<float>  pathlength;

  vector<float>  mtdTrackQualityMVA;

  vector<int> isMatchedToTP;
  vector<int> isMatchedToSimCluster;
  vector<int> isMatchedToGen;
    
};

#endif  //Validation_MtdValidation_MTDTrack_h
