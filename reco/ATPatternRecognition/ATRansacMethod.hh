/*******************************************************************
* Base class for ATRansacMethod                                    *
* Log: Class started 09-03-2019                                    *
* Author: G. Wilks                                                 *
********************************************************************/

#ifndef ATRANSACMETHOD_H
#define ATRANSACMETHOD_H

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

//System
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <future>
#include <cmath>

//ATTPCROOT
#include "ATPRA.hh"
#include "ATHit.hh"
#include "ATEvent.hh"
#include "ATPatternEvent.hh"
#include "ATDigiPar.hh"
#include "AtTpcMap.h"
#include "ATTrack.hh"
#include "TObject.h"

// FairRoot classes
#include "FairRootManager.h"
#include "FairLogger.h"
/*
//PCL
#include <pcl/common/common.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/io/pcd_io.h>
*/
/*
//trackfinder
#include "hc.hh"
#include "msd.hh"
#include "smoothenCloud.hh"
*/
#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

/*
struct Point {
  float x;
  float y;
  float z;
  float t;
  float A;
};
*/
namespace ATPATTERN{

struct ransac_params {
   float ratio;
   float iterations;
   float threshold;
};

struct PointHit {
   TVector3 position;
   ATHit* hit;
};


class ATRansacMethod : public ATPRA
{
     
  public:
      ATRansacMethod();
      ~ATRansacMethod();

      bool FindTracks(ATEvent &event, ATPatternEvent *patternEvent);
      std::vector<ATTrack> GetTrackCand();
      std::vector<TVector3> GetLine();
      
  private:
      
      void find_3Dline(TVector3 p1, TVector3 p2, TVector3& v0, TVector3& v);
      
      void find_distance(TVector3 v0, TVector3 v, TVector3 p, Double_t& distance);
      
      void eveToPointHit(ATEvent event, std::vector<ATPATTERN::PointHit>& points);
     
      void saveLine(TVector3 v0, TVector3 v);

      std::vector<ATTrack> fTrackCand; //Candidate tracks

      std::vector<TVector3> fcoefficients; //line coefficients
      
      ClassDef(ATRansacMethod, 1);

};

}//namespace

#endif
