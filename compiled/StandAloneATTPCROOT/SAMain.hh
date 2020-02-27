#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <vector>

#include "TClonesArray.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreePlayer.h"
#include "TTreeReaderValue.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStopwatch.h"

#include "FairRootManager.h"
#include "FairLogger.h"
#include "FairRun.h"
#include "FairRunAna.h"

#include "AtTpcPoint.h"
#include <fstream>
//#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
//#include <vector>
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
//#include "FairRootManager.h"
//#include "FairLogger.h"
#include "TMath.h"


/*
#include <pcl/console/parse.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/sample_consensus/sac_model_sphere.h>
#include <pcl/sample_consensus/sac_model_circle.h>
#include <boost/thread/thread.hpp>
#include <pcl/ModelCoefficients.h>
#include <pcl/segmentation/sac_segmentation.h>

*/


struct PointHit {
   TVector3 position;
   ATHit hit;
};


void find_3Dline(TVector3 p1, TVector3 p2, TVector3& v0, TVector3& v);

void find_distance(TVector3 v0, TVector3 v, TVector3 p, Double_t& distance);

void eveToPointHit(ATEvent *event, std::vector<PointHit>& points);

//void saveLine(TVector3 v0, TVector3 v);

void SetLine(double t, TVector3 v0, TVector3 v, double &x, double &y, double &z);
