#include "SAMain.hh"
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
#include "ATEvent.hh"
#include "ATTrack.hh"
#include "ATPad.hh"
#include "ATHit.hh"
#include "ATHoughSpace.hh"
#include "ATHoughSpaceLine.hh"
#include "ATHoughSpaceCircle.hh"
#include "ATPatternEvent.hh"
#include "ATRansacMethod.hh"
//ROOT
#include "TGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TMath.h"
#include "TF1.h"
#include "TAxis.h"

int main(int argc, char* argv[])
{
   TApplication app("app",&argc,argv);

   gSystem->Load("libATTPCReco.so");

   TStopwatch timer;
   timer.Start();
 
   FairRunAna* run = new FairRunAna(); //Forcing a dummy run
   TString FileName = "../output_proto.root";
   std::cout<<" Opening File : "<<FileName.Data()<<std::endl;
   TFile* file = new TFile(FileName.Data(),"READ");

   TTree* tree = (TTree*) file -> Get("cbmsim");
   Int_t nEvents = tree -> GetEntries();
   std::cout<<" Number of events : "<<nEvents<<std::endl;

   TTreeReader Reader1("cbmsim", file);
   TTreeReaderValue<TClonesArray> eventArray(Reader1, "ATEventH");

   TGraph2D * gr = new TGraph2D();
   int fiterations = 10000;
   float fthreshold = 20;
   float fratio = 0.35;
   TVector3 v0, v;

        for(Int_t i=0;i<nEvents;i++){
          //while (Reader1.Next()) {

              Reader1.Next();
              ATEvent* event = (ATEvent*) eventArray->At(0);

              if(i==20){//reading only event 20
              
	      std::vector<ATTrack> tracks;
             // std::vector<TVector3> coefficients; 


              if(event!=NULL){
                Int_t nHits = event->GetNumHits();
                std::cout<<" Event number "<<i<<" Number of hits "<<nHits<<"\n";
		
                std::vector<PointHit> points;
		eveToPointHit(event, points);

                for (int i=0; i<fiterations; i++){
              
                      //create new vector containing all of the points besides the two random points chosen for the line
                      //std::random_shuffle(points.begin(), points.end());
                      int x = std::rand()%points.size(), y = std::rand()%points.size();
                      
		      std::cout << x << " " << y << std::endl;
		      
			//TVector3 v0, v;
                      find_3Dline(points[x].position, points[y].position, v0, v);
                      std::vector<PointHit> inliers, outliers;
              //        std::cout<< "Point Selected " << points[0].position.X() << " " << points[0].position.Y() << " " << points[0].position.Z() << std::endl;
              
                      float num = 0.0;
              
                      //loop over the testing points to find distance from the drawn line
                      for (int j=0; j < (points.size()); j++){
                          if (j != x and j != y){
                          Double_t distance;
                          find_distance(v0, v, points[j].position, distance);
                      //    std::cout << "distance = " << distance << std::endl;
              
                              //compare distance from line to threshold distance
                              if (distance <= fthreshold){
                                 inliers.push_back(points[j]);
                                 num++;
              //              std::cout<<"This threshold is working"<<std::endl;
                              }//if
              
                              if (distance > fthreshold){
                                 outliers.push_back(points[j]);
                              }//if
                          }//if j
                      }//inner loop
              
                    float ratio = num/(float)(points.size()-2);
                   // std::cout<<"ratio found is " << ratio << ", num =  " << num << ", size of vector " << points.size() << std::endl;
                 //compare ratio of inliers to desired_ratio
                    if (ratio > fratio){
                        std::cout<<"ratio found is " << ratio << ", num =  " << num << ", size of vector " << points.size() << std::endl;
			inliers.push_back(points[x]);
			inliers.push_back(points[y]);
			ATTrack track_in;//, track_out;
                        for (int in=0; in<inliers.size(); in++) {track_in.AddHit(&inliers[in].hit);}//convert T3vector into athit
                        //for (int out=0; out<outliers.size(); out++) {track_out.AddHit(outliers[out].hit);}
                        //track_out.SetIsNoise(kTRUE);
                        tracks.push_back(track_in);
                        //tracks.push_back(track_out);
                        //ATPATTERN::ATRansacMethod::saveLine(v0, v);
                        //patternEvent->SetTrackCand(tracks);
                        //ATPATTERN::ATRansacMethod::SetTrackCand(tracks);
                        std::cout << "Track is filled!" <<std::endl;
                        break;
                    }//if
              
                 }//outer loop
              
                
              
             }//if event=!null
       	      gErrorIgnoreLevel=kFatal;
       
              std::vector<ATHit> *HitArray = tracks.at(0).GetHitArray();
              std::cout << "size of array : " << HitArray->size() << std::endl;     
            
              for(Int_t N=0;N<HitArray->size();N++){
              //std::cout << "I am finally working" << std::endl;
              ATHit hit = HitArray->at(N);
              TVector3 pos = hit.GetPosition();
              gr->SetPoint(N,pos.X(),pos.Y(),pos.Z());
             }

             gr->Draw("p0");

	     int n = 1000;
                double t0 = 0;
                double dt = 1000;
                TPolyLine3D *l = new TPolyLine3D(n);
                for (int j = -n; j <n;++j) {
                   double t = t0+ dt*j/n;
                   double x,y,z;
                   SetLine(t,v0,v,x,y,z);
                   l->SetPoint(j,x,y,z);
                   //std::cout<<" x : "<<x<<" y : "<<y<<"  z : "<<z<<std::endl;
                }
                l->SetLineColor(kRed);
                l->Draw("same");
	    	
}//if
}//for
   app.Run();


   return 0;

}
 

void eveToPointHit(ATEvent *event, std::vector<PointHit>& points){

  Int_t nHits = event->GetNumHits();

  for(Int_t iHit=0; iHit<nHits; iHit++){
        PointHit point;
        ATHit* hit = event->GetHit(iHit);
        point.hit = hit;
        Int_t PadNumHit = hit->GetHitPadNum();
        TVector3 position = hit->GetPosition();
        point.position = position;
        points.push_back(point);
//        std::cout<<"Point recieved from event"<<std::endl;
   }//loop

}//eveToPointHit


void find_3Dline(TVector3 p1, TVector3 p2, TVector3& v0, TVector3& v)
//function for finding the slope and intercept of line that goes through two random points
{
  v0 = p1;
  v.SetX(p2.X() - p1.X());
  v.SetY(p2.Y() - p1.Y());
  v.SetZ(p2.Z() - p1.Z());

}


void find_distance(TVector3 v0, TVector3 v, TVector3 p, Double_t& distance)
//function for finding the perpendicular line from the random line to each testing point
{
  TVector3 first = (p-v0);
 // std::cout << first.X() << " " << first.Y() << " " << first.Z() << std::endl;
  TVector3 second = (p-(v+v0));
 // std::cout << second.X() << " " << second.Y() << " " << second.Z() << std::endl;
  TVector3 top = first.Cross(second);
 // std::cout << top.X() << " " << top.Y() << " " << top.Z() << std::endl;
  distance = top.Mag()/v.Mag();
  


/*  Double_t t = -(v0.X()*p.X() + v0.Y()*p.Y() + v0.Z()*p.Z())/(v.X()*p.X() + v.Y()*p.Y() + v.Z()*p.Z());
  TVector3 p_norm;
  p_norm.SetX(v0.X() + t*v.X());
  p_norm.SetY(v0.Y() + t*v.Y());
  p_norm.SetZ(v0.Z() + t*v.Z());
//  std::cout<< "point " << p.X() << " " << p.Y() << " " << p.Z() << std::endl;
//  std::cout<< "normal point " << p_norm.X() << " " << p_norm.Y() << " " << p_norm.Z() << std::endl;
  distance = TMath::Sqrt(TMath::Power((p.X()-p_norm.X()),2.0)+TMath::Power((p.Y()-p_norm.Y()),2.0)+TMath::Power((p.Z()-p_norm.Z()),2.0));
//  std::cout<< "distance : " << distance << std::endl;
*/
}



void SetLine(double t, TVector3 v0, TVector3 v, double &x, double &y, double &z)
{
      // a parameteric line is define from 6 parameters but 4 are independent
      // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
      // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
      x = v0[0] + v[0]*t;
      y = v0[1] + v[1]*t;
      z = v0[2] + v[2]*t;
              
} 
