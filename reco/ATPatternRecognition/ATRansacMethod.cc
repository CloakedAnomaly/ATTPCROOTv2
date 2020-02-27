#include "ATRansacMethod.hh"

// FairRoot classes
#include "FairRuntimeDb.h"
#include "FairRun.h"

ClassImp(ATPATTERN::ATRansacMethod)


ATPATTERN::ATRansacMethod::ATRansacMethod()
{
//initialize data members
}

ATPATTERN::ATRansacMethod::~ATRansacMethod()
{

}

void ATPATTERN::ATRansacMethod::SetIterations(int iterations) {fiterations = iterations;}
void ATPATTERN::ATRansacMethod::SetRatio(float ratio) {fratio = ratio; }
void ATPATTERN::ATRansacMethod::SetThreshold(float threshold) {fthreshold = threshold; }

void ATPATTERN::ATRansacMethod::SetTrackCand(std::vector<ATTrack> tracks) {fTrackCand = tracks;}

std::vector<ATTrack> ATPATTERN::ATRansacMethod::GetTrackCand() {return fTrackCand;}

std::vector<TVector3> ATPATTERN::ATRansacMethod::GetLine() {return fcoefficients;}

bool ATPATTERN::ATRansacMethod::FindTracks(ATEvent *event, ATPatternEvent *patternEvent)
{
  //initialize ransac parameters
 /* ATPATTERN::ransac_params parameters;
 
  parameters.iterations = 200;
  parameters.ratio = 0.2;
  parameters.threshold = 10;
  */
  
  
  std::vector<ATTrack> tracks;
  //create and fill vector contatining positions for each point 
  //std::vector<TVector3> points;
  
  std::vector<ATPATTERN::PointHit> points;
  eveToPointHit(event, points);

  for (int i=0; i<fiterations; i++){
	
  	//create new vector containing all of the points besides the two random points chosen for the line
        std::random_shuffle(points.begin(), points.end());
		
	TVector3 v0, v;
	find_3Dline(points[0].position, points[1].position, v0, v);
	std::vector<ATPATTERN::PointHit> inliers, outliers;
//        std::cout<< "Point Selected " << points[0].position.X() << " " << points[0].position.Y() << " " << points[0].position.Z() << std::endl;

	float num = 0.0;
			
	//loop over the testing points to find distance from the drawn line
	for (int j=2; j <= (points.size()); j++){
				
	    Double_t distance;	
	    find_distance(v0, v, points[j].position, distance);
            //std::cout << "distance = " << distance << std::endl;				
		
                //compare distance from line to threshold distance
		if (distance <= fthreshold){
		   inliers.push_back(points[j]);
		   num++;
//		std::cout<<"This threshold is working"<<std::endl;
                }//if

 		if (distance > fthreshold){
                   outliers.push_back(points[j]);
 		}//if	
			
	}//inner loop

      float ratio = num/(float)(points.size()-2);
      std::cout<<"ratio found is " << ratio << ", num =  " << num << ", size of vector " << points.size() << std::endl;		
   //compare ratio of inliers to desired_ratio
      if (ratio > fratio){
          ATTrack track_in, track_out;
          for (int in=0; in<inliers.size(); in++) {track_in.AddHit(inliers[in].hit);}//convert T3vector into athit
          for (int out=0; out<outliers.size(); out++) {track_out.AddHit(outliers[out].hit);}
          track_out.SetIsNoise(kTRUE);
          tracks.push_back(track_in);
          tracks.push_back(track_out);
          ATPATTERN::ATRansacMethod::saveLine(v0, v);
          patternEvent->SetTrackCand(tracks); 
          ATPATTERN::ATRansacMethod::SetTrackCand(tracks);
          std::cout << "Track is filled!" <<std::endl;   
          break;
      }//if 
   
   }//outer loop

}//FindTracks function


void ATPATTERN::ATRansacMethod::eveToPointHit(ATEvent event, std::vector<ATPATTERN::PointHit>& points){

  Int_t nHits = event.GetNumHits();

  for(Int_t iHit=0; iHit<nHits; iHit++){
        ATPATTERN::PointHit point;
        ATHit* hit = event.GetHit(iHit);
        point.hit = hit;
        Int_t PadNumHit = hit->GetHitPadNum();
        TVector3 position = hit->GetPosition();
        point.position = position;
        points.push_back(point);
//        std::cout<<"Point recieved from event"<<std::endl;
   }//loop

}//eveToPointHit


void ATPATTERN::ATRansacMethod::find_3Dline(TVector3 p1, TVector3 p2, TVector3& v0, TVector3& v)
//function for finding the slope and intercept of line that goes through two random points
{
  v0 = p1;
  v.SetX(p2.X() - p1.X());
  v.SetY(p2.Y() - p1.Y());
  v.SetZ(p2.Z() - p1.Z());

}


void ATPATTERN::ATRansacMethod::find_distance(TVector3 v0, TVector3 v, TVector3 p, Double_t& distance)
//function for finding the perpendicular line from the random line to each testing point
{
  Double_t t = -(v0.X()*p.X() + v0.Y()*p.Y() + v0.Z()*p.Z())/(v.X()*p.X() + v.Y()*p.Y() + v.Z()*p.Z());
  TVector3 p_norm;
  p_norm.SetX(v0.X() + t*v.X());
  p_norm.SetY(v0.Y() + t*v.Y()); 
  p_norm.SetZ(v0.Z() + t*v.Z());
//  std::cout<< "point " << p.X() << " " << p.Y() << " " << p.Z() << std::endl;
//  std::cout<< "normal point " << p_norm.X() << " " << p_norm.Y() << " " << p_norm.Z() << std::endl; 
  distance = TMath::Sqrt(TMath::Power((p.X()-p_norm.X()),2.0)+TMath::Power((p.Y()-p_norm.Y()),2.0)+TMath::Power((p.Z()-p_norm.Z()),2.0));
//  std::cout<< "distance : " << distance << std::endl;
}

void ATPATTERN::ATRansacMethod::saveLine(TVector3 v0, TVector3 v)
//function to save the slope and initial point of the line for ransac
{
  std::vector<TVector3> coefficients; 
  coefficients.push_back(v0);
  coefficients.push_back(v);
  fcoefficients = coefficients; 

}
