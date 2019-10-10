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

std::vector<ATTrack> ATPATTERN::ATRansacMethod::GetTrackCand() {return fTrackCand;}

std::vector<TVector3> ATPATTERN::ATRansacMethod::GetLine() {return fcoefficients;}

bool ATPATTERN::ATRansacMethod::FindTracks(ATEvent &event, ATPatternEvent *patternEvent)
{
  //initialize ransac parameters
  ATPATTERN::ransac_params parameters;
 
  parameters.iterations = 0;
  parameters.ratio = 0;
  parameters.threshold = 0;
  
  std::vector<ATTrack> tracks;
  //create and fill vector contatining positions for each point 
  //std::vector<TVector3> points;
  
  std::vector<ATPATTERN::PointHit> points;
  eveToPointHit(event, points);

  for (int i=0; i<parameters.iterations; i++){
	
  	//create new vector containing all of the points besides the two random points chosen for the line
        std::random_shuffle(points.begin(), points.end());
		
	TVector3 v, v0;
	find_3Dline(points[0].position, points[1].position, v, v0);
	std::vector<ATPATTERN::PointHit> inliers, outliers;

	float num = 0.0;
			
	//loop over the testing points to find distance from the drawn line
	for (int j=2; j <= (points.size()); j++){
				
	    Double_t distance;	
	    find_distance(points[0].position, points[1].position, points[j].position, distance);
				
		//compare distance from line to threshold distance
		if (distance < parameters.threshold){
		   inliers.push_back(points[j]);
		   num++;
		}//if

 		if (distance > parameters.threshold){
                   outliers.push_back(points[j]);
 		}//if	
			
	}//inner loop

      float ratio = num/(float)(points.size()-2);
		
   //compare ratio of inliers to desired_ratio
      if (ratio > parameters.ratio){
          ATTrack track_in, track_out;
          for (int in=0; in<inliers.size(); in++) {track_in.AddHit(inliers[in].hit);}//convert T3vector into athit
          for (int out=0; out<outliers.size(); out++) {track_out.AddHit(outliers[out].hit);}
          track_out.SetIsNoise(kTRUE);
          tracks.push_back(track_in);
          tracks.push_back(track_out);
          ATPATTERN::ATRansacMethod::saveLine(v0, v);
          patternEvent->SetTrackCand(tracks);    
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
  Double_t t = -(v0.X()*p.X() + v0.Y()*p.Y() + v0.Z()*p.Z())/((v.X()-v0.X())*p.X() + (v.Y()-v0.Y())*p.Y() + (v.Z()-v0.Z())*p.Z());
  TVector3 p_norm;
  p_norm.SetX(v0.X() + t*v.X());
  p_norm.SetY(v0.Y() + t*v.Y()); 
  p_norm.SetZ(v0.Z() + t*v.Z());
  
  distance = sqrt(pow((p.X()-p_norm.X()),2)+pow((p.Y()-p_norm.Y()),2)+pow((p.Z()-p_norm.Z()),2));

}

void ATPATTERN::ATRansacMethod::saveLine(TVector3 v0, TVector3 v)
//function to save the slope and initial point of the line for ransac
{
  std::vector<TVector3> coefficients; 
  coefficients.push_back(v0);
  coefficients.push_back(v);
  fcoefficients = coefficients; 

}
