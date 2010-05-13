#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/HitConverter.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"

using std::vector;

MyHit HitConverter::convert(const TransientTrackingRecHit* recHit, float estimate){
 	MyHit hit;

 	GlobalPoint myPoint = recHit->globalPosition();
 	hit.SetX(myPoint.x());
 	hit.SetY(myPoint.y());
 	hit.SetZ(myPoint.z());
 	GlobalError myError = recHit->globalPositionError();
 	hit.dxx = myError.cxx();
 	hit.dxy = myError.cyx();
 	hit.dxz = myError.czx();
 	hit.dyy = myError.cyy();
 	hit.dyz = myError.czy();
 	hit.dzz = myError.czz();

 	// Composite hit count
 	hit.theCompositeCount = recHit->transientHits().size();

 	hit.theEstimate = estimate;

 	//hit.printHit();
 	return hit;
}


/*
MyHit MyEventConverter::myHitConverter(const TrackingRecHit& recHit){
  MyHit hit;
  // INSERT BEGIN
  hit.theEstimate = -1;
  // save only valid hits
  if (!recHit.isValid()) return hit;

  TransientTrackingRecHit::RecHitPointer myTTRecHit = TTRHBuilder->build(&recHit);
  GlobalPoint myPoint = myTTRecHit->globalPosition();
  hit.x = myPoint.x();
  hit.y = myPoint.y();
  hit.z = myPoint.z();
  GlobalError myError = myTTRecHit->globalPositionError();
  hit.dxx = myError.cxx();
  hit.dxy = myError.cyx();
  hit.dxz = myError.czx();
  hit.dyy = myError.cyy();
  hit.dyz = myError.czy();
  hit.dzz = myError.czz();

  // Composite hit count
  hit.theCompositeCount = myTTRecHit->transientHits().size();

  hit.theEstimate = 0; // currently not implemented

  // INSERT END
  hit.printHit();
  return hit;
}
*/

void HitConverter::addHits(std::vector<MyHit>& hits, const Trajectory& trajectory, int trackLabel){
        // Loop over measurements to get hits and hit estimates
        vector<TrajectoryMeasurement> meas = trajectory.measurements();
        vector<TrajectoryMeasurement>::const_iterator iMeasEnd = meas.end();
        vector<TrajectoryMeasurement>::const_iterator iMeas;
        for(iMeas = meas.begin(); iMeas != iMeasEnd; ++iMeas) {
                const TrajectoryMeasurement myMeasurement = *iMeas;
                const TransientTrackingRecHit* myMeasuredRecHit = &(*myMeasurement.recHit());
                if (myMeasuredRecHit->isValid()) {
		  MyHit hit = convert(myMeasuredRecHit, iMeas->estimate());
		  hit.trackAssociationLabel = trackLabel;
                  hits.push_back(hit);
                }
        }
}
