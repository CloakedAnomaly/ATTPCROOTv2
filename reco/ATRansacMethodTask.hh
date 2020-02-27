#ifndef ATRANSACMETHODTASK_H
#define ATRANSACMETHODTASK_H

// FAIRROOT classes
#include "FairTask.h"
#include "FairLogger.h"

// ATTPCROOT classes
#include "ATEvent.hh"
#include "ATProtoEvent.hh"
#include "ATDigiPar.hh"
#include "ATRansacMethod.hh"

// ROOT classes
#include "TClonesArray.h"

class ATRansacMethodTask : public FairTask {

  public:
    ATRansacMethodTask();
    ~ATRansacMethodTask();

    void SetPersistence(Bool_t value = kTRUE);
/*    void SetModelType(int model);
    void SetDistanceThreshold(Float_t threshold);
    void SetFullMode(); //Mode that calculates the RANSAC method for every potential line
    void SetMinHitsLine(Int_t nhits); //Set Mininum number of hits per line
*/

    virtual InitStatus Init();
    virtual void SetParContainers();
    virtual void Exec(Option_t *opt);
    virtual void Finish();

  private:
    FairLogger *fLogger;
    ATDigiPar *fPar;
    TClonesArray *fEventHArray;
    TClonesArray *fRansacArray;

    ATEvent *fEvent;

    Bool_t kIsPersistence;
/*    Bool_t kIsFullMode;
    int fRANSACModel;
    Float_t fRANSACThreshold;
    Int_t fMinHitsLine; // Minimum number of hits
*/



  ClassDef(ATRansacMethodTask, 1);
};

#endif
