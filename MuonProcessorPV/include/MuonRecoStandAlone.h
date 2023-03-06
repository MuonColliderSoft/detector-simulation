#ifndef MuonRecoStandAlone_h
#define MuonRecoStandAlone_h 1

#include "marlin/Processor.h"

#include "lcio.h"

//#include <EVENT/LCCollection.h>
//#include <EVENT/MCParticle.h>
////#include <EVENT/ReconstructedParticle.h>


//#include "TH2F.h"
//#include "TH1D.h"
//#include "TProfile2D.h"
//#include "TTree.h"
//#include "TFile.h"
//#include "TLorentzVector.h"
//#include <map>

using namespace lcio ;
using namespace marlin ;

//class TTree;


class MuonRecoStandAlone : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new MuonRecoStandAlone ; }
  
  MuonRecoStandAlone() ;
  
  MuonRecoStandAlone(const MuonRecoStandAlone&) = delete;
  MuonRecoStandAlone& operator=(const MuonRecoStandAlone&) = delete;
  
  // Initialisation - run at the beginning to start histograms, etc.
  virtual void init() ;
  
  // Called at the beginning of every run
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  // Run over each event - the main algorithm
  virtual void processEvent( LCEvent * evt ) ;
  
  // Run at the end of each event
  virtual void check( LCEvent * evt ) ;
  
  // Called at the very end for cleanup, histogram saving, etc.
  virtual void end() ;
  
protected:

  std::string m_muonColName="";
  std::string m_outputClusterCollection="";
  
  
  //std::string m_innerecColName="";
 // std::string m_innerbaColName="";
 // std::string m_outerbaColName="";
 // std::string m_outerecColName="";
 // std::string m_vertexbaColName="";
 // std::string m_vertexecColName="";
  std::vector<std::string> m_inputTrackerHitPlaneCollNames{} ;
  std::vector<std::string> m_outputTrackerHitPlaneCollNames{} ;
  double m_cutdRMuon{};
  double m_cutdRTracker{};
  int m_deltaModule{};
  int m_cutLayerEndcap{};
  int m_cutLayerBarrel{};
  // Run and event counters
  int m_eventNumber=0;
  int m_runNumber=0;
  int m_trk_flag{};
  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);
  
  
  
} ;

#endif



