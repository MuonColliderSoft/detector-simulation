#ifndef Analyzer_h
#define Analyzer_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>


#include "TH2F.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <map>

using namespace lcio ;
using namespace marlin ;

class TTree;


class Analyzer : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new Analyzer ; }
  
  Analyzer() ;

  Analyzer(const Analyzer&) = delete;
  Analyzer& operator=(const Analyzer&) = delete;
  
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
  double _magneticField = 0.0;
  std::string m_inputMCParticleCollection="";

  std::string m_inputRECOParticleCollection="";
  std::string m_inputRECOClusterCollection="";
  std::string m_muonColName="";

  std::string m_inputCTCollection="";
  std::string m_inputClusterCollection="";

  std::string m_rootFileName="";
  
  // Run and event counters
  int m_eventNumber=0;
  int m_runNumber=0;
  //=============================================
  //MCParticle
  //Primary particles  
  std::vector<int>*   m_mcpPDGID=NULL;
  std::vector<float>* m_mcpE=NULL;
  std::vector<float>* m_mcpPx=NULL;
  std::vector<float>* m_mcpPy=NULL;
  std::vector<float>* m_mcpPz=NULL;
  std::vector<float>* m_mcpTime=NULL;
  
 
  //PANDORA PARTICLES
  //Cluster 
  std::vector<float>*  m_E_pfoclus;
  std::vector<float>*  m_x_pfoclus;
  std::vector<float>*  m_y_pfoclus;
  std::vector<float>*  m_z_pfoclus;
  std::vector<float>*  m_Eecal_pfoclus;
  std::vector<float>*  m_Ehcal_pfoclus;
  std::vector<float>*  m_Eyoke_pfoclus;
  std::vector<int>*    m_nhitstot_pfoclus; 
  std::vector<int>*    m_nhitsmuon_pfoclus;
  std::vector<float>*  m_phi_pfoclus;
  std::vector<float>*  m_theta_pfoclus;
  //Pandora Particles
  std::vector<float>*  m_type_pfopart;
  std::vector<float>*  m_E_pfopart;
  std::vector<float>*  m_px_pfopart;
  std::vector<float>*  m_py_pfopart;
  std::vector<float>*  m_pz_pfopart;
  std::vector<float>*  m_dR_pfopart;
  std::vector<int>*    m_flag_pfopart;

  std::vector<float>*  m_E_SAclus;
  std::vector<float>*  m_x_SAclus;
  std::vector<float>*  m_y_SAclus;
  std::vector<float>*  m_z_SAclus;
  std::vector<int>*    m_nhitstot_SAclus; 
  std::vector<float>*  m_phi_SAclus;
  std::vector<float>*  m_theta_SAclus;

  std::vector<float>* m_trkPt;
  std::vector<float>* m_trkTheta;
  std::vector<float>* m_trkPhi;
  std::vector<float>* m_trkCharge;
  
  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);
  
  TFile * m_rootFile=NULL;
  TTree* m_outputTree=NULL;

} ;

#endif



