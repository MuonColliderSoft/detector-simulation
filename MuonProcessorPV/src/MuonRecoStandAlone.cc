#include "CalorimeterHitType.h"
#include "MuonRecoStandAlone.h"
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/LCCollection.h>
#include <UTIL/BitField64.h>
#include <algorithm>
#include "TVector3.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <set>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <cmath> 
#include <iostream>
#include <cstdlib>
#include <cstdint>

#include <UTIL/BitSet32.h>
using namespace lcio ;

MuonRecoStandAlone aMuonRecoStandAlone;
MuonRecoStandAlone::MuonRecoStandAlone() : Processor("MuonRecoStandAlone"){  
  _description = "MuonRecoStandAlone processor selects muon hits in a cone opened around the hits in the first muon layer and makes a cluster. Then it selects SimTrackerHits in the cone identified by this cluster.";
  registerInputCollection( LCIO::CALORIMETERHIT,
			   "MUON" ,  
			   "Name of the MUON Digi collection",
			   m_muonColName,
			   std::string("MUON"));  
  registerProcessorParameter( "DeltaRCutMuon" ,
			      "Maximum angular distance between the hits and the particle direction in the muon system" ,
			      m_cutdRMuon,
			      double(0.1));
  registerProcessorParameter( "DeltaRCutTracker" ,
			      "Maximum angular distance between the hits and the particle direction in the tracker system" ,
			      m_cutdRTracker,
			      double(0.1));
  
  registerProcessorParameter( "CutLayerEndcap" ,
			      "Minimum number of muon system layers in endcap with hits" ,
			      m_cutLayerEndcap,
			      int(5));
    registerProcessorParameter( "CutLayerBarrel" ,
			      "Minimum number of muon system layers in barrel with hits" ,
			      m_cutLayerBarrel,
			      int(4));
  registerProcessorParameter( "DeltaModule" ,
			      "Define module range for cluster" ,
			      m_deltaModule,
			      int(1));
  registerOutputCollection( LCIO::CLUSTER,
			    "StandAloneClusterCollectionName",
			    "Name of the StandAloneCluster (muon hits) output collection",
			    m_outputClusterCollection,
			    std::string("StandAloneCluster"));
  registerInputCollections(LCIO::TRACKERHITPLANE,
			   "TrackerHitPlaneInputCollections",
			   "Name of the tracker hit plane input collections",
			   m_inputTrackerHitPlaneCollNames,
			   {});
			   
  //registerInputCollections(LCIO::SIMTRACKERHIT,
  //			    "TrackerSimHitInputCollections",
  //			     "Name of the tracker simhit input collections",
  //			     m_inputTrackerSimHitsCollNames,
  //			     {} );
  registerProcessorParameter("TrackerHitPlaneOutputCollections",
			     "Name of the tracker hit plane output collections",
			     m_outputTrackerHitPlaneCollNames,
			     {} );
			     
  registerProcessorParameter( "Tracker_flag" ,
			      "Decide if you want to select hit in the tracker collections" ,
			      m_trk_flag,
			      int(0));
			    
  
}

void MuonRecoStandAlone::init() {  
  printParameters() ;
  m_runNumber = 0 ;
  m_eventNumber = 0 ;
}


void MuonRecoStandAlone::processRunHeader( LCRunHeader*) {
  ++m_runNumber ;
}

void MuonRecoStandAlone::processEvent( LCEvent* evt ) {
  std::cout << "------------Evento " <<evt->getEventNumber() <<std::endl;
  m_eventNumber=evt->getEventNumber();
  m_runNumber=evt->getRunNumber();
  Int_t m_cutLayer=0;
  //Create output cluster collection
  LCCollectionVec *cluster = new LCCollectionVec(LCIO::CLUSTER);
  cluster->setFlag( cluster->getFlag() | 1 << LCIO::CLBIT_HITS );
  //Tracker hits input and output collections
  
  const unsigned int nTrackerHitCol = m_inputTrackerHitPlaneCollNames.size();
  std::vector<LCCollection*> inputHitColls(nTrackerHitCol);
  std::vector<LCCollectionVec*> outputTrackerHitColls(nTrackerHitCol);
  if(m_trk_flag>0){
  for (unsigned int icol=0; icol<nTrackerHitCol ; ++icol) {
    try {
      inputHitColls[icol] = evt->getCollection(m_inputTrackerHitPlaneCollNames[icol]); 
    }
    catch( lcio::DataNotAvailableException& e ) {
      streamlog_out(WARNING) << m_inputTrackerHitPlaneCollNames[icol]
			     << " collection not available" << std::endl;
      continue;
    }
    
  std::string encoderString = inputHitColls[icol]->getParameters().getStringVal( "CellIDEncoding" );
  outputTrackerHitColls[icol] = new LCCollectionVec( inputHitColls[icol]->getTypeName() );
  outputTrackerHitColls[icol]->parameters().setValue( "CellIDEncoding", encoderString );
  LCFlagImpl lcFlag_sim(inputHitColls[icol]->getFlag());
  outputTrackerHitColls[icol]->setFlag(lcFlag_sim.getFlag());
  } //end nTrackerHitCollections
 }
  //Open the Muon Calorimeter Hit collection
  LCCollection * muon_hits =0;
  try{
    muon_hits = evt->getCollection(m_muonColName) ;
  }
  catch( lcio::DataNotAvailableException &e ){
    streamlog_out(WARNING) << m_muonColName << " collection not available" << std::endl;
    muon_hits = NULL;
  }
  
  if(muon_hits==NULL) return; 

  int mod_module =13- m_deltaModule;
  
  //------------------------------------------------------------
  // Create collection of hits per each layer of the muon system
  //------------------------------------------------------------
  UTIL::BitField64 b(muon_hits->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ));
  int n_layers = 7;
  std::vector<std::vector<std::pair<CalorimeterHit*,float>>> muonhit_layer;
  muonhit_layer.resize(7);

  for(int mh=0; mh<muon_hits->getNumberOfElements(); mh++){
    CalorimeterHit* mu_hit = dynamic_cast<CalorimeterHit*>(muon_hits->getElementAt(mh));
    lcio::long64 val = lcio::long64( mu_hit->getCellID0() & 0xffffffff ) |  ( lcio::long64( mu_hit->getCellID1() ) << 32) ;
    b.setValue(val);
    unsigned int layer = b["layer"];
    //unsigned int system = b["system"];
    muonhit_layer[layer].push_back(std::make_pair(mu_hit, mu_hit->getEnergy()));
  } //end loop on CalorimeterHit
  
  //-------------------------------------
  // We order the hits in the first layer
  //-------------------------------------
  std::sort(muonhit_layer[0].begin(),muonhit_layer[0].end(),[](const std::pair<CalorimeterHit*, float>& lhs, const std::pair<CalorimeterHit*, float>& rhs) {return lhs.second > rhs.second;});
  
  //-----------------------------------------------------------------------------------
  // Loop over the hits on the first layer as the starting seed for the standalone track
  //-----------------------------------------------------------------------------------
  std::vector<std::vector<int>> tmp;
  tmp.resize(n_layers);
  int n_SAmuons = -1;
  //std::vector<std::set<int> > hits_to_save(nTrackerHitCol);
  
  for (unsigned int h=0; h<muonhit_layer[0].size(); h++){
    ClusterImpl *  cl= new ClusterImpl();
    float cl_energy=0.;
    float cl_centroid[3]={0,0,0};
    tmp.clear();
    if (muonhit_layer[0].at(h).second==0) continue;
    TVector3 starting_hit(muonhit_layer[0].at(h).first->getPosition()[0], muonhit_layer[0].at(h).first->getPosition()[1], muonhit_layer[0].at(h).first->getPosition()[2]);
    int layerhit_count=0;
    lcio::long64 val1 = lcio::long64(muonhit_layer[0].at(h).first->getCellID0() & 0xffffffff ) |  ( lcio::long64( muonhit_layer[0].at(h).first->getCellID1() ) << 32) ;
    b.setValue(val1);
    int side_startinghit = b["side"];
    if(side_startinghit==0) m_cutLayer=m_cutLayerBarrel;
    else if (side_startinghit==-2 || side_startinghit==1 ) m_cutLayer=m_cutLayerEndcap;
      
    unsigned int module_startinghit = b["module"];
    int n_tmp[7]={0,0,0,0,0,0,0};
   
    for(int l=0; l<n_layers; l++){ //loop on all the layers
      std::vector<int> tmp2; 
      int flag_hitlayer=0;
      for(unsigned int hl=0; hl<muonhit_layer[l].size(); hl++){
	TVector3 all_hit(muonhit_layer[l].at(hl).first->getPosition()[0], muonhit_layer[l].at(hl).first->getPosition()[1], muonhit_layer[l].at(hl).first->getPosition()[2]);
	lcio::long64 val2 = lcio::long64(muonhit_layer[l].at(hl).first->getCellID0() & 0xffffffff ) |  ( lcio::long64(muonhit_layer[l].at(hl).first->getCellID1() ) << 32) ;
	b.setValue(val2);
	int side_allhit = b["side"];
	unsigned int module_allhit = b["module"];
	float dR_muonhit = all_hit.DeltaR(starting_hit);
	
	if (dR_muonhit< m_cutdRMuon && muonhit_layer[l].at(hl).second>0 && side_allhit == side_startinghit && (abs(static_cast<int>(module_allhit-module_startinghit))%mod_module)<(m_deltaModule+1)){  
	  tmp2.push_back(hl);
	  n_tmp[l]++;
	  flag_hitlayer=1; //here we store temporarly the used hits because we check later the number of hit layers 
	}	
      } //end loop hit on l-th layer
      tmp[l]=tmp2;
      if(flag_hitlayer==1) layerhit_count++;
    }//end loop layer
    //We check the number of hit layers and store only the hits in a track with at least m_cutLayer
    if(layerhit_count<m_cutLayer) continue;
    n_SAmuons++;
    //We add hits in the cluster and evaluate the centroid 
    for(int l=0; l<n_layers; l++){
      for(int thl=0; thl<n_tmp[l]; thl++) {
	muonhit_layer[l].at(tmp[l][thl]).second=0; 
	cl->addHit(muonhit_layer[l].at(tmp[l][thl]).first,muonhit_layer[l].at(tmp[l][thl]).first->getEnergy());
	cl_energy+=muonhit_layer[l].at(tmp[l][thl]).first->getEnergy();
	for (int k=0; k<3; k++){
	  cl_centroid[k]+=muonhit_layer[l].at(tmp[l][thl]).first->getPosition()[k]*muonhit_layer[l].at(tmp[l][thl]).first->getEnergy();
	}
      }
    }
    for(int k=0; k<3; k++){
      cl_centroid[k]=cl_centroid[k]/cl_energy;
    }
    cl->setPosition(cl_centroid);
    TVector3 cl_center(cl_centroid);
    cl->setIPhi(cl_center.Phi());
    cl->setITheta(cl_center.Theta());
    cl->setEnergy(cl_energy);
    cluster->addElement(cl);
    
    //--------------------------------------------------------------------------------------------
    // For each standalone muon track we created clusters with Vertex Inner and Outer Tracker hits
    //--------------------------------------------------------------------------------------------
    if(m_trk_flag>0){
    for (unsigned int icol=0; icol<inputHitColls.size(); ++icol){
      LCCollection* hit_col  =  inputHitColls[icol];
      if( !hit_col ) continue ;
      for (int ihit=0; ihit<hit_col->getNumberOfElements(); ++ihit){
	TrackerHitPlane* hitplane = dynamic_cast<TrackerHitPlane*>(hit_col->getElementAt(ihit));
	
	TVector3 trackerhit(hitplane->getPosition()[0],hitplane->getPosition()[1],hitplane->getPosition()[2] );
	//	std::cout << "eccoci in posizione: "<<hitplane->getPosition()[0] <<std::endl;
	
	float deltaR_trackerhit=trackerhit.DeltaR(cl_center);
	if ( deltaR_trackerhit < m_cutdRTracker){
	  TrackerHitPlaneImpl* hitplane_new = new TrackerHitPlaneImpl();
	  hitplane_new->setCellID0(hitplane->getCellID0());
	  hitplane_new->setCellID1(hitplane->getCellID1());
	  hitplane_new->setPosition(hitplane->getPosition());
	  hitplane_new->setEDep(hitplane->getEDep());
	  hitplane_new->setTime(hitplane->getTime());
	  //hitplane_new->setMomentum(hitplane->getMomentum());
	  //hitplane_new->setMCParticle(hitplane->getMCParticle());
	  //hitplane_new->setPathLength(hitplane->getPathLength());
	  int mask=0x1111;
	  int tracknumber=n_SAmuons;
	  hitplane_new->setQuality(mask&tracknumber);
	  hitplane_new->setU(hitplane->getU());
	  hitplane_new->setV(hitplane->getV());
	  hitplane_new->setdU(hitplane->getdU());
	  hitplane_new->setdV(hitplane->getdV());
	  //hitplane_new->setOverlay(hitplane->isOverlay());
	  // hitplane_new->setProducedBySecondary(hitplane->isProducedBySecondary());
	  outputTrackerHitColls[icol]->addElement(hitplane_new);
	  }
      } // ihit loop      
    } // icol loop
    } //close if trk_flag
  } //end loop hits on first layer
  
  //-------------------------------
  // New collections in the event 
  //-------------------------------
  if(m_trk_flag>0){
  for(unsigned int icol=0; icol<inputHitColls.size(); ++icol){
    evt->addCollection( outputTrackerHitColls[icol], m_outputTrackerHitPlaneCollNames[icol] ) ;
    }
  }
  evt->addCollection(cluster ,m_outputClusterCollection);
} //end void processEvent

void MuonRecoStandAlone::check(LCEvent*){
}

void MuonRecoStandAlone::end(){
std::cout << "MuonRecoStandAlone::end()  " << name() 
   	    << " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
     	    << std::endl ;
}
