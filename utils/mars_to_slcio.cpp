#include <fstream>
#include <stdio.h>
#include <math.h>
#include "lcio.h"
#include "IO/LCWriter.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/MCParticleImpl.h"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/QGSP_BERT_HP.hh"
#include "Geant4/G4ParticleTable.hh"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"

using namespace lcio;

/*
Convert the internal MARS particle ID to a PDG ID 
for extraction of the mass/charge information from the PDG table 
*/
int id_to_pdg(int particle_id) {
    switch (particle_id) {
        case 1: return 2212; break;
        case 2: return 2112; break;
        case 3: return 211; break;
        case 4: return -211; break;
        case 5: return 321; break;
        case 6: return -321; break;
        case 7: return -13; break;
        case 8: return 13; break;
        case 9: return 22; break;
        case 10: return 11; break;
        case 11: return -11; break;
        case 12: return -2212; break;
        case 13: return 111; break;
        case 14: return 1000010020; break;
        case 15: return 1000010030; break;
        case 16: return 1000020030; break;
        case 17: return 1000020040; break;
        case 18: return 14; break;
        case 19: return -14; break;
        case 20: return 12; break;
        case 21: return -12; break;
        case 22: return 130; break;
        case 23: return 310; break;
        case 24: return 311; break;
        case 25: return -311; break;
        case 26: return 3122; break;
        case 27: return -3122; break;
        case 28: return 3222; break;
        case 29: return 3212; break;
        case 30: return 3112; break;
        case 31: return -2112; break;
        case 32: return 3322; break;
        case 33: return 3312; break;
        case 34: return 3334; break;
        case 35: return 5112; break;
        case 36: return 5212; break;
        case 37: return 5222; break;
        case 38: return -3322; break;
        case 39: return -5132; break;
        case 40: return -5332; break;
    }
    return 0;
}


/*
Create a MCParticle instance with the given parameters
*/
MCParticleImpl* new_particle(G4ParticleTable *pTable, 
                             int NI, int JJ, float TOFF,
                             double X, double Y, double Z,
                             float PX, float PY, float PZ) {
    int pdg = id_to_pdg(JJ);
    // Returning nothing if particle is not present in the DPG MC particle naming scheme
    if (pdg == 0) return 0;
    MCParticleImpl* p = new MCParticleImpl();
    p->setPDG(pdg);
    p->setGeneratorStatus(1);
    p->setTime(TOFF*1e9);  // s -> ns
    const double position[3] = {X*10.0, Y*10.0, Z*10.0};  // cm -> mm
    p->setVertex(position);
    const float momentum[3] = {PX, PY, PZ};
    p->setMomentum(momentum);
    // Completing information about the particle using the GEANT4 particle table
    G4ParticleDefinition* particle = pTable->FindParticle(pdg);
    if (particle) {
        p->setMass(particle->GetPDGMass()/1000.0);  // MeV -> GeV
        p->setCharge(particle->GetPDGCharge());
    } else {
        printf("ERROR: No properties found found for a particle with PDG ID: %d\n", pdg);
    }
    return p;
}


/*
Create a single event to be stored to the output file
*/
LCEventImpl*  new_event(int eventId) {
    LCEventImpl*  event = new LCEventImpl();
    event->setRunNumber(0);
    event->setEventNumber(eventId);

    return event;
}


int main(int argc, char *argv[]) {

    if (argc < 3) {
        printf("USAGE: %s <MARS input file> <SLCIO output file> [N lines/event (1000)] [N lines max (-1)] [TOFF max [ns] (1e9)] [Neutron min Ekin [GeV] (0)] [N decays (5e5)] [W nominator (427780)]\n", argv[0]);
        return 1;
    }

    // Parsing the input arguments
    char *file_in = argv[1];
    char *file_out = argv[2];
    const int nLinesPerEvent = argc > 3 ? atoi(argv[3]) : 1000;
    int nLines_max = argc > 4 ? atoi(argv[4]) : -1;
    if (nLines_max == -1) nLines_max = 1e9;

    float toff_max = argc > 5 ? atof(argv[5]) : 1e9;
    float n_ekin_min = argc > 6 ? atof(argv[6]) : 0.0;
    
    float n_decay = argc > 7 ? atof(argv[7]) : 5e5;
    float W_nom = argc > 8 ? atof(argv[8]) : 427780;

    printf("Converting MARS data from file\n  %s\nto SLCIO file\n  %s\nwith N decays: %.0f\n", file_in, file_out, n_decay);
    printf("Reading max %d lines split by %d lines/event\n", nLines_max, nLinesPerEvent);
    
    // Calculating the correction to particle weights (from the original MARS2ROOT conversion script)
    float W_corr = W_nom/n_decay;

    // Initializing the LCIO writer
    LCWriter* lcWriter = LCFactory::getInstance()->createLCWriter();
    lcWriter->open(file_out, LCIO::WRITE_NEW);

    // Initializing the text reader of MARS output
    std::ifstream sIn;
    int nLines=0, nLinesEvent=0, iEvent=0;
    LCEventImpl*  event = new_event(iEvent);
    LCCollectionVec* vParticles = new LCCollectionVec(LCIO::MCPARTICLE);

    float Wp,Wp2;
    int NI,JJ;
    float X,Y,Z,PX,PY,PZ,TOFF,W,ZDEC,XORIG,YORIG,ZORIG,WORIG,EORIG,IORIG,KORIG;

    // Creating the particle table to complete particle information based on PDG ID
    G4RunManager* g4Manager = new G4RunManager();
    G4VUserPhysicsList* g4Physics = new QGSP_BERT_HP();
    g4Manager->SetUserInitialization(g4Physics);
    G4ParticleTable *g4Table = G4ParticleTable::GetParticleTable();

    // Creating the random number generator for adding weighted particle copies
    TRandom3* RNDM = new TRandom3(0);

    sIn.open(file_in);
    while(sIn.good()) {
        sIn >> NI >> JJ >> X >> Y >> Z >> PX >> PY >> PZ >> TOFF >> W >> ZDEC >> XORIG >> YORIG >> ZORIG >> WORIG >> EORIG >> IORIG >> KORIG;
        if (!sIn.good()) break;
        // Skipping particle if its time is outside the time region of interest
        if (TOFF*1e9 > toff_max) continue;
        // Skipping the particle if it is a neutron with too low kinetic energy (its hits will have too large time offset)
        if (JJ == 2 && sqrt(PX*PX + PY*PY + PZ*PZ) < n_ekin_min) continue;
        // Correcting the particle's weight
        W *= W_corr;
        // Testing if the particle can be created
        MCParticleImpl* p_ref = new_particle(g4Table, NI, JJ, TOFF, X,Y,Z, PX,PY,PZ);
        if (!p_ref) continue;

        ////////////////// Creating W copies of the particle randomly distributed in Phi
        // Converting the float weight to an integer with randomized +1 accounting for the fractional weight
        int nW = floor(W);
        float modint;
        if (RNDM->Rndm() < modf(W, &modint)) nW += 1;

        // Generating a random Phi rotation for each copy of the particle
        float rndm[nW];
        RNDM->RndmArray(nW, rndm);
        for (int iW=0; iW < nW; iW++) {
            TVector3 pos(X, Y, Z);
            TVector3 mom(PX, PY, PZ);
            float dPhi = rndm[iW] * TMath::TwoPi();
            pos.RotateZ(dPhi);
            mom.RotateZ(dPhi);
            MCParticleImpl* p = new_particle(g4Table, NI, JJ, TOFF, 
                                             pos.X(), pos.Y(), pos.Z(), 
                                             mom.X(), mom.Y(), mom.Z());
            vParticles->push_back(p);
        }
        // Updating the loop parameters
        nLines++;
        nLinesEvent++;
        if (nLinesEvent >= std::min(nLinesPerEvent, nLines_max)) {
            event->addCollection(vParticles, "MCParticle");
            lcWriter->writeEvent(event);
            delete vParticles;
            printf("Wrote event: %d\n", iEvent);
            // Creating the next event
            iEvent++;
            event = new_event(iEvent);
            vParticles = new LCCollectionVec(LCIO::MCPARTICLE);
            nLinesEvent = 0;
        }
        if (nLines_max > 0 && nLines >= nLines_max) break;
    }
    sIn.close();
    // Writing the final event
    if (nLinesEvent > 0) {
        event->addCollection(vParticles, "MCParticle");
        lcWriter->writeEvent(event);
	iEvent++;
    }

    lcWriter->close();

    printf("Finished writing %d events to file:\n  %s\n", iEvent, file_out);

    return 0;
}
