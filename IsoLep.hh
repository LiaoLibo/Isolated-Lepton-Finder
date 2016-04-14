#ifndef _IsoLep_hh_
#define _IsoLep_hh_

//#include <RConfigure.h>
#include <string>
#include <iostream>
#include <fstream>
#include "lcio.h"
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <TObject.h>
#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class IsoLep  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new IsoLep ; }

		IsoLep();

		~IsoLep() {};

		void init();

	    void processRunHeader( LCRunHeader* run); 

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		
		std::string _outputPFOsRemovedIsoLepCollection;
		std::string _outputIsoLepCollection;
		std::string _inputPFOsCollection;
		int _nArbor;
		double _mu_5_10_CorE;
		double _mu_5_10_cone;
		double _mu_10_15_CorE;
		double _mu_10_15_cone;
		double _mu_15_CorE;
		double _mu_15_cone;
		double _ele_5_10_CorE;
		double _ele_5_10_cone;
		double _ele_10_15_CorE;
		double _ele_10_15_cone;
		double _ele_15_CorE;
		double _ele_15_cone;
		double _IsoEnergy;
		double _coneEnergy;
		double _energy;
		double _pfo_energy;
		TVector3 _lep_momentum;
		TVector3 _pfo_momentum;

		int _Num_Lep, _Num_Cha, _Num_Neu;

		double _totalEnergy;
};

#endif
