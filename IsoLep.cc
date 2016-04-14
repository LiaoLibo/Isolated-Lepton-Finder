#include <IsoLep.hh>
#include <algorithm>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <marlin/VerbosityLevels.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		

#include <cmath>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace std; 

IsoLep a_IsoLep_instance;

IsoLep::IsoLep()
	: Processor("IsoLep")
{
	_description = "Print MC Truth" ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"InputPFOsCollection" ,
			"Input collection of ReconstructedParticles",
			_inputPFOsCollection,
			std::string("ArborPFOs"));

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"OutputCollectionWithoutIsolatedLepton",
			"Copy of input collection but without the isolated leptons",
			_outputPFOsRemovedIsoLepCollection,
			std::string("ArborPFOsWithoutIsoLep") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"OutputCollectionIsolatedLeptons",
			"Output collection of isolated leptons",
			_outputIsoLepCollection,
			std::string("Isolep") );

	_mu_5_10_cone =0.2;
	registerProcessorParameter( "ConeAngle_of_Muon_from_5_10" ,
			"A reasonable cone angle of isolated Muon between 5GeV and 10GeV" ,
			_mu_5_10_cone ,
			_mu_5_10_cone);
	
	_mu_5_10_CorE =0.2;
	registerProcessorParameter( "ConeEnergyorMuonEnergy_5_10" ,
			"Cone Energy or Muon's Energy between 5GeV and 10GeV" ,
			_mu_5_10_CorE ,
			_mu_5_10_CorE);

	_mu_10_15_cone =0.2;
	registerProcessorParameter( "ConeAngle_of_Muon_from_10_15" ,
			"A reasonable cone angle of isolated Muon between 10GeV and 15GeV" ,
			_mu_10_15_cone ,
			_mu_10_15_cone);
	
	_mu_10_15_CorE =0.2;
	registerProcessorParameter( "ConeEnergyorMuonEnergy_10_15" ,
			"Cone Energy or Muon's Energy between 10GeV and 15GeV" ,
			_mu_10_15_CorE ,
			_mu_10_15_CorE);

	_mu_15_cone =0.2;
	registerProcessorParameter( "ConeAngle_of_Muon_from_15" ,
			"A reasonable cone angle of isolated Muon bigger than 15GeV" ,
			_mu_15_cone ,
			_mu_15_cone);
	
	_mu_15_CorE =0.2;
	registerProcessorParameter( "ConeEnergyorMuonEnergy_15" ,
			"Cone Energy or Muon's Energy bigger than 15GeV" ,
			_mu_15_CorE ,
			_mu_15_CorE);

	_ele_5_10_cone =0.2;
	registerProcessorParameter( "ConeAngle_of_Electron_from_5_10" ,
			"A reasonable cone angle of isolated Electron between 5GeV and 10GeV" ,
			_ele_5_10_cone ,
			_ele_5_10_cone);
	
	_ele_5_10_CorE =0.2;
	registerProcessorParameter( "ConeEnergyorElectronEnergy_5_10" ,
			"Cone Energy or Electron's Energy between 5GeV and 10GeV" ,
			_ele_5_10_CorE ,
			_ele_5_10_CorE);

	_ele_10_15_cone =0.2;
	registerProcessorParameter( "ConeAngle_of_Electron_from_10_15" ,
			"A reasonable cone angle of isolated Electron between 10GeV and 15GeV" ,
			_ele_10_15_cone ,
			_ele_10_15_cone);
	
	_ele_10_15_CorE =0.2;
	registerProcessorParameter( "ConeEnergyorElectronEnergy_10_15" ,
			"Cone Energy or Electron's Energy between 10GeV and 15GeV" ,
			_ele_10_15_CorE ,
			_ele_10_15_CorE);

	_ele_15_cone =0.2;
	registerProcessorParameter( "ConeAngle_of_Electron_from_15" ,
			"A reasonable cone angle of isolated Electron bigger than 15GeV" ,
			_ele_15_cone ,
			_ele_15_cone);
	
	_ele_15_CorE =0.2;
	registerProcessorParameter( "ConeEnergyorElectronEnergy_15" ,
			"Cone Energy or Electron's Energy bigger than 15GeV" ,
			_ele_15_CorE ,
			_ele_15_CorE);

}

void IsoLep::init() {

	printParameters();
	
}

void IsoLep::processRunHeader( LCRunHeader* run) {
}

void IsoLep::processEvent( LCEvent * evtP ) 
{		
	if (evtP) 								
	{
		try
		{
			LCCollection* col_RecP = evtP->getCollection( _inputPFOsCollection ); 
			_nArbor = col_RecP->getNumberOfElements();

			//output collection without Isolepton;
			LCCollectionVec* otPFOsRemovedIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
			otPFOsRemovedIsoLepCol->setSubset(true) ;

			//output collection of IsoLeptons;
			LCCollectionVec* otIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
			otIsoLepCol->setSubset(true);

			_totalEnergy =0;

			std::vector<ReconstructedParticle*> Leptons;
			std::vector<ReconstructedParticle*> Non_LepFS;
			std::vector<ReconstructedParticle*> ChargedParticle;
			std::vector<ReconstructedParticle*> NeutralParticle;

			Leptons.clear();
			Non_LepFS.clear();
			ChargedParticle.clear();
			NeutralParticle.clear();

			for (int i=0; i<_nArbor; i++){
				ReconstructedParticle* a_lep = dynamic_cast<ReconstructedParticle*>(col_RecP->getElementAt(i));
				double PID = a_lep->getType();
				_totalEnergy += a_lep->getEnergy();
				if ( abs(PID) == 11 || abs(PID) ==13 )
				{
					Leptons.push_back(a_lep);
				}
				else 
				{
					otPFOsRemovedIsoLepCol->addElement(a_lep);
					Non_LepFS.push_back(a_lep);
					if (a_lep->getCharge()!=0) ChargedParticle.push_back(a_lep);
					else NeutralParticle.push_back(a_lep);
				}
			}

			_Num_Lep = Leptons.size();
			_Num_Cha = ChargedParticle.size();
			_Num_Neu = NeutralParticle.size();

			cout <<"TotalEnergy =="<<_totalEnergy<<"; "<< "Num_lep =="<<_Num_Lep<<"; "<<"Num_Charged =="<<_Num_Cha<<"; "<<"Num_Neutral =="<<_Num_Neu<<"; "<<endl;

			_energy =0;
			_lep_momentum.SetXYZ(0,0,0);
			double leptonID = 0;

			for (int i =0; i<_Num_Lep; i++){
				ReconstructedParticle* a_lep = Leptons[i];
				leptonID = a_lep->getType();
				_energy = a_lep->getEnergy();
				_lep_momentum = a_lep->getMomentum();
				_coneEnergy =0;
				_pfo_energy =0;
				_pfo_momentum.SetXYZ(0,0,0);
				double nconetrack =0;
				for ( int j =0; j<_nArbor; j++){
					ReconstructedParticle* pfos = dynamic_cast<ReconstructedParticle*>(col_RecP->getElementAt(j));
					if (a_lep == pfos) continue;
					_pfo_energy = pfos->getEnergy();
					if (_pfo_energy < 2) continue;
					_pfo_momentum = pfos->getMomentum();
					if (fabs(leptonID) == 13) // for Isolated Muon
					{
						if ( _energy >5 && _energy <10)
						{
							if (_pfo_momentum.Angle(_lep_momentum) < _mu_5_10_cone)
							{
								_coneEnergy += _pfo_energy;
								if (pfos->getCharge()!=0) nconetrack++;
							}
						}
						if ( _energy >10 && _energy <15)
						{
							if (_pfo_momentum.Angle(_lep_momentum) < _mu_10_15_cone)
							{
								_coneEnergy += _pfo_energy;
								if (pfos->getCharge()!=0) nconetrack++;
							}
						}
						if ( _energy >15)
						{
							if (_pfo_momentum.Angle(_lep_momentum) < _mu_15_cone)
							{
								_coneEnergy += _pfo_energy;
								if (pfos->getCharge()!=0) nconetrack++;
							}
						}
					}
					if (fabs(leptonID) == 11) // for Isolated Electron
					{
						if (_energy >5 && _energy <10)
						{
							if (_pfo_momentum.Angle(_lep_momentum) < _ele_5_10_cone)
							{
								_coneEnergy += _pfo_energy;
								if (pfos->getCharge()!=0) nconetrack++;
							}
						}
						if (_energy >10 && _energy <15)
						{
							if (_pfo_momentum.Angle(_lep_momentum) < _ele_10_15_cone)
							{
								_coneEnergy += _pfo_energy;
								if (pfos->getCharge()!=0) nconetrack++;
							}
						}
						if (_energy >15)
						{
							if (_pfo_momentum.Angle(_lep_momentum) < _ele_15_cone)
							{
								_coneEnergy += _pfo_energy;
								if (pfos->getCharge()!=0) nconetrack++;
							}
						}
					}
				}
				if (_energy >5 && _energy <10)
				{
					if ( (_coneEnergy/_energy < _mu_5_10_CorE && nconetrack ==0) || (_coneEnergy/_energy < _ele_5_10_CorE && nconetrack ==0)) otIsoLepCol->addElement(a_lep); //output Isolated Lepton
					else otPFOsRemovedIsoLepCol->addElement(a_lep);
				}
				else if (_energy >10 && _energy <15)
				{
					if ( (_coneEnergy/_energy < _mu_10_15_CorE && nconetrack ==0) || (_coneEnergy/_energy < _ele_10_15_CorE && nconetrack ==0)) otIsoLepCol->addElement(a_lep); //output Isolated Lepton
					else otPFOsRemovedIsoLepCol->addElement(a_lep);
				}
				else if (_energy >15)
				{
					if ( (_coneEnergy/_energy < _mu_15_CorE && nconetrack ==0) || (_coneEnergy/_energy < _ele_15_CorE && nconetrack ==0)) otIsoLepCol->addElement(a_lep); //output Isolated Lepton
					else otPFOsRemovedIsoLepCol->addElement(a_lep);
				}
				else otPFOsRemovedIsoLepCol->addElement(a_lep);
			}
			evtP->addCollection( otPFOsRemovedIsoLepCol, _outputPFOsRemovedIsoLepCollection.c_str() );
			evtP->addCollection( otIsoLepCol, _outputIsoLepCollection.c_str() );
		}
		catch (lcio::DataNotAvailableException err) { }
	}  	  
}	

void IsoLep::end()
{
}



