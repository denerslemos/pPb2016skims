#include "call_libraries.h" // call libraries from ROOT and C++
#include "uiclogo.h"	    // call UIC logo and initialization
#include "CATree.h" 	    // call re-cluster for WTA axis: see https://github.com/FHead/PhysicsMiniProjects/tree/master/JetSmallSystem/24622_Recluster
#include "ntrkoff.h"        // get Ntrk offline

/*
Main skim pPb data and MC

Written by Dener Lemos (dener.lemos@cern.ch)

--> Arguments
input_file: text file with a list of root input files: Forest or Skims
ouputfile: just a counting number to run on Condor
isMC: 0 for false --> data and > 0 for true --> MC
ntrkoff: 0 for no cut/selection or MC, 1 for MB [10,185], 2 for HM PD 1 to 6 [185,250] and 3 for HM PD 7 [250, inf]
*/
void pPbSkim(TString input_file, TString ouputfile, int isMC, int ntrkoff){

	bool is_MC; if(isMC == 0){is_MC = false;}else{is_MC = true;}

	float jetptmin = 30.0;
	float jetetamin = 4.0;

	if(is_MC){
		jetptmin = 0.0;
		jetetamin = 10000.;
	}

	TString outputFileName;
	outputFileName = Form("%s",ouputfile.Data());

	clock_t sec_start, sec_end;
	sec_start = clock(); // start timing measurement

	TDatime* date = new TDatime();

	printwelcome(true); // welcome message

	print_start(); // start timing print

	// Read the input file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << input_file.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(Form("root://osg-se.sprace.org.br/%s",file_chain.c_str()));}
	inputfile.close();
	// Maximum size of arrays
	const Int_t nMaxJet = 200;				// Maximum number of jets in an event
	const Int_t nMaxTrack = 2000;		// Maximum number of tracks in an event
	
	// Define trees to be read from the files
	const int nJetTrees = 4;
	TChain *heavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
	TChain *hltTree = new TChain("hltanalysis/HltTree");
	TChain *skimTree = new TChain("skimanalysis/HltTree");
	TChain *jetTree[nJetTrees];
	jetTree[0] = new TChain("ak4CaloJetAnalyzer/t");
	jetTree[1] = new TChain("ak4PFJetAnalyzer/t");
	jetTree[2] = new TChain("akCs4PFJetAnalyzer/t");
	jetTree[3] = new TChain("ak3PFJetAnalyzer/t");
	TChain *trackTree = new TChain("ppTrack/trackTree");
	TChain *genTrackTree;
	if(is_MC){genTrackTree = new TChain("HiGenParticleAna/hi");}
	TChain *particleFlowCandidateTree = new TChain("pfcandAnalyzer/pfTree");
	TChain *checkFlatteningTree = new TChain("checkflattening/tree");

	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		hltTree->Add(*listIterator);
		trackTree->Add(*listIterator);
		heavyIonTree->Add(*listIterator);
		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){jetTree[iJetType]->Add(*listIterator);}
		skimTree->Add(*listIterator);
		if(is_MC){genTrackTree->Add(*listIterator);}
		particleFlowCandidateTree->Add(*listIterator);
		checkFlatteningTree->Add(*listIterator);
	}
	
	// All the branches and leaves come in arrat of two, one for input and one for output
	
	// Branches for heavy ion tree
	TBranch *runBranch;						 	// Branch for run
	TBranch *eventBranch;					 	// Branch for event
	TBranch *lumiBranch;						// Branch for lumi
	TBranch *hiVzBranch;						// Branch for vertex z-position
	TBranch *hiHFplusBranch;					// Branch for HF+ energy deposity
	TBranch *hiHFminusBranch;					// Branch for HF- energy deposity
	TBranch *hiZDCplusBranch;					// Branch for ZDC+ energy deposity
	TBranch *hiZDCminusBranch;					// Branch for ZDC- energy deposity
	TBranch *ptHatBranch;						// Branch for pT hat
	TBranch *eventWeightBranch;		 			// Branch for pthat weight for MC

	// Leaves for heavy ion tree
	UInt_t run;					 // Run number
	ULong64_t event;			 // Event number
	UInt_t lumi;				 // Luminosity block
	Float_t vertexZ;			 // Vertex z-position
	Float_t hiHFplus;			 // transverse energy sum of HF+ tower;
	Float_t hiHFminus;			 // transverse energy sum of HF- tower;
	Float_t hiZDCplus;			 // energy deposit in ZDC+;
	Float_t hiZDCminus;			 // energy deposit in ZDC-;
	Float_t ptHat;				 // pT hat
	Float_t eventWeight;			 // jet weight in the tree

	// Branches for EP tree
	TBranch *eventPlaneAngleBranch;						// Branch for event plane angles
	TBranch *eventPLaneQBranch;							// Branch for Q-vector magnitude in an event plane
	TBranch *eventPlaneQxBranch;						// Branch for Q-vector x-component in an event plane
	TBranch *eventPlaneQyBranch;						// Branch for Q-vector y-component in an event plane
	TBranch *eventPlaneMultiplicityBranch;	 			// Branch for event plane multiplicity
	
	const int numberofEPleaves = 182;				    		// Event plane leaves	
	// Name of EP in pPb
	TString EPNames = "HFm1/D:HFp1/D:trackmid1/D:trackm1/D:trackp1/D:trackm122/D:trackm118/D:trackm114/D:trackm110/D:trackm106/D:trackm102/D:trackp102/D:trackp106/D:trackp110/D:trackp114/D:trackp118/D:trackp122/D:trackmid1mc/D:trackm1mc/D:trackp1mc/D:trackm122mc/D:trackm118mc/D:trackm114mc/D:trackm110mc/D:trackm106mc/D:trackm102mc/D:trackp102mc/D:trackp106mc/D:trackp110mc/D:trackp114mc/D:trackp118mc/D:trackp122mc/D:HFm1a/D:HFm1b/D:HFm1c/D:HFm1d/D:HFm1e/D:HFm1f/D:HFp1a/D:HFp1b/D:HFp1c/D:HFp1d/D:HFp1e/D:HFp1f/D:HFm2/D:HFp2/D:trackmid2/D:trackm2/D:trackp2/D:trackm222/D:trackm218/D:trackm214/D:trackm210/D:trackm206/D:trackm202/D:trackp202/D:trackp206/D:trackp210/D:trackp214/D:trackp218/D:trackp222/D:HFm2a/D:HFm2b/D:HFm2c/D:HFm2d/D:HFm2e/D:HFm2f/D:HFp2a/D:HFp2b/D:HFp2c/D:HFp2d/D:HFp2e/D:HFp2f/D:HFm3/D:HFp3/D:trackmid3/D:trackm3/D:trackp3/D:trackm322/D:trackm318/D:trackm314/D:trackm310/D:trackm306/D:trackm302/D:trackp302/D:trackp306/D:trackp310/D:trackp314/D:trackp318/D:trackp322/D:HFm3a/D:HFm3b/D:HFm3c/D:HFm3d/D:HFm3e/D:HFm3f/D:HFp3a/D:HFp3b/D:HFp3c/D:HFp3d/D:HFp3e/D:HFp3f/D:HFm4/D:HFp4/D:trackmid4/D:trackm4/D:trackp4/D:trackm422/D:trackm418/D:trackm414/D:trackm410/D:trackm406/D:trackm402/D:trackp402/D:trackp406/D:trackp410/D:trackp414/D:trackp418/D:trackp422/D:HFm4a/D:HFm4b/D:HFm4c/D:HFm43d/D:HFm4e/D:HFm4f/D:HFp4a/D:HFp4b/D:HFp4c/D:HFp4d/D:HFp4e/D:HFp4f/D:HFm5/D:HFp5/D:trackmid5/D:trackm5/D:trackp5/D:trackm522/D:trackm518/D:trackm514/D:trackm510/D:trackm506/D:trackm502/D:trackp502/D:trackp506/D:trackp510/D:trackp514/D:trackp518/D:trackp522/D:HFm6/D:HFp6/D:trackmid6/D:trackm6/D:trackp6/D:trackm622/D:trackm618/D:trackm614/D:trackm610/D:trackm606/D:trackm602/D:trackp602/D:trackp606/D:trackp610/D:trackp614/D:trackp618/D:trackp622/D:HFm7/D:HFp7/D:trackmid7/D:trackm7/D:trackp7/D:trackm722/D:trackm718/D:trackm714/D:trackm710/D:trackm706/D:trackm702/D:trackp702/D:trackp706/D:trackp710/D:trackp714/D:trackp718/D:trackp722/D"; 
	Double_t eventPlaneAngle[numberofEPleaves] = {0};			// Event plane angles
	Double_t eventPlaneQ[numberofEPleaves] = {0};				// Magnitude of Q-vector in event plane
	Double_t eventPlaneQx[numberofEPleaves] = {0};				// x-component of the Q-vector
	Double_t eventPlaneQy[numberofEPleaves] = {0};				// y-component of the Q-vector
	Double_t eventPlaneMultiplicity[numberofEPleaves] = {0};	// Particle multiplicity in an event plane

	// Branches for HLT tree
	// HLT
	TBranch *caloJetFilterBranch60;				 		// Branch for calo jet 60 filter bit
	TBranch *caloJetFilterBranch80;				 		// Branch for calo jet 80 filter bit
	TBranch *caloJetFilterBranch100;			 		// Branch for calo jet 100 filter bit
	TBranch *pfJetFilterBranch60;				 		// Branch for PF jet 60 filter bit
	TBranch *pfJetFilterBranch80;				 		// Branch for PF jet 80 filter bit
	TBranch *pfJetFilterBranch100;						// Branch for PF jet 100 filter bit
	TBranch *pfJetFilterBranch120;						// Branch for PF jet 120 filter bit
	
	// Leaves for the HLT tree
	Int_t caloJetFilterBit60;					// Filter bit for calorimeter jets 60
	Int_t caloJetFilterBit80;					// Filter bit for calorimeter jets 80
	Int_t caloJetFilterBit100;					// Filter bit for calorimeter jets 100
	Int_t pfJetFilterBit60;						// Filter bit for particle flow flow jets 60
	Int_t pfJetFilterBit80;						// Filter bit for particle flow jets 80
	Int_t pfJetFilterBit100;					// Filter bit for particle flow jets 100
	Int_t pfJetFilterBit120;					// Filter bit for particle flow jets 100

	// Branches for skim tree
	TBranch *primaryVertexBranch;						// Branch for primary vertex filter bit
	TBranch *beamScrapingBranch;				// Branch for beam scraping filter bit
	TBranch *hBHENoiseBranchLoose;				// Branch for HB/HE noise filter bit loose
	TBranch *hBHENoiseBranchTight;				// Branch for HB/HE noise filter bit tight
	TBranch *hfCoincidenceBranch;				// Branch for energy recorded one HF tower above threshold on each side
	TBranch *pVertexFilterCutdz1p0Branch;		// Branch for PU Filter default
	TBranch *pVertexFilterCutGplusBranch;		// Branch for PU Filter GPlus
	TBranch *pVertexFilterCutVtx1Branch;		// Branch for PU Filter 1 vertex only

	// Leaves for the skim tree
	Int_t primaryVertexFilterBit;				// Filter bit for primary vertex
	Int_t beamScrapingFilterBit;				// Filter bit for beam scraping
	Int_t hBHENoiseFilterLooseBit;	 			// Filter bit for HB/HE noise loose
	Int_t hBHENoiseFilterTightBit; 				// Filter bit for HB/HE noise tight
	Int_t hfCoincidenceFilterBit;				// Filter bit or energy recorded one HF tower above threshold on each side
	Int_t pVertexFilterCutdz1p0Bit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutGplusBit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutVtx1Bit;					// Filter bit for PU Filter
	
	// Branches for jet tree
	TBranch *nJetsBranch[nJetTrees];				// Branch for number of jets in an event
	TBranch *jetRawPtBranch[nJetTrees];				// Branch for raw jet pT
	TBranch *jetMaxTrackPtBranch[nJetTrees];		// Maximum pT for a track inside a jet
	TBranch *jetPhiBranch[nJetTrees];				// Branch for jet phi
	TBranch *jetPhiBranchWTA[nJetTrees];			// Branch for jet phi with WTA axis
	TBranch *jetEtaBranch[nJetTrees];				// Branch for jet eta
	TBranch *jetEtaBranchWTA[nJetTrees];			// Branch for jet eta with WTA axis

	TBranch *jetRefPtBranch[nJetTrees];				// Branch for reference generator level pT for a reconstructed jet
	TBranch *jetRefEtaBranch[nJetTrees];			// Branch for reference generator level eta for a reconstructed jet
	TBranch *jetRefPhiBranch[nJetTrees];			// Branch for reference generator level phi for a reconstructed jet
	TBranch *jetRefFlavorBranch[nJetTrees];			// Branch for flavor for the parton initiating the jet
	TBranch *jetRefFlavorForBBranch[nJetTrees];		// Branch for flavor for the parton initiating the jet
	TBranch *jetRefSubidBranch[nJetTrees];		    // Branch for jet subid

	TBranch *nGenJetsBranch[nJetTrees];				// Branch for the number of generator level jets in an event
	TBranch *genJetPtBranch[nJetTrees];				// Branch for the generator level jet pT
	TBranch *genJetEtaBranch[nJetTrees];			// Branch for the generetor level jet eta
	TBranch *genJetEtaBranchWTA[nJetTrees];			// Branch for the generetor level jet eta with WTA axis
	TBranch *genJetPhiBranch[nJetTrees];			// Branch for the generator level jet phi
	TBranch *genJetPhiBranchWTA[nJetTrees];			// Branch for the generator level jet phi with WTA axis
	TBranch *genJetSubidBranch[nJetTrees];                // Branch for the generator level jet subid
	TBranch *genJetMatchIndexBranch[nJetTrees];           // Branch for the generator level jet matched index
	
	// Leaves for jet tree
	Int_t nJets[nJetTrees];									// number of jets in an event
	Float_t jetRawPtArray[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the jets in an event
	Float_t jetMaxTrackPtArray[nJetTrees][nMaxJet] = {{0}}; // maximum track pT inside a jet for all the jets in an event
	Float_t jetPhiArray[nJetTrees][nMaxJet] = {{0}};		// phi of all the jets in an event
	Float_t jetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};		// phi of all the jets in an event	with WTA axis
	Float_t jetEtaArray[nJetTrees][nMaxJet] = {{0}};		// eta of all the jets in an event
	Float_t jetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};		// eta of all the jets in an event	with WTA axis

	Float_t jetRefPtArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Float_t jetRefEtaArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Float_t jetRefPhiArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Int_t jetRefFlavorArray[nJetTrees][nMaxJet] = {{0}};	// flavor for initiating parton for the reference gen jet
	Int_t jetRefFlavorForBArray[nJetTrees][nMaxJet] = {{0}};// heavy flavor for initiating parton for the reference gen jet
	Int_t jetRefSubidArray[nJetTrees][nMaxJet] = {{0}};     // jet subid

	Int_t nGenJets[nJetTrees];								// number of generator level jets in an event
	Float_t genJetPtArray[nJetTrees][nMaxJet] = {{0}};		// pT of all the generator level jets in an event
	Float_t genJetPhiArray[nJetTrees][nMaxJet] = {{0}};		// phi of all the generator level jets in an event
	Float_t genJetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};	// phi of all the generator level jets in an event with WTA axis
	Float_t genJetEtaArray[nJetTrees][nMaxJet] = {{0}};		// eta of all the generator level jets in an event
	Float_t genJetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};	// eta of all the generator level jets in an event with WTA axis
	Int_t genJetSubidArray[nJetTrees][nMaxJet] = {{0}};     // subid of all the generator level jets in an event
	Int_t genJetMatchIndexArray[nJetTrees][nMaxJet] = {{0}};// matched index of all the generator level jets in an event
	
	// Branches for track tree
	TBranch *nTracksBranch;									// Branch for number of tracks
	TBranch *trackPtBranch;									// Branch for track pT
	TBranch *trackPtErrorBranch;							// Branch for track pT error
	TBranch *trackPhiBranch;								// Branch for track phi
	TBranch *trackEtaBranch;								// Branch for track eta
	TBranch *trackHighPurityBranch;							// Branch for high purity of the track
	TBranch *trackVertexDistanceZBranch;			 		// Branch for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceZErrorBranch;				// Branch for error for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceXYBranch;					// Branch for track distance from primary vertex in xy-direction
	TBranch *trackVertexDistanceXYErrorBranch; 				// Branch for error for track distance from primary vertex in xy-direction
	TBranch *trackEnergyEcalBranch;							// Branch for track energy in ECal
	TBranch *trackEnergyHcalBranch;							// Branch for track energy in HCal
	TBranch *trackChargeBranch;								// Branch for track charge
	TBranch *PixelnHitsTrackBranch;							// Branch for number of valid pixel hits for the track
	
	// Leaves for the track tree
	Int_t nTracks;														// Number of tracks
	Float_t trackPtArray[nMaxTrack] = {0};								// Array for track pT
	Float_t trackPtErrorArray[nMaxTrack] = {0};							// Array for track pT errors
	Float_t trackPhiArray[nMaxTrack] = {0};								// Array for track phis
	Float_t trackEtaArray[nMaxTrack] = {0};								// Array for track etas
	Bool_t trackHighPurityArray[nMaxTrack] = {0};						// Array for the high purity of tracks
	Float_t trackVertexDistanceZArray[nMaxTrack] = {0};			 		// Array for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceZErrorArray[nMaxTrack] = {0};			// Array for error for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceXYArray[nMaxTrack] = {0};				// Array for track distance from primary vertex in xy-direction
	Float_t trackVertexDistanceXYErrorArray[nMaxTrack] = {0}; 			// Array for error for track distance from primary vertex in xy-direction
	Float_t trackEnergyEcalArray[nMaxTrack] = {0};						// Array for track energy in ECal
	Float_t trackEnergyHcalArray[nMaxTrack] = {0};						// Array for track energy in HCal
	Int_t trackChargeArray[nMaxTrack] = {0}; 										// Array for track charge
	UChar_t PixelnHitsTrackArray[nMaxTrack] = {0}; 							// Array for number of valid pixel hits for the track

	// Branches for generator level track tree
	TBranch *genTrackPtBranch;				 // Branch for generator level track pT:s
	TBranch *genTrackPhiBranch;				 // Branch for generator level track phis
	TBranch *genTrackEtaBranch;				 // Branch for generator level track etas
	TBranch *genTrackPdgBranch;				 // Branch for generator level track PDG code
	TBranch *genTrackChargeBranch;				 // Branch for generator level track charges
	TBranch *genTrackSubeventBranch;			 // Branch for generator level track subevent indices (0 = PYTHIA, (>0) = other MC)
	
	// Leaves for generator level track tree
	vector<float> *genTrackPtArray;			 // Array for generator level track pT
	vector<float> *genTrackPhiArray;		 // Array for generator level track phi
	vector<float> *genTrackEtaArray;		 // Array for generator level track eta
	vector<int>	 *genTrackPdgArray;		 // Array for generator level track PDG code
	vector<int>	 *genTrackChargeArray;		 // Array for generator level track charges
	vector<int>	 *genTrackSubeventArray;	 // Array for generator level track subevent indices (0 = PYTHIA, (>0) = other MC)
	
	// Branches for particle flow candidate ID tree
	TBranch *particleFlowCandidatePtBranch;		// Branch for particle flow candidate pT
	TBranch *particleFlowCandidateEtaBranch;	 // Branch for particle flow candidate eta
	TBranch *particleFlowCandidatePhiBranch;	 // Branch for particle flow candidate phi
	
	// Leaves for particle flow candidate tree
	vector<float> *particleFlowCandidatePtVector;		 // Vector for particle flow candidate pT
	vector<float> *particleFlowCandidateEtaVector;		// Vector for particle flow candidate eta
	vector<float> *particleFlowCandidatePhiVector;		// Vector for particle flow candidate phi

	// ========================================== //
	// Read all the branches from the input trees //
	// ========================================== //
	
	// Connect the branches of the heavy ion tree
	heavyIonTree->SetBranchStatus("*",0); // remove all branchs to read it fast
	heavyIonTree->SetBranchStatus("run",1);
	heavyIonTree->SetBranchAddress("run",&run,&runBranch);
	heavyIonTree->SetBranchStatus("evt",1);
	heavyIonTree->SetBranchAddress("evt",&event,&eventBranch);
	heavyIonTree->SetBranchStatus("lumi",1);
	heavyIonTree->SetBranchAddress("lumi",&lumi,&lumiBranch);
	heavyIonTree->SetBranchStatus("vz",1);
	heavyIonTree->SetBranchAddress("vz",&vertexZ,&hiVzBranch);
	heavyIonTree->SetBranchStatus("hiHFplus",1);
	heavyIonTree->SetBranchAddress("hiHFplus",&hiHFplus,&hiHFplusBranch);
	heavyIonTree->SetBranchStatus("hiHFminus",1);
	heavyIonTree->SetBranchAddress("hiHFminus",&hiHFminus,&hiHFminusBranch);
	heavyIonTree->SetBranchStatus("hiZDCplus",1);
	heavyIonTree->SetBranchAddress("hiZDCplus",&hiZDCplus,&hiZDCplusBranch);
	heavyIonTree->SetBranchStatus("hiZDCminus",1);
	heavyIonTree->SetBranchAddress("hiZDCminus",&hiZDCminus,&hiZDCminusBranch);
	
	// Event plane
	checkFlatteningTree->SetBranchStatus("epang",1);
	checkFlatteningTree->SetBranchAddress("epang",&eventPlaneAngle,&eventPlaneAngleBranch);
	checkFlatteningTree->SetBranchStatus("q",1);
	checkFlatteningTree->SetBranchAddress("q",&eventPlaneQ,&eventPLaneQBranch);
	checkFlatteningTree->SetBranchStatus("qx",1);
	checkFlatteningTree->SetBranchAddress("qx",&eventPlaneQx,&eventPlaneQxBranch);
	checkFlatteningTree->SetBranchStatus("qy",1);
	checkFlatteningTree->SetBranchAddress("qy",&eventPlaneQy,&eventPlaneQyBranch);
	checkFlatteningTree->SetBranchStatus("mult",1);
	checkFlatteningTree->SetBranchAddress("mult",&eventPlaneMultiplicity,&eventPlaneMultiplicityBranch);

	// ptHat and event weight only for MC
	if(is_MC){
		heavyIonTree->SetBranchStatus("pthat",1);
		heavyIonTree->SetBranchAddress("pthat",&ptHat,&ptHatBranch);
		heavyIonTree->SetBranchStatus("weight",1);
		heavyIonTree->SetBranchAddress("weight",&eventWeight,&eventWeightBranch);
	}
	
	// Connect the branches to the HLT tree
	hltTree->SetBranchStatus("*",0);
	
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet60_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v3",&caloJetFilterBit60,&caloJetFilterBranch60);
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet80_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v3",&caloJetFilterBit80,&caloJetFilterBranch80);
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet100_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v3",&caloJetFilterBit100,&caloJetFilterBranch100);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet60_Eta5p1_v4",1);

	hltTree->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4",&pfJetFilterBit60,&pfJetFilterBranch60);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet80_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3",&pfJetFilterBit80,&pfJetFilterBranch80);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet100_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3",&pfJetFilterBit100,&pfJetFilterBranch100);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet120_Eta5p1_v2",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet120_Eta5p1_v2",&pfJetFilterBit120,&pfJetFilterBranch120);


	// Connect the branches to the skim tree
	skimTree->SetBranchStatus("*",0);
	skimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
	skimTree->SetBranchAddress("pPAprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
	skimTree->SetBranchStatus("pBeamScrapingFilter",1);
	skimTree->SetBranchAddress("pBeamScrapingFilter",&beamScrapingFilterBit,&beamScrapingBranch);
	skimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
	skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&hBHENoiseFilterLooseBit,&hBHENoiseBranchLoose);
	skimTree->SetBranchStatus("HBHENoiseFilterResultRun2Tight",1);
	skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Tight",&hBHENoiseFilterTightBit,&hBHENoiseBranchTight);
	skimTree->SetBranchStatus("phfCoincFilter",1);
	skimTree->SetBranchAddress("phfCoincFilter", &hfCoincidenceFilterBit, &hfCoincidenceBranch);
	skimTree->SetBranchStatus("pVertexFilterCutdz1p0",1);
	skimTree->SetBranchAddress("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, &pVertexFilterCutdz1p0Branch);
	skimTree->SetBranchStatus("pVertexFilterCutGplus",1);
	skimTree->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplusBit, &pVertexFilterCutGplusBranch);
	skimTree->SetBranchStatus("pVertexFilterCutVtx1",1);
	skimTree->SetBranchAddress("pVertexFilterCutVtx1", &pVertexFilterCutVtx1Bit, &pVertexFilterCutVtx1Branch);


	// Same branch names for all jet collections
	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
		// Connect the branches to the jet tree
		jetTree[iJetType]->SetBranchStatus("*",0);

		jetTree[iJetType]->SetBranchStatus("nref",1);
		jetTree[iJetType]->SetBranchAddress("nref",&nJets[iJetType],&nJetsBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("rawpt",1);
		jetTree[iJetType]->SetBranchAddress("rawpt",&jetRawPtArray[iJetType],&jetRawPtBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("trackMax",1);
		jetTree[iJetType]->SetBranchAddress("trackMax",&jetMaxTrackPtArray[iJetType],&jetMaxTrackPtBranch[iJetType]);
		
		// Jet phi with E-scheme, WTA axes calculated later
		jetTree[iJetType]->SetBranchStatus("jtphi",1);
		jetTree[iJetType]->SetBranchAddress("jtphi",&jetPhiArray[iJetType],&jetPhiBranch[iJetType]);
		
		// Jet eta with E-scheme, WTA axes calculated later
		jetTree[iJetType]->SetBranchStatus("jteta",1);
		jetTree[iJetType]->SetBranchAddress("jteta",&jetEtaArray[iJetType],&jetEtaBranch[iJetType]);
	
		// If we are looking at Monte Carlo, connect the reference pT and parton arrays
		if(is_MC){
			// Matched jet variables
			jetTree[iJetType]->SetBranchStatus("refpt",1);
			jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refpt",1);
			jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refpt",1);
			jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refeta",1);
			jetTree[iJetType]->SetBranchAddress("refeta",&jetRefEtaArray[iJetType],&jetRefEtaBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refphi",1);
			jetTree[iJetType]->SetBranchAddress("refphi",&jetRefPhiArray[iJetType],&jetRefPhiBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refparton_flavor",1);
			jetTree[iJetType]->SetBranchAddress("refparton_flavor",&jetRefFlavorArray[iJetType],&jetRefFlavorBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refparton_flavorForB",1);
			jetTree[iJetType]->SetBranchAddress("refparton_flavorForB", &jetRefFlavorForBArray[iJetType], &jetRefFlavorForBBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("subid",1);
			jetTree[iJetType]->SetBranchAddress("subid", &jetRefSubidArray[iJetType], &jetRefSubidBranch[iJetType]);
			
			// Gen jet variables
			jetTree[iJetType]->SetBranchStatus("ngen",1);
			jetTree[iJetType]->SetBranchAddress("ngen",&nGenJets[iJetType],&nGenJetsBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genpt",1);
			jetTree[iJetType]->SetBranchAddress("genpt",&genJetPtArray[iJetType],&genJetPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genphi",1);
			jetTree[iJetType]->SetBranchAddress("genphi",&genJetPhiArray[iJetType],&genJetPhiBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("geneta",1);
			jetTree[iJetType]->SetBranchAddress("geneta",&genJetEtaArray[iJetType],&genJetEtaBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genmatchindex",1);
			jetTree[iJetType]->SetBranchAddress("genmatchindex",&genJetMatchIndexArray[iJetType],&genJetMatchIndexBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("gensubid",1);
			jetTree[iJetType]->SetBranchAddress("gensubid",&genJetSubidArray[iJetType],&genJetSubidBranch[iJetType]);
			
		}
		
	} // Loop over different jet collections


	// Connect the branches to the track tree
	trackTree->SetBranchStatus("*",0);

	trackTree->SetBranchStatus("nTrk",1);
	trackTree->SetBranchAddress("nTrk",&nTracks,&nTracksBranch);
	trackTree->SetBranchStatus("highPurity",1);
	trackTree->SetBranchAddress("highPurity",&trackHighPurityArray,&trackHighPurityBranch);
	trackTree->SetBranchStatus("trkPt",1);
	trackTree->SetBranchAddress("trkPt",&trackPtArray,&trackPtBranch);
	trackTree->SetBranchStatus("trkPtError",1);
	trackTree->SetBranchAddress("trkPtError",&trackPtErrorArray,&trackPtErrorBranch);
	trackTree->SetBranchStatus("trkPhi",1);
	trackTree->SetBranchAddress("trkPhi",&trackPhiArray,&trackPhiBranch);
	trackTree->SetBranchStatus("trkEta",1);
	trackTree->SetBranchAddress("trkEta",&trackEtaArray,&trackEtaBranch);
	trackTree->SetBranchStatus("trkDz1",1);
	trackTree->SetBranchAddress("trkDz1",&trackVertexDistanceZArray,&trackVertexDistanceZBranch);
	trackTree->SetBranchStatus("trkDzError1",1);
	trackTree->SetBranchAddress("trkDzError1",&trackVertexDistanceZErrorArray,&trackVertexDistanceZErrorBranch);
	trackTree->SetBranchStatus("trkDxy1",1);
	trackTree->SetBranchAddress("trkDxy1",&trackVertexDistanceXYArray,&trackVertexDistanceXYBranch);
	trackTree->SetBranchStatus("trkDxyError1",1);
	trackTree->SetBranchAddress("trkDxyError1",&trackVertexDistanceXYErrorArray,&trackVertexDistanceXYErrorBranch);
	trackTree->SetBranchStatus("trkNPixelHit",1);
	trackTree->SetBranchAddress("trkNPixelHit",&PixelnHitsTrackArray,&PixelnHitsTrackBranch);
	trackTree->SetBranchStatus("pfEcal",1);
	trackTree->SetBranchAddress("pfEcal",&trackEnergyEcalArray,&trackEnergyEcalBranch);
	trackTree->SetBranchStatus("pfHcal",1);
	trackTree->SetBranchAddress("pfHcal",&trackEnergyHcalArray,&trackEnergyHcalBranch);
	trackTree->SetBranchStatus("trkCharge",1);
	trackTree->SetBranchAddress("trkCharge",&trackChargeArray,&trackChargeBranch);
	
	// Generator level tracks only in Monte Carlo
	if(is_MC){
		
		// Connect the branches to generator level track tree
		genTrackTree->SetBranchStatus("*",0);
		genTrackTree->SetBranchStatus("pt",1);
		genTrackTree->SetBranchAddress("pt",&genTrackPtArray,&genTrackPtBranch);
		genTrackTree->SetBranchStatus("phi",1);
		genTrackTree->SetBranchAddress("phi",&genTrackPhiArray,&genTrackPhiBranch);
		genTrackTree->SetBranchStatus("eta",1);
		genTrackTree->SetBranchAddress("eta",&genTrackEtaArray,&genTrackEtaBranch);
		genTrackTree->SetBranchStatus("pdg",1);
		genTrackTree->SetBranchAddress("pdg",&genTrackPdgArray,&genTrackPdgBranch);
		genTrackTree->SetBranchStatus("chg",1);
		genTrackTree->SetBranchAddress("chg",&genTrackChargeArray,&genTrackChargeBranch);
		genTrackTree->SetBranchStatus("sube",1);
		genTrackTree->SetBranchAddress("sube",&genTrackSubeventArray,&genTrackSubeventBranch);
		
	}
	
	// Connect the branches to the particle flow candidate tree if requested
	particleFlowCandidateTree->SetBranchStatus("*",0);
	particleFlowCandidateTree->SetBranchStatus("pfPt",1);
	particleFlowCandidateTree->SetBranchAddress("pfPt",&particleFlowCandidatePtVector,&particleFlowCandidatePtBranch);
	particleFlowCandidateTree->SetBranchStatus("pfPhi",1);
	particleFlowCandidateTree->SetBranchAddress("pfPhi",&particleFlowCandidatePhiVector,&particleFlowCandidatePhiBranch);
	particleFlowCandidateTree->SetBranchStatus("pfEta",1);
	particleFlowCandidateTree->SetBranchAddress("pfEta",&particleFlowCandidateEtaVector,&particleFlowCandidateEtaBranch);

	// ========================================== //
	//			 Define output trees
	// ========================================== //
	
	// Copy the heavy ion tree to the output
	TTree *heavyIonTreeOutput = new TTree("HiTree","");
	// Connect the branches of the heavy ion tree
	heavyIonTreeOutput->Branch("run",&run,"run/i");
	heavyIonTreeOutput->Branch("evt",&event,"evt/l");
	heavyIonTreeOutput->Branch("lumi",&lumi,"lumi/i");
	heavyIonTreeOutput->Branch("vz",&vertexZ,"vz/F");
	heavyIonTreeOutput->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
	heavyIonTreeOutput->Branch("hiHFminus",&hiHFminus,"hiHFminus/F");
	heavyIonTreeOutput->Branch("hiZDCplus",&hiZDCplus,"hiZDCplus/F");
	heavyIonTreeOutput->Branch("hiZDCminus",&hiZDCminus,"hiZDCminus/F");
	
	// Event plane
	TTree *checkFlatteningTreeOutput = new TTree("tree","");
	/*
	checkFlatteningTreeOutput->Branch("epang",&eventPlaneAngle,Form("%s",EPNames.Data()));
	checkFlatteningTreeOutput->Branch("q",&eventPlaneQ,Form("%s",EPNames.Data()));
	checkFlatteningTreeOutput->Branch("qx",&eventPlaneQx,Form("%s",EPNames.Data()));
	checkFlatteningTreeOutput->Branch("qy",&eventPlaneQy,Form("%s",EPNames.Data()));
	checkFlatteningTreeOutput->Branch("mult",&eventPlaneMultiplicity,Form("%s",EPNames.Data()));
	*/
	Float_t epang_HFm2, epang_HFp2,epang_HFm3, epang_HFp3, epang_HFm4, epang_HFp4, epang_HFm5, epang_HFp5, epang_HFm6, epang_HFp6;	
	Float_t q_HFm2, q_HFp2,q_HFm3, q_HFp3, q_HFm4, q_HFp4, q_HFm5, q_HFp5, q_HFm6, q_HFp6;
	Float_t qx_HFm2, qx_HFp2,qx_HFm3, qx_HFp3, qx_HFm4, qx_HFp4, qx_HFm5, qx_HFp5, qx_HFm6, qx_HFp6;
	Float_t qy_HFm2, qy_HFp2,qy_HFm3, qy_HFp3, qy_HFm4, qy_HFp4, qy_HFm5, qy_HFp5, qy_HFm6, qy_HFp6;
	Float_t mult_HFm2, mult_HFp2,mult_HFm3, mult_HFp3, mult_HFm4, mult_HFp4, mult_HFm5, mult_HFp5, mult_HFm6, mult_HFp6;

	checkFlatteningTreeOutput->Branch("epang_HFm2",&epang_HFm2,"epang_HFm2/F");
	checkFlatteningTreeOutput->Branch("epang_HFp2",&epang_HFp2,"epang_HFp2/F");
	checkFlatteningTreeOutput->Branch("epang_HFm3",&epang_HFm3,"epang_HFm3/F");
	checkFlatteningTreeOutput->Branch("epang_HFp3",&epang_HFp3,"epang_HFp3/F");
	checkFlatteningTreeOutput->Branch("epang_HFm4",&epang_HFm4,"epang_HFm4/F");
	checkFlatteningTreeOutput->Branch("epang_HFp4",&epang_HFp4,"epang_HFp4/F");
	checkFlatteningTreeOutput->Branch("epang_HFm5",&epang_HFm5,"epang_HFm5/F");
	checkFlatteningTreeOutput->Branch("epang_HFp5",&epang_HFp5,"epang_HFp5/F");
	checkFlatteningTreeOutput->Branch("epang_HFm6",&epang_HFm6,"epang_HFm6/F");
	checkFlatteningTreeOutput->Branch("epang_HFp6",&epang_HFp6,"epang_HFp6/F");

	checkFlatteningTreeOutput->Branch("q_HFm2",&q_HFm2,"q_HFm2/F");
	checkFlatteningTreeOutput->Branch("q_HFp2",&q_HFp2,"q_HFp2/F");
	checkFlatteningTreeOutput->Branch("q_HFm3",&q_HFm3,"q_HFm3/F");
	checkFlatteningTreeOutput->Branch("q_HFp3",&q_HFp3,"q_HFp3/F");
	checkFlatteningTreeOutput->Branch("q_HFm4",&q_HFm4,"q_HFm4/F");
	checkFlatteningTreeOutput->Branch("q_HFp4",&q_HFp4,"q_HFp4/F");
	checkFlatteningTreeOutput->Branch("q_HFm5",&q_HFm5,"q_HFm5/F");
	checkFlatteningTreeOutput->Branch("q_HFp5",&q_HFp5,"q_HFp5/F");
	checkFlatteningTreeOutput->Branch("q_HFm6",&q_HFm6,"q_HFm6/F");
	checkFlatteningTreeOutput->Branch("q_HFp6",&q_HFp6,"q_HFp6/F");

	checkFlatteningTreeOutput->Branch("qx_HFm2",&qx_HFm2,"qx_HFm2/F");
	checkFlatteningTreeOutput->Branch("qx_HFp2",&qx_HFp2,"qx_HFp2/F");
	checkFlatteningTreeOutput->Branch("qx_HFm3",&qx_HFm3,"qx_HFm3/F");
	checkFlatteningTreeOutput->Branch("qx_HFp3",&qx_HFp3,"qx_HFp3/F");
	checkFlatteningTreeOutput->Branch("qx_HFm4",&qx_HFm4,"qx_HFm4/F");
	checkFlatteningTreeOutput->Branch("qx_HFp4",&qx_HFp4,"qx_HFp4/F");
	checkFlatteningTreeOutput->Branch("qx_HFm5",&qx_HFm5,"qx_HFm5/F");
	checkFlatteningTreeOutput->Branch("qx_HFp5",&qx_HFp5,"qx_HFp5/F");
	checkFlatteningTreeOutput->Branch("qx_HFm6",&qx_HFm6,"qx_HFm6/F");
	checkFlatteningTreeOutput->Branch("qx_HFp6",&qx_HFp6,"qx_HFp6/F");

	checkFlatteningTreeOutput->Branch("qy_HFm2",&qy_HFm2,"qy_HFm2/F");
	checkFlatteningTreeOutput->Branch("qy_HFp2",&qy_HFp2,"qy_HFp2/F");
	checkFlatteningTreeOutput->Branch("qy_HFm3",&qy_HFm3,"qy_HFm3/F");
	checkFlatteningTreeOutput->Branch("qy_HFp3",&qy_HFp3,"qy_HFp3/F");
	checkFlatteningTreeOutput->Branch("qy_HFm4",&qy_HFm4,"qy_HFm4/F");
	checkFlatteningTreeOutput->Branch("qy_HFp4",&qy_HFp4,"qy_HFp4/F");
	checkFlatteningTreeOutput->Branch("qy_HFm5",&qy_HFm5,"qy_HFm5/F");
	checkFlatteningTreeOutput->Branch("qy_HFp5",&qy_HFp5,"qy_HFp5/F");
	checkFlatteningTreeOutput->Branch("qy_HFm6",&qy_HFm6,"qy_HFm6/F");
	checkFlatteningTreeOutput->Branch("qy_HFp6",&qy_HFp6,"qy_HFp6/F");

	checkFlatteningTreeOutput->Branch("mult_HFm2",&mult_HFm2,"mult_HFm2/F");
	checkFlatteningTreeOutput->Branch("mult_HFp2",&mult_HFp2,"mult_HFp2/F");
	checkFlatteningTreeOutput->Branch("mult_HFm3",&mult_HFm3,"mult_HFm3/F");
	checkFlatteningTreeOutput->Branch("mult_HFp3",&mult_HFp3,"mult_HFp3/F");
	checkFlatteningTreeOutput->Branch("mult_HFm4",&mult_HFm4,"mult_HFm4/F");
	checkFlatteningTreeOutput->Branch("mult_HFp4",&mult_HFp4,"mult_HFp4/F");
	checkFlatteningTreeOutput->Branch("mult_HFm5",&mult_HFm5,"mult_HFm5/F");
	checkFlatteningTreeOutput->Branch("mult_HFp5",&mult_HFp5,"mult_HFp5/F");
	checkFlatteningTreeOutput->Branch("mult_HFm6",&mult_HFm6,"mult_HFm6/F");
	checkFlatteningTreeOutput->Branch("mult_HFp6",&mult_HFp6,"mult_HFp6/F");
	
	// ptHat and event weight only for MC
	if(is_MC){
		heavyIonTreeOutput->Branch("pthat",&ptHat,"pthat/F");
		heavyIonTreeOutput->Branch("weight",&eventWeight,"weight/F");
	}
	
	// Copy the HLT tree to the output
	TTree *hltTreeOutput = new TTree("HltTree","");	
	// Connect the branches of the HLT tree
	hltTreeOutput->Branch("HLT_PAAK4CaloJet60_Eta5p1_v3",&caloJetFilterBit60,"HLT_PAAK4CaloJet60_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4CaloJet80_Eta5p1_v3",&caloJetFilterBit80,"HLT_PAAK4CaloJet80_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4CaloJet100_Eta5p1_v3",&caloJetFilterBit100,"HLT_PAAK4CaloJet100_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet60_Eta5p1_v4",&pfJetFilterBit60,"HLT_PAAK4PFJet60_Eta5p1_v4/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet80_Eta5p1_v3",&pfJetFilterBit80,"HLT_PAAK4PFJet80_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet100_Eta5p1_v3",&pfJetFilterBit100,"HLT_PAAK4PFJet100_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet120_Eta5p1_v2",&pfJetFilterBit120,"HLT_PAAK4PFJet120_Eta5p1_v2/I");

	// Copy the skim tree to the output
	TTree *skimTreeOutput = new TTree("HltTree","");
	skimTreeOutput->Branch("pPAprimaryVertexFilter",&primaryVertexFilterBit,"pPAprimaryVertexFilter/I");
	skimTreeOutput->Branch("pBeamScrapingFilter",&beamScrapingFilterBit,"pBeamScrapingFilter/I");
	skimTreeOutput->Branch("HBHENoiseFilterResultRun2Loose",&hBHENoiseFilterLooseBit,"HBHENoiseFilterResultRun2Loose/I");
	skimTreeOutput->Branch("HBHENoiseFilterResultRun2Tight",&hBHENoiseFilterTightBit,"HBHENoiseFilterResultRun2Tight/I");
	skimTreeOutput->Branch("phfCoincFilter", &hfCoincidenceFilterBit, "phfCoincFilter/I");
	skimTreeOutput->Branch("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, "pVertexFilterCutdz1p0/I");
	skimTreeOutput->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplusBit,"pVertexFilterCutGplus/I");
	skimTreeOutput->Branch("pVertexFilterCutVtx1",&pVertexFilterCutVtx1Bit,"pVertexFilterCutVtx1/I");


 	// Copy the jet trees to the output
	TTree *jetTreeOutput[nJetTrees];
	
	// Leaves for jet tree
	Int_t nJetsOutput[nJetTrees];										// number of jets in an event
	Float_t jetPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};				// phi of all the jets in an event
	Float_t jetPhiArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};			// phi of all the jets in an event	with WTA axis
	Float_t jetEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};				// eta of all the jets in an event
	Float_t jetEtaArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};			// eta of all the jets in an event	with WTA axis
	Float_t jetRawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// raw jet pT for all the jets in an event
	Float_t jetMaxTrackPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// maximum track pT inside a jet for all the jets in an event

	Float_t jetRefPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level pT for a reconstructed jet
	Float_t jetRefEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level eta for a reconstructed jet
	Float_t jetRefPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level phi for a reconstructed jet
	Int_t jetRefFlavorArrayOutput[nJetTrees][nMaxJet] = {{0}};			// flavor for initiating parton for the reference gen jet
	Int_t jetRefFlavorForBArrayOutput[nJetTrees][nMaxJet] = {{0}};		// heavy flavor for initiating parton for the reference gen jet
	Int_t jetRefSubidArrayOutput[nJetTrees][nMaxJet] = {{0}};           // jet subid


	Int_t nGenJetsOutput[nJetTrees];								 	// number of generator level jets in an event
	Float_t genJetPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// pT of all the generator level jets in an event
	Float_t genJetPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};			// phi of all the generator level jets in an event
	Float_t genJetPhiArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the generator level jets in an event with WTA axis
	Float_t genJetEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};			// eta of all the generator level jets in an event
	Float_t genJetEtaArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the generator level jets in an event with WTA axis
	Int_t genJetSubidArrayOutput[nJetTrees][nMaxJet] = {{0}};     		// subid of all the generator level jets in an event
	Int_t genJetMatchIndexArrayOutput[nJetTrees][nMaxJet] = {{0}};		// matched index of all the generator level jets in an event


	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
		jetTreeOutput[iJetType] = new TTree("t","");
		
		jetTreeOutput[iJetType]->Branch("nref",&nJetsOutput[iJetType],"nref/I");
		jetTreeOutput[iJetType]->Branch("rawpt",&jetRawPtArrayOutput[iJetType],"rawpt[nref]/F");
		jetTreeOutput[iJetType]->Branch("trackMax",&jetMaxTrackPtArrayOutput[iJetType],"trackMax[nref]/F");

		// Jet eta with E-scheme and WTA axes
		jetTreeOutput[iJetType]->Branch("jtphi",&jetPhiArrayOutput[iJetType],"jtphi[nref]/F");
		jetTreeOutput[iJetType]->Branch("WTAphi",&jetPhiArrayWTAOutput[iJetType],"WTAphi[nref]/F");
		
		// Jet phi with E-scheme and WTA axes
		jetTreeOutput[iJetType]->Branch("jteta",&jetEtaArrayOutput[iJetType],"jteta[nref]/F");
		jetTreeOutput[iJetType]->Branch("WTAeta",&jetEtaArrayWTAOutput[iJetType],"WTAeta[nref]/F");
		
		// If we are looking at Monte Carlo, connect the reference pT and parton arrays
		if(is_MC){
			jetTreeOutput[iJetType]->Branch("refpt",&jetRefPtArrayOutput[iJetType],"refpt[nref]/F");
			jetTreeOutput[iJetType]->Branch("refeta",&jetRefEtaArrayOutput[iJetType],"refeta[nref]/F");
			jetTreeOutput[iJetType]->Branch("refphi",&jetRefPhiArrayOutput[iJetType],"refphi[nref]/F");
			jetTreeOutput[iJetType]->Branch("refparton_flavor", &jetRefFlavorArrayOutput[iJetType], "refparton_flavor[nref]/I");
			jetTreeOutput[iJetType]->Branch("refparton_flavorForB", &jetRefFlavorForBArrayOutput[iJetType], "refparton_flavorForB[nref]/I");
			jetTreeOutput[iJetType]->Branch("subid", &jetRefSubidArrayOutput[iJetType], "subid[nref]/I");
		 
			jetTreeOutput[iJetType]->Branch("ngen",&nGenJetsOutput[iJetType],"ngen/I");
			jetTreeOutput[iJetType]->Branch("genpt",&genJetPtArrayOutput[iJetType],"genpt[ngen]/F");
			
			// Gen jet phi for e-scheme and WTA axes
			jetTreeOutput[iJetType]->Branch("genphi",&genJetPhiArrayOutput[iJetType],"genphi[ngen]/F");
			jetTreeOutput[iJetType]->Branch("WTAgenphi",&genJetPhiArrayWTAOutput[iJetType],"WTAgenphi[ngen]/F");
			
			// Gen jet eta for e-scheme and WTA axes
			jetTreeOutput[iJetType]->Branch("geneta",&genJetEtaArrayOutput[iJetType],"geneta[ngen]/F");
			jetTreeOutput[iJetType]->Branch("WTAgeneta",&genJetEtaArrayWTAOutput[iJetType],"WTAgeneta[ngen]/F");

			// Gen match and subid
			jetTreeOutput[iJetType]->Branch("genmatchindex",&genJetMatchIndexArrayOutput[iJetType],"genmatchindex[ngen]/F");
			jetTreeOutput[iJetType]->Branch("gensubid",&genJetSubidArrayOutput[iJetType],"gensubid[ngen]/F");
		
		} // Branches only for MC

	} // Jet type loop
	
	// Copy the track trees to the output
	TTree *trackTreeOutput = new TTree("trackTree","");
	Int_t nTracksOutput;										 // Number of tracks
	Float_t trackPtOutput[nMaxTrack] = {0};		 				 // Array for track pT:s
	Float_t trackPtErrorOutput[nMaxTrack] = {0};				 // Array for track pT errors
	Float_t trackPhiOutput[nMaxTrack] = {0};					 // Array for track phis
	Float_t trackEtaOutput[nMaxTrack] = {0};					 // Array for track etas
	Bool_t trackHighPurityOutput[nMaxTrack] = {0};				 // Array for the high purity of tracks
	Float_t trackVertexDistanceZOutput[nMaxTrack] = {0};		 // Array for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceZErrorOutput[nMaxTrack] = {0};	 // Array for error for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceXYOutput[nMaxTrack] = {0};		 // Array for track distance from primary vertex in xy-direction
	Float_t trackVertexDistanceXYErrorOutput[nMaxTrack] = {0};	 // Array for error for track distance from primary vertex in xy-direction
	Float_t trackEnergyEcalOutput[nMaxTrack] = {0};							// Array for track energy in ECal
	Float_t trackEnergyHcalOutput[nMaxTrack] = {0};							// Array for track energy in HCal
	UChar_t PixelnHitsTrackOutput[nMaxTrack] = {0};				 // Array for number of hits for the track
	Int_t trackChargeOutput[nMaxTrack] = {0};					 // Array for track charge
	
	trackTreeOutput->Branch("nTrk",&nTracksOutput,"nTrk/I");
	trackTreeOutput->Branch("trkPt",&trackPtOutput,"trkPt[nTrk]/F");
	trackTreeOutput->Branch("trkPtError",&trackPtErrorOutput,"trkPtError[nTrk]/F");
	trackTreeOutput->Branch("trkPhi",&trackPhiOutput,"trkPhi[nTrk]/F");
	trackTreeOutput->Branch("trkEta",&trackEtaOutput,"trkEta[nTrk]/F");
	trackTreeOutput->Branch("highPurity",&trackHighPurityOutput,"highPurity[nTrk]/O");
	trackTreeOutput->Branch("trkDz1",&trackVertexDistanceZOutput,"trkDz1[nTrk]/F");
	trackTreeOutput->Branch("trkDzError1",&trackVertexDistanceZErrorOutput,"trkDzError1[nTrk]/F");
	trackTreeOutput->Branch("trkDxy1",&trackVertexDistanceXYOutput,"trkDxy1[nTrk]/F");
	trackTreeOutput->Branch("trkDxyError1",&trackVertexDistanceXYErrorOutput,"trkDxyError1[nTrk]/F");
	trackTreeOutput->Branch("trkNPixelHit",&PixelnHitsTrackOutput,"trkNPixelHit[nTrk]/b");
	trackTreeOutput->Branch("pfEcal",&trackEnergyEcalOutput,"pfEcal[nTrk]/F");
	trackTreeOutput->Branch("pfHcal",&trackEnergyHcalOutput,"pfHcal[nTrk]/F");
	trackTreeOutput->Branch("trkCharge",&trackChargeOutput,"trkCharge[nTrk]/I");
	
	// Generator level tracks only in Monte Carlo
	TTree *genTrackTreeOutput = new TTree("hi","");
	std::vector<float> *genTrackPtVector = new std::vector<float>(); genTrackPtVector->clear();
	std::vector<float> *genTrackPhiVector = new std::vector<float>(); genTrackPhiVector->clear();
	std::vector<float> *genTrackEtaVector = new std::vector<float>(); genTrackEtaVector->clear();
	std::vector<int> *genTrackPdgVector = new std::vector<int>(); genTrackPdgVector->clear();
	std::vector<int> *genTrackChargeVector = new std::vector<int>(); genTrackChargeVector->clear();
	std::vector<int> *genTrackSubeventVector = new std::vector<int>(); genTrackSubeventVector->clear();
	// Connect the branches to generator level track tree
	if(is_MC){
		genTrackTreeOutput->Branch("pt","vector<float>", &genTrackPtVector);
		genTrackTreeOutput->Branch("phi","vector<float>", &genTrackPhiVector);
		genTrackTreeOutput->Branch("eta","vector<float>", &genTrackEtaVector);
		genTrackTreeOutput->Branch("pdg","vector<int>", &genTrackPdgVector);
		genTrackTreeOutput->Branch("chg","vector<int>", &genTrackChargeVector);
		genTrackTreeOutput->Branch("sube","vector<int>", &genTrackSubeventVector);
	}

	// ========================================== //
	//				Loop over all events 					//
	// ========================================== //
	
	int nEvents = heavyIonTree->GetEntries();
	cout << "There are " << nEvents << " events" << endl;
	
	bool passTrackCuts;
	bool passJetCuts;
	int iTrackOutput;
	int iJetOutput;

	for(int iEvent = 0; iEvent < nEvents; iEvent++) {
		
		if( iEvent % 1000 == 0 )	std::cout << "iEvent: " << iEvent <<	" of " << nEvents << std::endl;

		// ========================================== //
		//	Read the event to input trees	      //
		// ========================================== //
		
		heavyIonTree->GetEntry(iEvent);
		hltTree->GetEntry(iEvent);
		skimTree->GetEntry(iEvent);
		trackTree->GetEntry(iEvent);
		if(is_MC) genTrackTree->GetEntry(iEvent);
		particleFlowCandidateTree->GetEntry(iEvent);
		checkFlatteningTree->GetEntry(iEvent);

		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
			jetTree[iJetType]->GetEntry(iEvent);
		}

		int multiplicity = get_Ntrkoff(nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);

		bool multsel = true;
		if(ntrkoff==0){if(iEvent==0){cout << "No multiplicity cut" << endl;}}
		if(ntrkoff==1){if(iEvent==0){cout << "MB: [0,185]" << endl;} if(multiplicity >= 185){multsel=false;}}
		if(ntrkoff==2){if(iEvent==0){cout << "HM 1 to 6: [185,250]" << endl;} if(multiplicity < 185 || multiplicity >= 250){multsel=false;}}
		if(ntrkoff==3){if(iEvent==0){cout << "HM 7: [250,inf]" << endl;} if(multiplicity < 250){multsel=false;}}

		if(multsel==false) continue;		

		heavyIonTreeOutput->Fill();
		hltTreeOutput->Fill();
		skimTreeOutput->Fill();

		//Event plane (just what we want EP from 2 to 6)
		epang_HFm2 = (float) eventPlaneAngle[44];
		epang_HFp2 = (float) eventPlaneAngle[45];
		epang_HFm3 = (float) eventPlaneAngle[73];
		epang_HFp3 = (float) eventPlaneAngle[74];
		epang_HFm4 = (float) eventPlaneAngle[102];
		epang_HFp4 = (float) eventPlaneAngle[103];
		epang_HFm5 = (float) eventPlaneAngle[131];
		epang_HFp5 = (float) eventPlaneAngle[132];
		epang_HFm6 = (float) eventPlaneAngle[148];
		epang_HFp6 = (float) eventPlaneAngle[149];

		q_HFm2 = (float) eventPlaneQ[44];
		q_HFp2 = (float) eventPlaneQ[45];
		q_HFm3 = (float) eventPlaneQ[73];
		q_HFp3 = (float) eventPlaneQ[74];
		q_HFm4 = (float) eventPlaneQ[102];
		q_HFp4 = (float) eventPlaneQ[103];
		q_HFm5 = (float) eventPlaneQ[131];
		q_HFp5 = (float) eventPlaneQ[132];
		q_HFm6 = (float) eventPlaneQ[148];
		q_HFp6 = (float) eventPlaneQ[149];

		qx_HFm2 = (float) eventPlaneQx[44];
		qx_HFp2 = (float) eventPlaneQx[45];
		qx_HFm3 = (float) eventPlaneQx[73];
		qx_HFp3 = (float) eventPlaneQx[74];
		qx_HFm4 = (float) eventPlaneQx[102];
		qx_HFp4 = (float) eventPlaneQx[103];
		qx_HFm5 = (float) eventPlaneQx[131];
		qx_HFp5 = (float) eventPlaneQx[132];
		qx_HFm6 = (float) eventPlaneQx[148];
		qx_HFp6 = (float) eventPlaneQx[149];
		
		qy_HFm2 = (float) eventPlaneQy[44];
		qy_HFp2 = (float) eventPlaneQy[45];
		qy_HFm3 = (float) eventPlaneQy[73];
		qy_HFp3 = (float) eventPlaneQy[74];
		qy_HFm4 = (float) eventPlaneQy[102];
		qy_HFp4 = (float) eventPlaneQy[103];
		qy_HFm5 = (float) eventPlaneQy[131];
		qy_HFp5 = (float) eventPlaneQy[132];
		qy_HFm6 = (float) eventPlaneQy[148];
		qy_HFp6 = (float) eventPlaneQy[149];

		mult_HFm2 = (float) eventPlaneMultiplicity[44];
		mult_HFp2 = (float) eventPlaneMultiplicity[45];
		mult_HFm3 = (float) eventPlaneMultiplicity[73];
		mult_HFp3 = (float) eventPlaneMultiplicity[74];
		mult_HFm4 = (float) eventPlaneMultiplicity[102];
		mult_HFp4 = (float) eventPlaneMultiplicity[103];
		mult_HFm5 = (float) eventPlaneMultiplicity[131];
		mult_HFp5 = (float) eventPlaneMultiplicity[132];
		mult_HFm6 = (float) eventPlaneMultiplicity[148];
		mult_HFp6 = (float) eventPlaneMultiplicity[149];

		checkFlatteningTreeOutput->Fill();

    	// Fill jet histograms using basic jet cuts
		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
			iJetOutput = 0;
			nJetsOutput[iJetType] = nJets[iJetType];
		
			for(int iJet = 0; iJet < nJets[iJetType]; iJet++){ // loop over reco jets

				passJetCuts = true;

				vector<Node *> NodesWTAScheme;
				NodesWTAScheme.clear();

				double jetR = 0.4;
				if( iJetType == 3 ) jetR = 0.3;
				Float_t jetPhiWTA = -999;
				Float_t jetEtaWTA = -999;

				// recluster and find WTA axis
				if( iJetType == 0){
					for(int itrk = 0; itrk < nTracks; itrk++) {
     					// Set particle kinematics
     					double deltaphi = jetPhiArray[iJetType][iJet] - trackPhiArray[itrk];
     					while(deltaphi > (TMath::Pi())){deltaphi += -2*TMath::Pi();}
     					while(deltaphi < (-1.0*TMath::Pi())){deltaphi += 2*TMath::Pi();}
     					double deltaeta = jetEtaArray[iJetType][iJet] - trackEtaArray[itrk];
     					double deltaR = sqrt(pow(deltaphi,2) + pow(deltaeta,2));
     					if(deltaR > jetR) continue;
      					FourVector P;
      					P.SetPtEtaPhiMass(trackPtArray[itrk], trackEtaArray[itrk], trackPhiArray[itrk], 0);
      					// Add into the node object vector
      					NodesWTAScheme.push_back(new Node(P));
					}					
				}else{
					for(int pfi = 0; pfi < particleFlowCandidatePtVector->size(); pfi++) {
     					// Set particle kinematics
     					double deltaphi = jetPhiArray[iJetType][iJet] - particleFlowCandidatePhiVector->at(pfi);
     					while(deltaphi > (TMath::Pi())){deltaphi += -2*TMath::Pi();}
     					while(deltaphi < (-1.0*TMath::Pi())){deltaphi += 2*TMath::Pi();}
     					double deltaeta = jetEtaArray[iJetType][iJet] - particleFlowCandidateEtaVector->at(pfi);
     					double deltaR = sqrt(pow(deltaphi,2) + pow(deltaeta,2));
     					if(deltaR > jetR) continue;
      					FourVector P;
      					P.SetPtEtaPhiMass(particleFlowCandidatePtVector->at(pfi), particleFlowCandidateEtaVector->at(pfi), particleFlowCandidatePhiVector->at(pfi), 0);
      					// Add into the node object vector
      					NodesWTAScheme.push_back(new Node(P));
   					}
				}	

				// Do the reclustering!
   				if(NodesWTAScheme.size()>0){
   					BuildCATree(NodesWTAScheme, -1, WTAScheme);
   					//BuildCATree(NodesWTAScheme, -1, EScheme);
  					jetPhiWTA = NodesWTAScheme[0]->P.GetPhi();
  					jetEtaWTA = NodesWTAScheme[0]->P.GetEta();
  					delete NodesWTAScheme[0];
  					NodesWTAScheme.clear();
				}


				// Apply very basic jet cuts
				if(jetRawPtArray[iJetType][iJet] < jetptmin) passJetCuts = false;    // Minumum pT cut of 30 GeV
				if(fabs(jetEtaArray[iJetType][iJet]) > jetetamin && fabs(jetEtaWTA) > jetetamin) passJetCuts = false;    // Maximum eta cut of 2.1
			
				// Fill the jet arrays with reconstructed jets
				if(passJetCuts){
					jetRawPtArrayOutput[iJetType][iJetOutput] = jetRawPtArray[iJetType][iJet];
					jetPhiArrayOutput[iJetType][iJetOutput] = jetPhiArray[iJetType][iJet];
					jetPhiArrayWTAOutput[iJetType][iJetOutput] = jetPhiWTA;
					jetEtaArrayOutput[iJetType][iJetOutput] = jetEtaArray[iJetType][iJet];
					jetEtaArrayWTAOutput[iJetType][iJetOutput] = jetEtaWTA;
					jetMaxTrackPtArrayOutput[iJetType][iJetOutput] = jetMaxTrackPtArray[iJetType][iJet];
				
					if(is_MC){
						jetRefPtArrayOutput[iJetType][iJetOutput] = jetRefPtArray[iJetType][iJet];
						jetRefEtaArrayOutput[iJetType][iJetOutput] = jetRefEtaArray[iJetType][iJet];
						jetRefPhiArrayOutput[iJetType][iJetOutput] = jetRefPhiArray[iJetType][iJet];
						jetRefFlavorArrayOutput[iJetType][iJetOutput] = jetRefFlavorArray[iJetType][iJet];
						jetRefFlavorForBArrayOutput[iJetType][iJetOutput] = jetRefFlavorForBArray[iJetType][iJet];
						jetRefSubidArrayOutput[iJetType][iJetOutput] = jetRefSubidArray[iJetType][iJet];						
					}
					iJetOutput++;
				} else {nJetsOutput[iJetType]--;}
			} // Reconstructed jet loop
		
			if(is_MC){
			
				iJetOutput = 0;
				nGenJetsOutput[iJetType] = nGenJets[iJetType];
			
				for(int iJet = 0; iJet < nGenJets[iJetType]; iJet++){
				
					passJetCuts = true;

					vector<Node *> NodesWTASchemeGen;
					NodesWTASchemeGen.clear();

					double jetRGen = 0.4;
					if( iJetType == 3 ) jetRGen = 0.3;
					Float_t jetPhiWTAGen = -999;
					Float_t jetEtaWTAGen = -999;

					for(int gpi = 0; gpi < genTrackPtArray->size(); gpi++) {
     					// Set particle kinematics
     					double deltaphi = genJetPhiArray[iJetType][iJet] - genTrackPhiArray->at(gpi);
     					while(deltaphi > (TMath::Pi())){deltaphi += -2*TMath::Pi();}
     					while(deltaphi < (-1.0*TMath::Pi())){deltaphi += 2*TMath::Pi();}
     					double deltaeta = genJetEtaArray[iJetType][iJet] - genTrackEtaArray->at(gpi);
     					double deltaR = sqrt(pow(deltaphi,2) + pow(deltaeta,2));
     					if(deltaR > jetRGen) continue;
      					FourVector P;
      					P.SetPtEtaPhiMass(genTrackPtArray->at(gpi), genTrackEtaArray->at(gpi), genTrackPhiArray->at(gpi), 0);
      					// Add into the node object vector
      					NodesWTASchemeGen.push_back(new Node(P));
   					}

   					// Do the reclustering!
   					if(NodesWTASchemeGen.size()>0){
   						BuildCATree(NodesWTASchemeGen, -1, WTAScheme);
   						//BuildCATree(NodesWTASchemeGen, -1, EScheme);
  						jetPhiWTAGen = NodesWTASchemeGen[0]->P.GetPhi();
  						jetEtaWTAGen = NodesWTASchemeGen[0]->P.GetEta();
  						delete NodesWTASchemeGen[0];
  						NodesWTASchemeGen.clear();
					}

					// Apply very basic jet cuts
					if(genJetPtArray[iJetType][iJet] < jetptmin) passJetCuts = false;    // Minimum pT cut of 30 GeV
					if(fabs(genJetEtaArray[iJetType][iJet]) > jetetamin && fabs(jetEtaWTAGen) > jetetamin) passJetCuts = false;    // Maximum eta cut of 2.1
				
					// Fill the jet arrays with generated jets
					if(passJetCuts){
				
					genJetPtArrayOutput[iJetType][iJetOutput] = genJetPtArray[iJetType][iJet];
					genJetPhiArrayOutput[iJetType][iJetOutput] = genJetPhiArray[iJetType][iJet];
					genJetPhiArrayWTAOutput[iJetType][iJetOutput] = jetPhiWTAGen;
					genJetEtaArrayOutput[iJetType][iJetOutput] = genJetEtaArray[iJetType][iJet];
					genJetEtaArrayWTAOutput[iJetType][iJetOutput] = jetEtaWTAGen;
					genJetSubidArrayOutput[iJetType][iJetOutput] = genJetSubidArray[iJetType][iJet];
					genJetMatchIndexArrayOutput[iJetType][iJetOutput] = genJetMatchIndexArray[iJetType][iJet];

					iJetOutput++;

					} else {nGenJetsOutput[iJetType]--;}
					} // Generator level jet loop
			} // If for filling generator jet loop
			jetTreeOutput[iJetType]->Fill();
    	} // Loop over jet collections

	    // Reco track loop
    	nTracksOutput = nTracks;
    	iTrackOutput = 0;
    	for(int iTrack = 0; iTrack < nTracks; iTrack++){
      
    		passTrackCuts = true;
      
    		// Do basic track cuts for the reconstructed tracks
    		if(trackHighPurityArray[iTrack] != 1) passTrackCuts = false;
    		if(fabs(trackPtErrorArray[iTrack]/trackPtArray[iTrack]) > 0.15) passTrackCuts = false;
    		if(fabs(trackVertexDistanceZArray[iTrack]/trackVertexDistanceZErrorArray[iTrack]) > 5.0) passTrackCuts = false;
    		if(fabs(trackVertexDistanceXYArray[iTrack]/trackVertexDistanceXYErrorArray[iTrack]) > 5.0) passTrackCuts = false;
    		if(fabs(trackEtaArray[iTrack]) >= 2.4) passTrackCuts = false;  //acceptance of the tracker
    		if(trackPtArray[iTrack] <= 0.4) passTrackCuts = false;   // Minimum track pT
      
    		if(passTrackCuts){
    			trackPtOutput[iTrackOutput] = trackPtArray[iTrack];
    			trackPtErrorOutput[iTrackOutput] = trackPtErrorArray[iTrack];
    			trackPhiOutput[iTrackOutput] = trackPhiArray[iTrack];
        		trackEtaOutput[iTrackOutput] = trackEtaArray[iTrack];
        		trackHighPurityOutput[iTrackOutput] = trackHighPurityArray[iTrack];
        		trackVertexDistanceZOutput[iTrackOutput] = trackVertexDistanceZArray[iTrack];
        		trackVertexDistanceZErrorOutput[iTrackOutput] = trackVertexDistanceZErrorArray[iTrack];
        		trackVertexDistanceXYOutput[iTrackOutput] = trackVertexDistanceXYArray[iTrack];
        		trackVertexDistanceXYErrorOutput[iTrackOutput] = trackVertexDistanceXYErrorArray[iTrack];
        		trackEnergyEcalOutput[iTrackOutput] = trackEnergyEcalArray[iTrack];
        		trackEnergyHcalOutput[iTrackOutput] = trackEnergyHcalArray[iTrack];
        		trackChargeOutput[iTrackOutput] = trackChargeArray[iTrack];
        		PixelnHitsTrackOutput[iTrackOutput] = PixelnHitsTrackArray[iTrack];
        		iTrackOutput++;
      		} else {nTracksOutput--;} 
    	}
    	trackTreeOutput->Fill();


    	// Gen track loop
    	if(is_MC){
    		for(int iTrack = 0; iTrack < genTrackPtArray->size(); iTrack++){
    			// Cut away low pT tracks and tracks with eta outside of tracker acceptance
			if(TMath::Abs(genTrackEtaArray->at(iTrack)) >= 2.4) continue; //acceptance of the tracker
			if(genTrackPtArray->at(iTrack) <= 0.4) continue;   // Minimum track pT
		        // Fill the output vectors with gen particles surviving the cuts
        		genTrackPtVector->push_back(genTrackPtArray->at(iTrack));
        		genTrackPhiVector->push_back(genTrackPhiArray->at(iTrack));
        		genTrackEtaVector->push_back(genTrackEtaArray->at(iTrack));
       			genTrackPdgVector->push_back(genTrackPdgArray->at(iTrack));
        		genTrackChargeVector->push_back(genTrackChargeArray->at(iTrack));
        		genTrackSubeventVector->push_back(genTrackSubeventArray->at(iTrack));
      		}
      
      		genTrackTreeOutput->Fill();
    		// Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
       		genTrackPtVector->clear();
      		genTrackPhiVector->clear();
      		genTrackEtaVector->clear();
      		genTrackPdgVector->clear();
      		genTrackChargeVector->clear();
      		genTrackSubeventVector->clear();

    	} // Filling gen tracks for MC

   		// Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
    	//particleFlowCandidatePtOutputVector->clear();
    	//particleFlowCandidatePhiOutputVector->clear();
    	//particleFlowCandidateEtaOutputVector->clear();

	} // End loop over events

	// Write the skimmed trees to the output file
  	TFile *outputFile = new TFile(outputFileName, "RECREATE");

	gDirectory->mkdir("hiEvtAnalyzer");
	gDirectory->cd("hiEvtAnalyzer");
	heavyIonTreeOutput->Write();

	gDirectory->cd("../");
	gDirectory->mkdir("hltanalysis");
	gDirectory->cd("hltanalysis");
	hltTreeOutput->Write();
  
	gDirectory->cd("../");
	gDirectory->mkdir("skimanalysis");
	gDirectory->cd("skimanalysis");
	skimTreeOutput->Write();

	gDirectory->cd("../");
	gDirectory->mkdir("checkflattening");
	gDirectory->cd("checkflattening");
	checkFlatteningTreeOutput->Write();

	gDirectory->cd("../");
	const char *jetDirectories[] = {"ak4CaloJetAnalyzer","ak4PFJetAnalyzer","akCs4PFJetAnalyzer","ak3PFJetAnalyzer"};
	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		gDirectory->mkdir(jetDirectories[iJetType]);
		gDirectory->cd(jetDirectories[iJetType]);
		jetTreeOutput[iJetType]->Write();
		gDirectory->cd("../");
	} // Loop over jet types
  
  	gDirectory->mkdir("ppTrack");
  	gDirectory->cd("ppTrack");
	trackTreeOutput->Write();
	gDirectory->cd("../");

	// Generator particles only present in MC
	if(is_MC){
		gDirectory->mkdir("HiGenParticleAna");
		gDirectory->cd("HiGenParticleAna");
		genTrackTreeOutput->Write();
		gDirectory->cd("../");
	}
  
	outputFile->Close();

	cout << endl;
	cout << "------------------------------------- SKIMMER DONE --------------------------------------" << endl;
	cout << endl;


	sec_end = clock(); // stop time counting
	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

	print_stop(); // Print time, date and hour when it stops

}

int main(int argc, char** argv){
				TString firstArgument(argv[1]);
				TString outfile(argv[2]);
				int mc = atoi(argv[3]);
				int ntrkoffline = atoi(argv[4]);
				pPbSkim(firstArgument,outfile,mc,ntrkoffline);
}
