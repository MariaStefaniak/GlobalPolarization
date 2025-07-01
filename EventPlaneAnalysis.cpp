#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2.h>
#include <iostream>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <vector>
#include <string>
#include <TRegexp.h>
#include <TMath.h>

class Event {
    public:
        ULong64_t       outGPSTIME;
        ULong64_t       outEVENTNUMBER;
        Float_t         outPVX;   
        Float_t         outPVY;  
        Float_t         outPVZ;  
        UInt_t          outRUNNUMBER;
        Int_t           outnBackTracks;
        Int_t           outnVeloClusters;
        Int_t           outnVeloTracks;
        Int_t           outnEcalClusters;
        Double_t        outQx_back[2];
        Double_t        outQy_back[2];
        Double_t        outQx_for[2][4];
        Double_t        outQy_for[2][4];
        Int_t           out_Qmulti[4];
    
        Event() : outGPSTIME(0), outEVENTNUMBER(0), outRUNNUMBER(0){}
};
/*
std::vector<std::string> GetRootFilesInDirectory(const std::string& dirPath) {
    std::vector<std::string> fileNames;

    TSystemDirectory dir("dir", dirPath.c_str());
    TList* files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Could not open directory: " << dirPath << std::endl;
        return fileNames;
    }

    TIter next(files);
    TSystemFile* file;
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
            fileNames.push_back(dirPath + "/" + fname.Data());
        }
    }

    delete files;
    return fileNames;
}*/

std::vector<std::string> GetFilteredRootFiles(const std::string& dirPath) {
    std::vector<std::string> fileNames;

    TSystemDirectory dir("dir", dirPath.c_str());
    TList* files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Could not open or read directory: " << dirPath << std::endl;
        return fileNames;
    }

    TIter next(files);
    TSystemFile* file;

    // Match: 00274156_00000###_1.tuple_pbpb2024.root
    TRegexp pattern("^00274156_00000[0-4][0-9][0-9]_1\\.tuple_pbpb2024\\.root$");

    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(pattern)) {
            std::string fullpath = dirPath + "/" + fname.Data();
            std::cout << "Adding file: " << fullpath << std::endl;  // DEBUG
            fileNames.push_back(fullpath);
        }
    }

    delete files;
    return fileNames;
}

void EventPlaneAnalysis(){

    const int maxTracks = 10000; 
    std::string eosDir = "/eos/lhcb/grid/prod/lhcb/anaprod/lhcb/LHCb/Lead24/TUPLE_PBPB2024.ROOT/00274156/0000";
    std::vector<std::string> fileNames = GetFilteredRootFiles(eosDir);  
  //  std::string eosDir = "/eos/lhcb/grid/prod/lhcb/anaprod/lhcb/LHCb/Lead24/TUPLE_PBPB2024.ROOT/00274156/0000";
  //  std::vector<std::string> fileNames = GetRootFilesInDirectory(eosDir);
   /* std::vector<std::string> fileNames = {
       // "/Volumes/Mike_disc/Maria/PbPb/VELO/velo1/pbpb_velo_1.root",
       // "/Volumes/Mike_disc/Maria/PbPb/VELO/velo1/pbpb_velo_2.root"
       "/eos/lhcb/grid/prod/lhcb/anaprod/lhcb/LHCb/Lead24/TUPLE_PBPB2024.ROOT/00274156/0000/00274156_00000020_1.tuple_pbpb2024.root"
        // Add more files here
    };*/
    // Create output file and tree
    TFile* outFile = new TFile("June30/event_plane_pbpb_0_499.root", "RECREATE");
    TTree* outTree = new TTree("EventPlaneTuple", "Event Plane");

    // Create one Event object and connect it to the tree
    Event* evt = new Event();
    outTree->Branch("event", &evt);

    for (const auto& fileName : fileNames) { // loop over files:

        TFile *file = TFile::Open(fileName.c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Could not open file: " << fileName << std::endl;
            continue;
        }

        // Navigate to the subdirectory
        TDirectory* dir = (TDirectory*)file->Get("EventTuplePV");
        if (!dir) {
            std::cerr << "Directory 'EventTuplePV' not found in file: " << fileName << std::endl;
            file->Close();
            continue;
        }

        // Get the tree from the directory
        TTree* tree = (TTree*)dir->Get("EventTuplePV"); // replace with your actual tree name
        if (!tree) {
            std::cerr << "Tree not found in directory 'EventTuplePV' in file: " << fileName << std::endl;
            file->Close();
            continue;
        }
        std::cerr << "file: " << fileName << " opened" << std::endl;

        ULong64_t       GPSTIME;
        ULong64_t       EVENTNUMBER;
        Float_t         PVX[maxTracks];   
        Float_t         PVY[maxTracks];   
        Float_t         PVZ[maxTracks];   
        UInt_t          RUNNUMBER;
        Float_t         VELOTRACK_BIPCHI2[maxTracks];   
        Float_t         VELOTRACK_ETA[maxTracks];  
        Float_t         VELOTRACK_ISBACKWARD[maxTracks];  
        Float_t         VELOTRACK_NVPHITS[maxTracks];   
        Float_t         VELOTRACK_PHI[maxTracks];   
        Float_t         VELOTRACK_PX[maxTracks];   
        Float_t         VELOTRACK_PY[maxTracks];  
        Float_t         VELOTRACK_PZ[maxTracks];  
        Float_t         VELOTRACK_RHO[maxTracks];   
        Int_t           nBackTracks;
        Int_t           nPVs;
        Int_t           nVeloClusters;
        Int_t           nVeloTracks;
        Int_t           nEcalClusters;

        //Setting branches:
        tree->SetBranchAddress("GPSTIME", &GPSTIME);
        tree->SetBranchAddress("EVENTNUMBER", &EVENTNUMBER);
        tree->SetBranchAddress("PVX", &PVX);
        tree->SetBranchAddress("PVY", &PVY);
        tree->SetBranchAddress("PVZ", &PVZ);
        tree->SetBranchAddress("RUNNUMBER", &RUNNUMBER);
        tree->SetBranchAddress("VELOTRACK_BIPCHI2", &VELOTRACK_BIPCHI2);
        tree->SetBranchAddress("VELOTRACK_ETA", &VELOTRACK_ETA);
        tree->SetBranchAddress("VELOTRACK_ISBACKWARD", &VELOTRACK_ISBACKWARD);
        tree->SetBranchAddress("VELOTRACK_NVPHITS", &VELOTRACK_NVPHITS);
        tree->SetBranchAddress("VELOTRACK_PHI", &VELOTRACK_PHI);
        tree->SetBranchAddress("VELOTRACK_PX", &VELOTRACK_PX);
        tree->SetBranchAddress("VELOTRACK_PY", &VELOTRACK_PY);
        tree->SetBranchAddress("VELOTRACK_PZ", &VELOTRACK_PZ);
        tree->SetBranchAddress("VELOTRACK_RHO", &VELOTRACK_RHO);
        tree->SetBranchAddress("nBackTracks", &nBackTracks);
        tree->SetBranchAddress("nPVs", &nPVs);
        tree->SetBranchAddress("nVeloClusters", &nVeloClusters);
        tree->SetBranchAddress("nVeloTracks", &nVeloTracks);
        tree->SetBranchAddress("nEcalClusters", &nEcalClusters);


        // Loop over entries
        Long64_t nEvents = tree->GetEntries();
        for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {
            tree->GetEntry(iEvent);

            // EVENT CUTS:
            // - only one collsion per event:
            if(nPVs!=1)continue;
            // SMOG or Collider (should be all PbPb but...)
            if(nBackTracks < 10) continue;
            if(PVZ[0]<-100 || PVZ[0]>100) continue;
            if(nVeloTracks<15) continue;


            int nPrimaryForTracks = 0; int nPrimaryBackTracks = 0; 
            int nPrimaryForTracks_1 = 0; int nPrimaryForTracks_2 = 0; int nPrimaryForTracks_3 = 0;

                 //Event plane Q vectors:
            double Qx_back[2], Qy_back[2], Qx_for[2][4], Qy_for[2][4];
            for(int i = 0; i < 2; i++){
                Qx_back[i] = 0; Qy_back[i] = 0;
                for(int j = 0; j < 4; j++){
                    Qx_for[i][j] = 0; Qy_for[i][j] = 0;
                }
            }

            for(int iTrack = 0; iTrack < nVeloTracks; iTrack++){
                
                // ============  quality checks ================
                if(VELOTRACK_BIPCHI2[iTrack]>1.5)continue;
                if(VELOTRACK_ISBACKWARD[iTrack]==1) {
                    VELOTRACK_ETA[iTrack] = VELOTRACK_ETA[iTrack]*(-1);
                    VELOTRACK_PHI[iTrack] = VELOTRACK_PHI[iTrack]+TMath::Pi();
                } /// if it is backward switch eta -> -eta, and phi -> phi+pi
                if(VELOTRACK_ISBACKWARD[iTrack]!=1) nPrimaryForTracks++;
                if(VELOTRACK_ISBACKWARD[iTrack]==1) nPrimaryBackTracks++;

                // Q - vectors:

                double w1=1;//VELOTRACK_ETA[iTrack];      
                double w2=1;            
                //backward Q:
                if(VELOTRACK_ETA[iTrack] < -0.5){
                    Qx_back[0] =  Qx_back[0] + (w1 * cos(1 * VELOTRACK_PHI[iTrack]));
                    Qx_back[1] =  Qx_back[1] + (w2 * cos(2 * VELOTRACK_PHI[iTrack]));

                }
                // forwards Q eta 1:
                if(VELOTRACK_ETA[iTrack] > 0.5 && VELOTRACK_ETA[iTrack] <= 2.5 ){ // 1-2.5 -4-6
                    Qx_for[0][0] =  Qx_for[0][0] + (w1 * cos(1 * VELOTRACK_PHI[iTrack]));
                    Qx_for[1][0] =  Qx_for[1][0] + (w2 * cos(2 * VELOTRACK_PHI[iTrack]));

                    nPrimaryForTracks_1++;
                }
                // forwards Q eta 2:
                if(VELOTRACK_ETA[iTrack] > 2.5 && VELOTRACK_ETA[iTrack] <= 4.0 ){
                    Qx_for[0][1] =  Qx_for[0][1] + (w1 * cos(1 * VELOTRACK_PHI[iTrack]));
                    Qx_for[1][1] =  Qx_for[1][1] + (w2 * cos(2 * VELOTRACK_PHI[iTrack]));
                    
                    nPrimaryForTracks_2++;
                }
                // forwards Q eta 3:
                if(VELOTRACK_ETA[iTrack] > 4.0 && VELOTRACK_ETA[iTrack] <= 6.0 ){
                    Qx_for[0][2] =  Qx_for[0][2] + (w1 * cos(1 * VELOTRACK_PHI[iTrack]));
                    Qx_for[1][2] =  Qx_for[1][2] + (w2 * cos(2 * VELOTRACK_PHI[iTrack]));

                    nPrimaryForTracks_3++;
                }
                // forwards Q eta 3:
                if(VELOTRACK_ETA[iTrack] > 0.5 && VELOTRACK_ETA[iTrack] <= 6.0 ){
                    Qx_for[0][3] =  Qx_for[0][3] + (w1 * cos(1 * VELOTRACK_PHI[iTrack]));
                    Qx_for[1][3] =  Qx_for[1][3] + (w2 * cos(2 * VELOTRACK_PHI[iTrack]));
                 
                    nPrimaryForTracks++; 
                }

                
            }// end of iTrack loop
            
            // check if we have enough tracks to get both EP
            if(nPrimaryForTracks<5)     continue;
            if(nPrimaryBackTracks<5)    continue;
            if(nPrimaryForTracks_1<5)   continue;
            if(nPrimaryForTracks_2<5)   continue;
            if(nPrimaryForTracks_3<5)   continue;

        
            evt->outGPSTIME         = GPSTIME;
            evt->outEVENTNUMBER     = EVENTNUMBER;
            evt->outPVX             = PVX[0];
            evt->outPVY             = PVY[0]; 
            evt->outPVZ             = PVZ[0];  
            evt->outRUNNUMBER       = RUNNUMBER;
            evt->outnBackTracks     = nBackTracks;
            evt->outnVeloClusters   = nVeloClusters;
            evt->outnVeloTracks     = nVeloTracks;
            evt->outnEcalClusters   = nEcalClusters;
            evt->out_Qmulti[0]    = nPrimaryForTracks_1;
            evt->out_Qmulti[1]    = nPrimaryForTracks_2;
            evt->out_Qmulti[2]    = nPrimaryForTracks_3;
            evt->out_Qmulti[3]    = nPrimaryBackTracks;

            
            for(int i = 0; i < 2; i++){
                evt->outQx_back[i] = Qx_back[i]; 
                evt->outQy_back[i] = Qy_back[i]; 
                for(int j = 0; j < 4; j++){
                    evt->outQx_for[i][j] = Qx_for[i][j]; 
                    evt->outQy_for[i][j] = Qy_for[i][j];
                }
            }
            outTree->Fill(); 
        }// end of iEvent loop
        file->Close();
        delete file;
    } // end of file loop

   // Save and close
    outFile->cd();
    outTree->Write();
    outFile->Close();
    std::cout << "Done. Saved to event_plane_pbpb.root" << std::endl;
    delete evt; 

}