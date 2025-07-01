void CalculateResolution(){
    TFile *inFile = new TFile("Psi12_PbPbVELO_nrCent3.root","READ");
    TH2D *hResolutionIn = (TH2D*)inFile->Get("hResolutionCentrality");
    
    int nrCent = 3;
    double Resolution1[nrCent];
    double Resolution2[nrCent];

    for (int i = 0; i < nrCent; i++){
        double r1= hResolutionIn->GetBinContent(1,i+1);
        double r2= hResolutionIn->GetBinContent(2,i+1);
        double nr = hResolutionIn->GetBinContent(3,i+1);
        cout << r2 << "  " << nr << endl;
        Resolution1[i] = sqrt(2*r1/nr) *100;
        Resolution2[i] = sqrt(2*r2/nr) *100;
    }

    int centBins[4] = {14  , 126 ,   270 ,   1000};
    double x_bins[3] = {centBins[1] + (centBins[1]-centBins[0])*0.5,centBins[2] + (centBins[2]-centBins[1])*0.5,  centBins[3] + (centBins[3]-centBins[2])*0.5,};
   
    TGraph *grResolution1 = new TGraph(3,x_bins, Resolution1);
    TGraph *grResolution2 = new TGraph(3,x_bins, Resolution2);
 
    grResolution1->SetMarkerStyle(8);
    grResolution2->SetMarkerStyle(8);
    grResolution2->GetXaxis()->SetTitle("nrVeloTracks");
    grResolution2->GetYaxis()->SetTitle("Resolution #Psi_{2}[%]");
    grResolution2->GetYaxis()->SetRangeUser(0, 100);
    grResolution1->GetXaxis()->SetTitle("nrVeloTracks");
    grResolution1->GetYaxis()->SetTitle("Resolution #Psi_{1}[%]");
    grResolution1->GetYaxis()->SetRangeUser(0, 100);

    TCanvas *c7 = new TCanvas();
    c7->Divide(2,1);
    c7->cd(1);
    grResolution1->Draw("ap");
    c7->cd(2);
    grResolution2->Draw("ap");

    cout << "x: " << x_bins[0]      << " , " << x_bins[1]      << ",  " << x_bins[2]      << endl;
    cout << "R1 " << Resolution1[0] << " , " << Resolution1[1] << " , " << Resolution1[2] << endl;
    cout << "R2 " << Resolution2[0] << " , " << Resolution2[1] << " , " << Resolution2[2] << endl;
 

}