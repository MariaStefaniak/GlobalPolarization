void ExecuteEPcalculations(){
    for(int ii = 1; ii <=1; ii++){
        gROOT->ProcessLine(".L calculateEventPlane.cpp+");
        gROOT->ProcessLine(Form("calculateEventPlane(%d)",ii));
    }
}

