void ExecuteEPcalculations(){
    for(int ii = 1; ii <=3; ii++){
        gROOT->ProcessLine(".L calculateEventPlane.cpp+");
        gROOT->ProcessLine(Form("calculateEventPlane(%d)",ii));
    }
}

