#include "TH1D.h"
#include "TGraph.h"
#include "selector.h"
#include "constants.h"
#include <iostream>
#include <string>
using namespace std;


Double_t selector::stretch_calib(Double_t value,       // Value of the distribution
				 string method,       // "e5" for ELEC5 or "zu" for Zufo
				 string distribution, // Fmax or deltaZ
				 Double_t eta,
				 Double_t et)
{
   // If not prompt photon MC - return initial value
  if (Data || (!Data && mc_type == "mc_bg_norad")) return value;

//  cout<<"Debug:: Calibration is in process.."<<endl;
  // Construct graph name
  string graph_name;
  /*******************
    Method recognition:
  ********************/
       if (method=="e5") graph_name+="elec5_";
  else if (method=="zu") graph_name+="zufo_";
  else {cout<<"Wrong method. Return -1"<<endl; return -1;}

  /*******************
    Distrubution
  ********************/
       if (distribution=="fmax") graph_name+="fmax_";
  else if (distribution=="deltaz") graph_name+="dz_";
  else {cout<<"Wrong distribution. Return -1"<<endl; return -1;}

  /*******************
    Pseudorapidity range
  ********************/
  graph_name+= "calib_eta_";
       if ((eta>=-0.7)&&(eta<-0.3))  graph_name+="-0.7_-0.3_";
  else if ((eta>=-0.3)&&(eta< 0.1))  graph_name+="-0.3_0.1_";
  else if ((eta>= 0.1)&&(eta< 0.5))  graph_name+="0.1_0.5_";
  else if ((eta>= 0.5)&&(eta< 0.9))  graph_name+="0.5_0.9_";
  else {/*cout<<"Wrong pseudorapidity. Return initial value"<<endl;*/ return value;}

  /*******************
    Transverse energy
  ********************/
  graph_name+= "Et_";
       if ((et>=4. )&&(et<6. )) graph_name+=" 4_ 6";
  else if ((et>=6. )&&(et<8. )) graph_name+=" 6_ 8";
  else if ((et>=8. )&&(et<10.)) graph_name+=" 8_10";
  else if ((et>=10.)&&(et<15.)) graph_name+="10_15";
  else {/*cout<<"Wrong transverse energy. Return initial value"<<endl;*/ return value;}

  // Opening calibration data file
  TGraph *gCalib = (TGraph*)(fCalib)->Get(graph_name.c_str());
  Double_t return_value = gCalib->Eval(value);
  if (!gCalib)
  delete gCalib;

/*
  Double_t x,y;
  for (int iii=0;iii<200;iii++)
  {
    gCalib->GetPoint(iii,x,y);
    cout<<"Point "<<iii<<"  "<<x<<"  "<<y<<endl;
  }
*/
/*
  Int_t argc = 0;
  char **argv = 0;
  TApplication theApp("App", &argc, argv);
  TBrowser *b = 0;
  b = new TBrowser("RunMC Browser", fCalib, "Histograms", 10, 20, 600 , 300);
//  gCalib->Draw();
  theApp.Run();
*/
  return return_value;
}

