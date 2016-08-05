#include "TH1D.h"
#include "selector.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


Double_t selector::BPRES_mips(TVector3 v_photon, Double_t conerad)
{

  TVector3 v_bpres_mip;

  Double_t bpres_mips_e = 0.;
  for(Int_t i=0; i<Nbpchn; i++)
    {
      v_bpres_mip.SetXYZ(Bpxyz[i][0]-Xvtx, Bpxyz[i][1]-Yvtx, Bpxyz[i][2]-Zvtx);
      Double_t dr = v_photon.DeltaR(v_bpres_mip);
      if(check_cuts) 
        if (nodebugmode) cout << "bpre: "<< i << " "<<  Bpmip[i] << " " << bpres_mips_e <<" " << dr << endl;
      if(dr < conerad) {
  //  cout << "Bpmip[" << i << "] = " << Bpmip[i] << endl; 
  bpres_mips_e += Bpmip[i];
      }
    }

  //  if(!Data) {
    //    cout << "emip before dobpre: " << bpres_mips_e << endl; 
  //    bpres_mips_e = doBpreReadout(10);//bpres_mips_e);
    //    cout << "emip after dobpre: " << bpres_mips_e << endl; 
  //  }
  return bpres_mips_e;
}

void  selector::getBpreReadout() // A.Volynets
{
// read landay readout
   string name = "input/landin.dat";
   ifstream myin1(name.c_str(), ios::in); //only check at the beginning
   if (myin1)
   {
//   read landay
     if (nodebugmode) cout << "Global::getBpreReadout(): read " << "input/landin.dat" << endl;
     for (Int_t i=0; i<100; i++)
     {
       myin1 >> Land[i];
     // G4cout << " Land " << i << " " << Land[i] << G4endl;
     }
   } else
   {
     if (nodebugmode) cout << "Global::getBpreReadout(): Cannot open READOUT file " << endl;
   }
   if (nodebugmode) cout << "Land:" << endl;
   for(Int_t i=0; i<100; i++)
     if (nodebugmode) cout << i << " " << Land[i] << endl;
   myin1.close();
}

Double_t selector::doBpreReadout(Double_t mips) // A.Volynets
{
  // start smearing
   Double_t array[100];
   Int_t   irans[100];

   gRandom->RndmArray(100, array);
   for (Int_t i=0; i<100; i++)
   {
     //   G4cout << array[i] << G4endl;
     irans[i]=(int)(array[i]*100)+1;
     if (irans[i] < 1)   irans[i]=1;
     if (irans[i] > 100) irans[i]=100;
   }

   Double_t vv=0;
   Double_t vcc[100];
   if (mips>100)
   {
     for (Int_t i = 0; i < 100; i++)
     {
       vcc[i] = Land[ irans[i] ];
       vv = vv + vcc[i];
     }
   } else
   if (mips>1)
   {
     for (Int_t i = 0; i < (int)mips; i++)
     {
       vcc[i] = Land[ irans[i] ];
       vv = vv + vcc[i];
     }
   } else  vv=mips;
   return  vv;
}

void selector::getRunsBPRE() // A.Volynets
{
    string name="input/badruns2004_2007.d";
    ifstream myfile;
    myfile.open(name.c_str(), ios::in);
    Int_t irun=0; Int_t ibad=0;

    if (!myfile)
    {
      cerr << "Global::getRunsBPRE():  Can't open input file:  " << name << endl;
      exit(1);
    } else
    {
       cout << "-> Read file=" << name << endl;
    }

    while ( myfile >> irun >> ibad )
    {
       // cout << irun << "  " << ibad << endl;
        BadBpre.insert ( make_pair ( irun, ibad ) );
    }
    myfile.close();
}
