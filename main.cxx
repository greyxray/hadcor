//HADRONISATION CORRECTION CODE
#include <iostream>
#include <fstream>
#include <vector>


#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atoi */

using namespace std;

bool nodebugmode = false;
int filenum = 1;
#include <TH1D.h>
#include <TChain.h>
#include <TString.h>
#include "constants.h"
#include "selector.h"
#include "runinfo.h"
extern void fill_chain(TChain* chain, TString period, Bool_t data, TString mc_type, Bool_t use_clustered);

int main(int argc, char *argv[])
{
  bool nodebugmode = false;
  for (Int_t i = 0; i < argc; i++) 
    cout << "argv[" << i << "] = " << argv[i] << endl;
  Bool_t use_corr = kTRUE;
  TString s_use_corr = "1";
  Bool_t use_2ndcorr = kFALSE;
  TString s_use_2ndcorr = "1";
  Bool_t data = kTRUE;//kFALSE;
  TString s_data = "data";
  Bool_t use_clustered = kFALSE;
  TString s_use_clustered = "1";
  TString mc_type = "mc_bg_rad";
  TString mc_corr_type = "lepto_corr";
  if((argc == 2) && ((TString)argv[1] == "-h" || (TString)argv[1] == "--help") ) 
    {
    cout << "using: main <period> <data> <use_corr> <use_2ndcorr> <use_clustered>" << endl;
    cout << "  <period>: 0405e, 06e, 0607p" << endl;
    cout << "  <data>: data, mc_prph, mc_bg_rad, mc_bg_norad" << endl;
    //    cout << "  <mc_type_corr>: lepto_corr, ariadne_corr" << endl;
    cout << "  <use_corr>: 1, 0" << endl;
    cout << "  <use_2ndcorr>: 1, 0" << endl;
    cout << "  <use_clustered>: 1, 0" << endl;

    exit(0);
    }
  
  TString period = "989900";
  if (argc>1) period = (TString)argv[1];
  if (argc>2) 
  {
    s_data = (TString)argv[2];
    if (s_data == "data") data = kTRUE;
    else if (s_data == "mc_prph")
    {
      data = kFALSE;
      mc_type = "mc_prph";
    }
    else if (s_data == "mc_bg_rad")
    {
      data = kFALSE;
      mc_type = "mc_bg_rad";
    }
    else if (s_data == "mc_bg_norad")
		{
		  data = kFALSE;
		  mc_type = "mc_bg_norad";
		}
    else    
		{
		  cout << "unknown data type: " << s_data << endl;
		  exit(-1);
		}
  }
  if (argc>3) 
  {
    s_use_corr = (TString)argv[3];
    if (s_use_corr=="1") use_corr = kTRUE;
    else if(s_use_corr=="0") use_corr = kFALSE;
    else
  	{
  	  cout << "unknown parameter for use_corr (1 or 0): " << s_use_corr << endl;
  	  exit(-1);
  	}
  }
  if (argc>4) 
  {
    s_use_2ndcorr = (TString)argv[4];
    if (s_use_2ndcorr=="1") use_2ndcorr = kTRUE;
    else if (s_use_2ndcorr=="0") use_2ndcorr = kFALSE;
    else
    {
      cout << "unknown parameter for use_2ndcorr (1 or 0): " << s_use_2ndcorr << endl;
      exit(-1);
    }
  }
  if(argc>5) 
  {
    s_use_clustered = (TString)argv[5];
    if(s_use_clustered=="1") {
      use_clustered = kTRUE;
    } 
    else if(s_use_clustered=="0") {
	   use_clustered = kFALSE;
      } 
    else
  	{
  	  cout << "unknown parameter for use_clustered (1 or 0): " << s_use_clustered << endl;
  	  exit(-1);
  	}
  }
  if (argc > 6) 
  {
    filenum = atoi (argv[6]);
    cout << "filenum: " << filenum << " should: " << argv[6] << endl;
  }

  TString s_chain;
  s_chain = "orange";
  TChain* ch = new TChain(s_chain);
  fill_chain(ch, period, data, mc_type, use_clustered);
  //  ch->Add("0.root");
  //  ch->Add("1.root");
  //  ch->Add("62107_1.root");
  selector PromptPhotonPlusJetDIS;
  PromptPhotonPlusJetDIS.Init(ch, period, data, mc_type, mc_corr_type, use_corr, use_2ndcorr, use_clustered);
  if (argc>7)
  {
    PromptPhotonPlusJetDIS.q2_sufix =(TString)argv[7];//Form

    if  (PromptPhotonPlusJetDIS.q2_sufix.Contains("q2_"))  PromptPhotonPlusJetDIS.q2_sufix  = "_" + PromptPhotonPlusJetDIS.q2_sufix;

    if (PromptPhotonPlusJetDIS.q2_sufix.Contains("q2_lt_30"))       PromptPhotonPlusJetDIS.q2_cut_high  = 30;
    else if  (PromptPhotonPlusJetDIS.q2_sufix.Contains("q2_gt_30")) PromptPhotonPlusJetDIS.q2_cut_low  = 30;
  }
  PromptPhotonPlusJetDIS.Begin();
  //  for(Int_t i=0; i<ch->GetEntries(); i++) 
  PromptPhotonPlusJetDIS.Process();
  //  PromptPhotonPlusJetDIS.Terminate();
  return 0;
}
