#include <stdlib.h>
#include <iostream>
using namespace std;

#include <TChain.h>
#include "04pNtuplesNames.h"
#include "0405eBGntuplesNames.h" //warning
#include "0405e_2BGntuplesNames.h" //warning
//#include "0405e_3BGntuplesNames.h" //warning  THIS MAKES NO SENCE - INSIDE DIFFERENT RUN
//#include "0405ePHPDIR_PYTHIAntuplesNames.h"//warning
#include "0405eNtuplesNames.h"
#include "0607pBGntuplesNames.h"
#include "0607pNtuplesNames.h"
#include "06eBGntuplesNames.h"
#include "06e_2BGntuplesNames.h"
#include "06eNtuplesNames.h"
#include "0405ePRPHntuplesNames.h"//warning
//#include "0405ePRPHPHPntuplesNames.h"//warning
#include "06ePRPHntuplesNames.h"
#include "PrivateNtuplesNames.h"
void fill(TChain* chain, TString location, TString* str, Int_t n);

void fill_chain(TChain* chain, TString run_period, Bool_t data, TString mc_type, Bool_t use_clustered)
{
  	if(mc_type == "mc_prph")
	{
		if(run_period == "0405e") 
		{
			   //not working yet cout << "filenum2: " << filenum << " should: " << argv[6] << endl;
			// Ian MC - new Funnel
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSBF25.H25125.8B.PRPH_DIS_PYT64_2005_01.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_02.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_03.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_04.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_05.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_06.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_07.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_08.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_09.root");
				chain->Add("/nfs/dust/zeus/group/glushenko/new_zmcsp/samples/Ian/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_10.root");
			// Ian MC - old Funnel
				//without EM-conserv: num_had = 69727 num_part = 62370
				//with EM-conserv: num_had = 39647 num_part = 46993
				
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSBF25.H25125.8B.PRPH_DIS_PYT64_2005_01_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_02_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_03_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_04_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_05_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_06_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_07_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_08_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_09_1.root");
				// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/Ian/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_10_1.root");
				
			
			// My MC
				//Old
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/25000_ev/ZEUSMC.SDSBF25.H25125.8B.PRPH_DIS_PYT64_2005_01_1.root");                        
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/25000_ev/ZEUSMC.SDSBF25.H25125.8B.PRPH_DIS_PYT64_2005_01_2.root");			
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_02_1.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_02_2.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_03_1.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_03_2.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_04_1.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_04_2.root");
					// //-//
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_05_1.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_05_2.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_06_1.root");
					// //problem chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_06_2.root");
					// //-//
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_07_1.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_07_2.root");
					// //-//
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_08_1.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_08_2.root");
					// //-//
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_09_1.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_09_2.root");
					// //-//
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_10_1.root");
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/250000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_10_2.root");
				//----//
				//New mv mc_prph0405e_parton.root mc_prph0405e_parton_new_8.root 
				// hadd mc_prph0405e_parton_new_new_with_EM_conserv.root  mc_prph0405e_parton_new_1.root    mc_prph0405e_parton_new_4.root  mc_prph0405e_parton_new_6-7.root  mc_prph0405e_parton_new_2_1.root  mc_prph0405e_parton_new_5.root  mc_prph0405e_parton_new_8.root
					//no EM check num_had = 4467 num_part = 4004
					//with EM check num_had = 2307 num_part = 2720
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSBF25.H25125.8B.PRPH_DIS_PYT64_2005_01_1.root");// 1		
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_02_1.root");// 2
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_04_1.root");// 4
					//  chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_05_1.root");// 5
					//  chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_06_1.root");// 6
					//  chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_07_1.root");// 6
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_08_1.root");// 7
					//chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSNE25.H25125.8B.PRPH_DIS_PYT64_2005_03_1.root");//no
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_09_1.root");// NO
					// chain->Add("/nfs/dust/zeus/group/glushenko/zmcsp/samples/100000_ev/ZEUSMC.SDSME25.H25125.8B.PRPH_DIS_PYT64_2005_10_1.root");// NO
				
			
			cout << chain->GetEntries() << " events chained" << endl;

			//Old MC
			/*
			if (!use_clustered)
			  fill(chain, str_location_prph0405e, rootFile_prph0405e, numberOfRootFiles_prph0405e);
			else
			  fill(chain, str_location_prph0405eprivate, rootFile_prph0405eprivate, numberOfRootFiles_prph0405eprivate);
			  */      
		} 
	    else if(run_period == "06e") 
	    {
			if (!use_clustered)
			  fill(chain, str_location_prph06e, rootFile_prph06e, numberOfRootFiles_prph06e);
			else
			  fill(chain, str_location_prph06eprivate, rootFile_prph06eprivate, numberOfRootFiles_prph06eprivate);      
	    } 
	    else 
	    {
	    	cerr << "mc_prph: unknown period: " << run_period << endl; 
	    	exit(-1);
	    }   
	}
  	else
    {
		if(run_period == "0405e")
		{
			if(data)
			{
			  if(!use_clustered)
				fill(chain, str_location0405e, rootFile0405e, numberOfRootFiles0405e);
			  else
				fill(chain, str_location0405eprivate, rootFile0405eprivate, numberOfRootFiles0405eprivate);
			} 
			else
			{
				if(mc_type == "mc_bg_rad" || mc_type == "mc_bg_norad")
				{
					if(!use_clustered)
					{
						fill(chain, str_location0405e_bg, rootFile0405e_bg, numberOfRootFiles0405e_bg);
						fill(chain, str_location0405e_2bg, rootFile0405e_2bg, numberOfRootFiles0405e_2bg);//warning
						//fill(chain, str_location0405e_3bg, rootFile0405e_3bg, numberOfRootFiles0405e_3bg);//warning
					}
					else
						fill(chain, str_location0405eleptoprivate, rootFile0405eleptoprivate, numberOfRootFiles0405eleptoprivate);		
				} 
				else
				{
					cout << "unknown mc type: " << mc_type << endl;
					exit(-1);
				}
			}
		} 
		else if(run_period == "06e")
		{
			if(data)
			{
				if(!use_clustered)
				fill(chain, str_location06e, rootFile06e, numberOfRootFiles06e);
				else
				fill(chain, str_location06eprivate, rootFile06eprivate, numberOfRootFiles06eprivate);
			} 
			else
			{
				if(mc_type == "mc_bg_rad" || mc_type == "mc_bg_norad")
				{
					if(!use_clustered) 
					{
						fill(chain, str_location06e_bg, rootFile06e_bg, numberOfRootFiles06e_bg);
						fill(chain, str_location06e_2bg, rootFile06e_2bg, numberOfRootFiles06e_2bg);
					}
					else
						fill(chain, str_location06eleptoprivate, rootFile06eleptoprivate, numberOfRootFiles06eleptoprivate);
				} 
				else
				{
					cout << "unknown mc type: " << mc_type << endl;
					exit(-1);
				}
			}
		}
      	else if(run_period == "0607p")
		{
			if(data)
			{
				if(!use_clustered)
					fill(chain, str_location0607p, rootFile0607p, numberOfRootFiles0607p);
				else
					fill(chain, str_location0607pprivate, rootFile0607pprivate, numberOfRootFiles0607pprivate);
			} 
			else
			{
				if(mc_type == "mc_bg_rad" || mc_type == "mc_bg_norad")
				{
					if(!use_clustered)
						fill(chain, str_location0607p_bg, rootFile0607p_bg, numberOfRootFiles0607p_bg);
					else
						fill(chain, str_location0607pleptoprivate, rootFile0607pleptoprivate, numberOfRootFiles0607pleptoprivate);
				} 
				else
				{
					cout << "unknown mc type: " << mc_type << endl;
					exit(-1);
				}
			}
		}
      	else if(run_period == "04p")
		{
			if(data)
			{
				if(!use_clustered)
					fill(chain, str_location04p, rootFile04p, numberOfRootFiles04p);
				else
					fill(chain, str_location04pprivate, rootFile04pprivate, numberOfRootFiles04pprivate);
			} 
		}
	    else if(run_period == "all")
		{
			if(data) 
			{
				if(!use_clustered) 
				{
					fill(chain, str_location04p, rootFile04p, numberOfRootFiles04p);
					fill(chain, str_location0405e, rootFile0405e, numberOfRootFiles0405e);
					fill(chain, str_location06e, rootFile06e, numberOfRootFiles06e);
					fill(chain, str_location0607p, rootFile0607p, numberOfRootFiles0607p);
				}
				else 
				{
					fill(chain, str_location04pprivate, rootFile04pprivate, numberOfRootFiles04pprivate);
					fill(chain, str_location0405eprivate, rootFile0405eprivate, numberOfRootFiles0405eprivate);
					fill(chain, str_location06eprivate, rootFile06eprivate, numberOfRootFiles06eprivate);
					fill(chain, str_location0607pprivate, rootFile0607pprivate, numberOfRootFiles0607pprivate);
				}
			}
			else
			{
				if(mc_type == "mc_bg_rad" || mc_type == "mc_bg_norad")
				{
					if(!use_clustered) 
					{
						fill(chain, str_location0405e_bg, rootFile0405e_bg, numberOfRootFiles0405e_bg);
						fill(chain, str_location06e_bg, rootFile06e_bg, numberOfRootFiles06e_bg);
						fill(chain, str_location0607p_bg, rootFile0607p_bg, numberOfRootFiles0607p_bg);
					}
					else
						fill(chain, str_location0607pleptoprivate, rootFile0607pleptoprivate, numberOfRootFiles0607pleptoprivate);
				} 
				else
				{
					cout << "unknown mc type: " << mc_type << endl;
					exit(-1);
				}
			}
		}
   		else
		{ 
			cout << "error. Unknown period: " << run_period << endl; 
			exit(-1);
		} 
	}  //  chain->Add("selectedEvents_lepto_989900.root");
  //chain->Add("1.root");
	cout << chain->GetEntries() << " events chained" << endl;
}


void fill(TChain* chain, TString location, TString* str, Int_t n)
{
  for(Int_t i = 0; i < n; i++)
    
    {
      TString file = location + str[i];
      chain->Add(file);
      cout << file << " added to chain in fill" << endl;
    }
    cout << chain->GetEntries() << " events chained" << endl;
}
