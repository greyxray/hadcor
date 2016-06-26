#include <iostream>
#include <fstream>
#include <vector>
#include "selector.h"

using namespace std;

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
//#include "Dijets.h"
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"
using KtJet::KtLorentzVector;
using KtJet::KtEvent;

Bool_t selector::SelectHadronLevel(Bool_t take_det_event)
{
  //  cout << "...going into hadr. selection" << endl;
  Bool_t take_hevent = kTRUE;
  Bool_t here_is_true_prph = kTRUE;
  Bool_t is_true_prph_candidate = kTRUE;
  Bool_t here_is_true_jet = kFALSE;
  Bool_t take_true_jet = kTRUE;
  Double_t accomp_jet_et = 0.;
  Double_t accomp_jet_eta = -999.;
  Double_t accomp_jet_phi = -999.;
  hadron_level_jet_cont_photon_index = -1;
  hadron_level_ephoton_over_ejet = -1;
 
  //event selection
  //  if(Eventnr == 7085 || Eventnr == 8575 || Eventnr == 13318 || Eventnr == 10212 || Eventnr == 13437)
  //    check_cuts = kTRUE;

  if(check_cuts) cout << "going into hadron level Eventnr = " << Eventnr << endl;
  hist.mccorel_q2_y->Fill(Mc_q2, Mc_y, wtx);
  hist.mccorel_q2_x->Fill(Mc_q2, Mc_x, wtx);
  hist.mccorel_q2_y_noweight->Fill(Mc_q2, Mc_y);
  hist.mccorel_q2_x_noweight->Fill(Mc_q2, Mc_x);
  if(Mc_q2 < 10.) {
    take_hevent = kFALSE;
    if(check_cuts)
      cout << "rejected by cut on true q2 = " << Mc_q2 << endl;
  }
  TVector3 v_true_electron(Mc_pfsl[0],Mc_pfsl[1],Mc_pfsl[2]);
  v_true_scattered_electron->SetPxPyPzE(Mc_pfsl[0],Mc_pfsl[1],Mc_pfsl[2],Mc_pfsl[3]);
  if ((v_true_electron.Theta()*180.0/TMath::Pi() < 139.8 )||
      (v_true_electron.Theta()*180.0/TMath::Pi() > 180.0 ))  {
    //      (v_true_electron.Theta()*180.0/TMath::Pi() > 171.9 ))  { //warning
    take_hevent = kFALSE;
    if(check_cuts)
      cout << "rejected by cut on true electron theta = " << v_true_electron.Theta()*180.0/TMath::Pi() << endl;
  }
  if(Mc_pfsl[3] < 10.)
    {
      take_hevent = kFALSE;
      if(check_cuts)
	cout << "rejected by cut on true electron energy = " << Mc_pfsl[3] << endl;
    }

  if(take_hevent) {
    hist.dis_Q2_true->Fill(Mc_q2, wtx);
    hist.dis_electron_e_true->Fill(Mc_pfsl[3], wtx);
    hist.dis_electron_theta_true->Fill(v_true_electron.Theta()*180.0/TMath::Pi(), wtx);
    hist.dis_x_true->Fill(Mc_x, wtx);
    hist.dis_y_true->Fill(Mc_y, wtx);
  }
  //find true photon
  Int_t index_true_photon = -1;
  Int_t index_jet = -1;//this variable is to compare with Sanja, Nazar, Ian, Natasha
  if(!Data && (mc_type == "mc_bg_rad")/* || mc_type == "mc_bg_norad")*/)
    {
      for (Int_t i=0; i<Npart; i++)
	{
	  if (Part_id[i]==5)
	    {
	      index_jet=Part_jetid[i]-1;
	      index_true_photon=i;
	      break;
	    }
	}
    }
  if(!Data && mc_type == "mc_prph")
    {
      for (Int_t i=0; i<Npart; i++)
	{
	  if (Part_id[i]==12)
	    {
	      index_jet=Part_jetid[i]-1;
	      index_true_photon=i;
	      break;
	    }
	}
    }
  if(index_true_photon < 0)
    {
      here_is_true_prph = kFALSE;
      if(check_cuts)
	cout << "there is no photon in generated level" << endl;
    }
  is_true_prph_candidate = here_is_true_prph;
  if(take_hevent && here_is_true_prph) {
    TVector3 v_true_photon(Part_p[index_true_photon][0],Part_p[index_true_photon][1],Part_p[index_true_photon][2]);
    v_true_prompt_photon->SetPxPyPzE(Part_p[index_true_photon][0],Part_p[index_true_photon][1],Part_p[index_true_photon][2], Part_p[index_true_photon][3]);
    v_true_jet_cont_prompt_photon->SetPxPyPzE(Ktrjets[index_jet][0], Ktrjets[index_jet][1], Ktrjets[index_jet][2], Ktrjets[index_jet][3]);
    Double_t true_photon_et  = Part_p[index_true_photon][3]*TMath::Sin(v_true_photon.Theta());
    Double_t true_photon_eta = v_true_photon.Eta();
    if (true_photon_et<4. || true_photon_et>15.) {
      is_true_prph_candidate = kFALSE;
      if(check_cuts)
	cout << "rejected by tue_photon_et = " << true_photon_et << endl;
    }
    if (true_photon_eta<-0.7 || true_photon_eta>0.9) {
      is_true_prph_candidate = kFALSE;
      if(check_cuts)
	cout << "rejected by cut on true photon eta = " << true_photon_eta << endl;
    }

    //    if (Ktrnjets<1) {
    //      is_true_prph_candidate = kFALSE;
    //      if(check_cuts)
    //	cout << "rejected by Ktrnjets<1: " << Ktrnjets << endl;
    //    }
    //    if (index_jet<0) {
    //      is_true_prph_candidate = kFALSE;
    //      if(check_cuts)
    //	cout << "rejected by index_jet<0: " << index_jet << endl;
    //    }
    //    if ((Part_p[index_true_photon][3]/Ktrjets[index_jet][3])<0.9) {
    //      is_true_prph_candidate = kFALSE;
    //      if(check_cuts)
    //	cout << "rejected by photon/jet energy ratio: " << Part_p[index_true_photon][3]/Ktrjets[index_jet][3] << endl;
    //    }
  } //if(take_hevent && here_is_true_prph)

  //  if(det_cross_sec_bin_et>=0)
  //    cout << "hadron level: here_is_true_prph = " << here_is_true_prph << endl;
  if(take_hevent && is_true_prph_candidate) {
  vector<KtJet::KtLorentzVector> input_hadrons;
  vector<KtJet::KtLorentzVector> true_jets;
  Int_t index_true_electron = -1;
  Int_t index_photon_vector = -1;
  Int_t index_photon_jet = -1;
  Double_t ephoton_over_ejet = -1;
  for(Int_t i=0; i<Npart; i++)
    {
      if ((Part_prt[i]==23)||(Part_prt[i]==24))
	{
	  index_true_electron = i;
	  break;
	}
    }
  KtJet::KtLorentzVector r;
  for (Int_t i=0; i<Npart; i++)
    {
      // Skip DIS electron
      Double_t M2 = Part_p[i][3]*Part_p[i][3] - Part_p[i][2]*Part_p[i][2]-Part_p[i][1]*Part_p[i][1]-Part_p[i][0]*Part_p[i][0];
	
      if (i==index_true_electron) continue;
      if (TMath::Abs(Part_p[i][3]) < 1.e-3) continue;
      if(Part_p[i][3] <= 0) continue;
      //	if(M2 < 0 ) continue;
      //      if(TMath::Abs(Part_p[i][2]) > Part_p[i][3]) continue;
      if (i==index_true_photon) index_photon_vector = input_hadrons.size();


      // create a KtJet with TRUE particles and put it onto
      // back of the input_particles vector
      if (TMath::Abs(Part_p[i][2])<Part_p[i][3])
	r = KtJet::KtLorentzVector(Part_p[i][0],Part_p[i][1],Part_p[i][2],Part_p[i][3]);
      else
	r = KtJet::KtLorentzVector(Part_p[i][0],Part_p[i][1],Part_p[i][2], 
				   TMath::Sqrt(Part_p[i][0]*Part_p[i][0] + Part_p[i][1]*Part_p[i][1] + Part_p[i][2]*Part_p[i][2] ));
	 
      

      input_hadrons.push_back(r);
    }
  double rparameter = 1.0;
  //  cout << "find hadr. level jets: " << endl;
  KtJet::KtEvent ev(input_hadrons,3,2,1,rparameter);
  //  cout << "OK" << endl;
  //
  // Get vector of jets
  //
  //  cout << "get jets..." << endl;
  true_jets = ev.getJetsEt();
  hadron_level_jets = true_jets;
  //  cout << "... OK" << endl;

  if((mc_type=="mc_prph" || mc_type=="mc_prph_rad" ) && index_true_photon<0)
    {
      cout<<"Warning in ReClusterize: index of TRUE photon in FMCKIN table < 0."<<endl;
      //	exit(-1);
      //      return jetsEtTRUE.size();
    }
  for (UInt_t i=0; i<true_jets.size(); i++)
    {
      if (true_jets[i].contains(input_hadrons[index_photon_vector]))
	{
	  index_photon_jet = i;
	  hadron_level_jet_cont_photon_index = i;
	  ephoton_over_ejet = input_hadrons[index_photon_vector].e() / true_jets[i].e();
	  hadron_level_ephoton_over_ejet = ephoton_over_ejet;
	  //	       << ", eta_photon = " << input_hadrons[index_photon_vector].eta() << endl;
	  break;
	}
    }
  //cout << "contains - ok" << endl;
  if (index_photon_jet<0)  { 
    is_true_prph_candidate = kFALSE;
    if(check_cuts)
      cout<<"Index of photon in jet < 0 - bad!!"<<endl; 
  }
  if(ephoton_over_ejet < 0.9) { 
    is_true_prph_candidate = kFALSE;
    if(check_cuts)
      cout << "rejected, because ephoton_over_ejet at hadron level is " << ephoton_over_ejet << endl;
  }
  here_is_true_prph = is_true_prph_candidate;
  if(take_hevent && here_is_true_prph) {    
    Double_t max_et_jet = -999.;
    Int_t index_of_accomp_jet = -1;

    for (UInt_t j=0; j<true_jets.size(); j++)
      {
	take_true_jet = kTRUE;
	// Calculate whether the Zu photon is in the current jet
	if (j==index_photon_jet) {
	  take_true_jet = kFALSE;
	  if(check_cuts) cout << "jet rejected, because it is a prompt photon jet" << endl;
	}
	// phase space
	if (true_jets[j].et() < ET_JET_CUT) {
	  take_true_jet = kFALSE;
	  if(check_cuts) cout << "jet rejected, because of its et: " << true_jets[j].et()<< endl;
	}  

	if (true_jets[j].eta() < -1.5 || true_jets[j].eta() > 1.8) {
	  take_true_jet = kFALSE;
	  if(check_cuts) cout << "jet rejected, because of its eta: " << true_jets[j].eta()<< endl;
	}
	
	if(take_true_jet)
	  {
	    here_is_true_jet = kTRUE;
	    count_nat++;

	    if(true_jets[j].et() > max_et_jet) {
	      max_et_jet = true_jets[j].et();
	      index_of_accomp_jet = j;
	      v_true_acc_jet->SetPxPyPzE(true_jets[j].px(), true_jets[j].py(), true_jets[j].pz(), true_jets[j].e());
	    }
	  }
      }//for j over true jets
    //Fill histograms

    if(here_is_true_jet && !inclusive_prph_analysis) {
      accomp_jet_et = true_jets[index_of_accomp_jet].et();
      accomp_jet_eta = true_jets[index_of_accomp_jet].eta();
      accomp_jet_phi = true_jets[index_of_accomp_jet].phi();
    }
    if(inclusive_prph_analysis)
      here_is_true_jet = kTRUE;
    if(take_hevent && here_is_true_prph && here_is_true_jet)
      {
	hist.had_prph_e->Fill(input_hadrons[index_photon_vector].et(), wtx);
	hist.had_jet_e->Fill(accomp_jet_et, wtx);
	hist.had_Q2->Fill(Mc_q2, wtx);
	hist.had_x->Fill(Mc_x, wtx);

	hist.had_cross_et->Fill(input_hadrons[index_photon_vector].et(), wtx);
	hist.had_cross_eta->Fill(input_hadrons[index_photon_vector].eta(), wtx);
	hist.had_cross_Q2->Fill(Mc_q2, wtx);
	hist.had_cross_x->Fill(Mc_x, wtx);
	hist.had_cross_et_jet->Fill(accomp_jet_et, wtx);
	hist.had_cross_et_jet2->Fill(accomp_jet_et, wtx);
	hist.had_cross_eta_jet->Fill(accomp_jet_eta, wtx);

	had_et = input_hadrons[index_photon_vector].et();
	had_eta = input_hadrons[index_photon_vector].eta();
	had_et_jet = accomp_jet_et;
	had_eta_jet = accomp_jet_eta;
	had_x = Mc_x;
	had_Q2 = Mc_q2;

	hist.prof_had_cross_et->Fill(input_hadrons[index_photon_vector].et(), input_hadrons[index_photon_vector].et(), wtx);
	hist.prof_had_cross_eta->Fill(input_hadrons[index_photon_vector].eta(), input_hadrons[index_photon_vector].eta(), wtx);
	hist.prof_had_cross_Q2->Fill(Mc_q2, Mc_q2, wtx);
	hist.prof_had_cross_x->Fill(Mc_x, Mc_x, wtx);
	hist.prof_had_cross_et_jet->Fill(accomp_jet_et, accomp_jet_et, wtx);
	hist.prof_had_cross_eta_jet->Fill(accomp_jet_eta, accomp_jet_eta, wtx);

	if(take_det_event) {
	  if(hist.had_cross_et->FindBin(input_hadrons[index_photon_vector].et()) == det_cross_sec_bin_et) {
	    hist.hd_cross_et->Fill(input_hadrons[index_photon_vector].et(), wtx);
	    hist.prof_hd_cross_et->Fill(input_hadrons[index_photon_vector].et(), input_hadrons[index_photon_vector].et(), wtx);
	  }
	  if(hist.had_cross_eta->FindBin(input_hadrons[index_photon_vector].eta()) == det_cross_sec_bin_eta) {
	    hist.hd_cross_eta->Fill(input_hadrons[index_photon_vector].eta(), wtx);
	    hist.prof_hd_cross_eta->Fill(input_hadrons[index_photon_vector].eta(), input_hadrons[index_photon_vector].eta(), wtx);
	  }
	  if(hist.had_cross_Q2->FindBin(Mc_q2) == det_cross_sec_bin_q2) {
	    hist.hd_cross_Q2->Fill(Mc_q2, wtx);
	    hist.prof_hd_cross_Q2->Fill(Mc_q2, Mc_q2, wtx);
	  }
	  if(hist.had_cross_x->FindBin(Mc_x) == det_cross_sec_bin_x) {
	    hist.hd_cross_x->Fill(Mc_x, wtx);
	    hist.prof_hd_cross_x->Fill(Mc_x, Mc_x, wtx);
	  }
	  if(hist.had_cross_et_jet->FindBin(accomp_jet_et) == det_cross_sec_bin_et_jet) {
	    hist.hd_cross_et_jet->Fill(accomp_jet_et, wtx);
	    hist.prof_hd_cross_et_jet->Fill(accomp_jet_et, accomp_jet_et, wtx);
	  }
	  if(hist.had_cross_eta_jet->FindBin(accomp_jet_eta) == det_cross_sec_bin_eta_jet) {
	    hist.hd_cross_eta_jet->Fill(accomp_jet_eta, wtx);
	    hist.prof_hd_cross_eta_jet->Fill(accomp_jet_eta, accomp_jet_eta, wtx);
	  }
	}// if take_det_event
      } //    if(take_hevent && here_is_true_prph && here_is_true_jet)    
  }//if(take_hevent && here_is_true_prph) to find jets and cut on prph kinematics

  //cout << "exiting hadron selection " << endl;
  //  check_cuts = kFALSE;
  if((take_hevent && here_is_true_prph && here_is_true_jet)) {
    cout << "hadr level, Eventnr = " << Eventnr << ", N of hadrons = " << Npart << ", N of jets = " << true_jets.size() << ": " << endl;
    for(Int_t i=0; i<true_jets.size(); i++) {
      cout << "ktjet #" << i << ": e = " << true_jets[i].e() << ", et = " << true_jets[i].et() << ", eta = " << true_jets[i].eta() << endl;;
    }
    cout << "found hadron level jet containing true photon: et = " << true_jets[index_photon_jet].et()
	 << ", eta = " << true_jets[index_photon_jet].eta() << " et_photon = " << input_hadrons[index_photon_vector].et()
	 << "e_photon = " << input_hadrons[index_photon_vector].e()
	 << "pz_photon = " << input_hadrons[index_photon_vector].pz() << endl;
    cout << "hadron level Eventnr = " << Eventnr << ": q2 = " << Mc_q2 << ", et_photon = " << input_hadrons[index_photon_vector].et()
	 << ", eta_photon = " << input_hadrons[index_photon_vector].eta()
	 << ", et_accomp_jet = " << accomp_jet_et 
	 << ", eta_accomp_jet = " << accomp_jet_eta 
	 << ", ephoton_over_ejet = " << ephoton_over_ejet << endl;
    cout << "================================HADRON LEVEL EVENT " << Eventnr << " ACCEPTED============================================" << endl;
  }
  }// if(take_hevent && here_is_true_prph), in order to time economy
  //  check_cuts = kFALSE;
  //  cout << "...going out hadr. selection" << endl;
  return (take_hevent && here_is_true_prph && here_is_true_jet);
}
