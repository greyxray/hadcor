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

Bool_t selector::SelectPartonLevel(Bool_t take_det_event, Bool_t take_had_event)
{
  //  cout << "going into part. selection" << endl;
  if(take_had_event) check_cuts = kTRUE;
  //
  // general histograms, filled for every event
  //===========================================
  //
  // event momentum at hadron level:
  Double_t had_px_sum = 892.48*TMath::Sin(2.7e-4);
  Double_t had_py_sum = 892.48*TMath::Sin(1.3e-4);
  Double_t had_pz_sum = -892.48;
  Double_t had_e_sum = - 27.52 - TMath::Sqrt(920*920+0.9382*0.9382);
  Double_t had_et_sum = 0;
  //  cout << "in SelectPartons: number of hadrons is " << Npart << endl;
  for(Int_t i=0; i<Npart; i++)
    {
      had_px_sum += Part_p[i][0];
      had_py_sum += Part_p[i][1];
      had_pz_sum += Part_p[i][2];
      had_e_sum += Part_p[i][3];
      TLorentzVector v; v.SetPxPyPzE(Part_p[i][0], Part_p[i][1], Part_p[i][2], Part_p[i][3]);
      had_et_sum += v.Et();
    }
  //  hist.had_sum_e->Fill(had_e_sum, wtx);
  //  hist.had_sum_px->Fill(had_px_sum, wtx);
  //  hist.had_sum_py->Fill(had_py_sum, wtx);
  //  hist.had_sum_pz->Fill(had_pz_sum, wtx);
  //  hist.had_sum_et->Fill(had_et_sum, wtx);
  //  cout << "event #" << Eventnr << ": " << had_e_sum << " " << had_px_sum << " " << had_py_sum << " " << had_pz_sum << " " << had_et_sum << endl;
  //
  // event momentum at parton level:
  Double_t part_px_sum = 892.48*TMath::Sin(2.7e-4);
  Double_t part_py_sum = 892.48*TMath::Sin(1.3e-4);
  Double_t part_pz_sum = -892.48;
  Double_t part_e_sum = - 27.52 - TMath::Sqrt(920*920+0.9382*0.9382);
  Double_t part_et_sum = 0;
  //  cout << "number of partons is " << nppart << endl;
  for(Int_t i=0; i<nppart; i++)
    {
      //      hist.part_idpart->Fill(idpart[i], wtx);
      part_px_sum += ppart[i][0];
      part_py_sum += ppart[i][1];
      part_pz_sum += ppart[i][2];
      part_e_sum += ppart[i][3];
      TLorentzVector v; v.SetPxPyPzE(ppart[i][0], ppart[i][1], ppart[i][2], ppart[i][3]);
      part_et_sum += v.Et();
    }
  //  hist.part_sum_e->Fill(part_e_sum, wtx);
  //  hist.part_sum_px->Fill(part_px_sum, wtx);
  //  hist.part_sum_py->Fill(part_py_sum, wtx);
  //  hist.part_sum_pz->Fill(part_pz_sum, wtx);
  //  hist.part_sum_et->Fill(part_et_sum, wtx);
  //  hist.part_nppart->Fill(nppart, wtx);
  Double_t m_q2_true = -1.*(bosene*bosene - bospx*bospx - bospy*bospy - bospz*bospz);
  Double_t m_x_true = m_q2_true / (2.*m_Proton_Energy*(bosene-bospz));
  Double_t m_y_true = m_q2_true / (4. *m_Lepton_Energy * m_Proton_Energy * m_x_true);
  hist.part_q2->Fill(m_q2_true, wtx);
  hist.part_y->Fill(m_y_true, wtx);
  hist.h2d_etevent_part_had->Fill(part_et_sum, had_et_sum, wtx);
  hist.had_q2->Fill(Mc_q2, wtx);
  hist.had_y->Fill(Mc_y, wtx);

  Bool_t take_pevent = kTRUE;
  Bool_t here_is_true_prph = kTRUE;
  Bool_t is_true_prph_candidate = kTRUE;
  Bool_t here_is_true_jet = kFALSE;
  Bool_t take_true_jet = kTRUE;
  Double_t ephoton_over_ejet = -1.;
  Double_t accomp_jet_et = -1;
  Double_t accomp_jet_eta = -1;
  if(Mc_q2 < 10.) {
    take_pevent = kFALSE;
    if(check_cuts)
      cout << "rejected by cut on true q2 = " << Mc_q2 << endl;
  }
  TVector3 v_true_electron(Mc_pfsl[0],Mc_pfsl[1],Mc_pfsl[2]);
  if ((v_true_electron.Theta()*180.0/TMath::Pi() < 139.8 )||
      (v_true_electron.Theta()*180.0/TMath::Pi() > 180.0 ))  {
    //      (v_true_electron.Theta()*180.0/TMath::Pi() > 171.9 ))  { //warning
    take_pevent = kFALSE;
    if(check_cuts)
      cout << "rejected by cut on true electron theta = " << v_true_electron.Theta()*180.0/TMath::Pi() << endl;
  }
  if(Mc_pfsl[3] < 10.)
    {
      take_pevent = kFALSE;
      if(check_cuts)
	cout << "rejected by cut on true electron energy = " << Mc_pfsl[3] << endl;
    }

  if(nppart < 1)
    {
      take_pevent = kFALSE;
      if(check_cuts)
	cout << "in this event nppart < 1" << endl;
    }

  //find true photon
  Int_t index_true_photon_hadlevel = -1;
  Int_t index_jet = -1;//this variable is to compare with Sanja, Nazar, Ian, Natasha
  if(!Data && (mc_type == "mc_bg_rad")/* || mc_type == "mc_bg_norad")*/)
    {
      for (Int_t i=0; i<Npart; i++)
	{
	  if (Part_id[i]==5)
	    {
	      index_jet=Part_jetid[i]-1;
	      index_true_photon_hadlevel=i;
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
	      index_true_photon_hadlevel=i;
	      break;
	    }
	}
    }
  if(index_true_photon_hadlevel < 0)
    {
      here_is_true_prph = kFALSE;
      if(check_cuts)
	cout << "there is no photon in hadron level" << endl;
    }
  Int_t index_true_photon_partlevel = -1;
/*   for(Int_t i=0; i<nppart; i++) { */
/*     cout << Eventnr << " parton " << i << ": idpart = " << idpart[i] << ", " << ppart[i][0] << " " << ppart[i][1] << " " << ppart[i][2] << " " << ppart[i][3] << endl; */
/*   } */
/*   for(Int_t i=0; i<Npart; i++) { */
/*     cout << Eventnr << " hadron " << i << ": Part_prt = " << Part_prt[i] << ", " << Part_p[i][0] << " " << Part_p[i][1] << " " << Part_p[i][2] << " " << Part_p[i][3] << endl; */
/*   } */
  //find the photon among partons
  for(Int_t i=0; i<nppart; i++) {
    if(idpart[i] == 29) {
      index_true_photon_partlevel = i;
      break;
    }
  } 
  if(index_true_photon_partlevel < 0)
    {
      here_is_true_prph = kFALSE;
      if(check_cuts)
	cout << "there is no photon in parton level" << endl;
    }

  is_true_prph_candidate = here_is_true_prph;
  Int_t index_true_electron_hadlevel = -1;
  Int_t index_true_electron_partlevel = -1;
  Int_t index_photon_vector = -1;
  Int_t index_photon_jet = -1;
  if(take_pevent && here_is_true_prph) {
    //find the electron at hadron level
    for(Int_t i=0; i<Npart; i++)
      {
    	if ((Part_prt[i]==23)||(Part_prt[i]==24))
    	  {
    	    index_true_electron_hadlevel = i;
    	    break;
    	  }
      }
    //find the electron at parton level
    for(Int_t i=0; i<nppart; i++)
      {
    	if ((idpart[i]==23)||(idpart[i]==24))
    	  {
    	    index_true_electron_partlevel = i;
    	    break;
    	  }
      }
/*     cout << "true electron: " << v_true_scattered_electron->Px() << " " */
/* 	 << v_true_scattered_electron->Py() << " " */
/* 	 << v_true_scattered_electron->Pz() << " " */
/* 	 << v_true_scattered_electron->E()  */
/* 	 << ", id = " << Part_id[index_true_electron_hadlevel] */
/* 	 << ", prt = " << Part_prt[index_true_electron_hadlevel]  */
/* 	 << ", e from part_p = " << Part_p[index_true_electron_hadlevel][3] << endl; */
/*     cout << "true photon: " << v_true_prompt_photon->Px() << " " */
/* 	 << v_true_prompt_photon->Py() << " " */
/* 	 << v_true_prompt_photon->Pz() << " " */
/* 	 << v_true_prompt_photon->E() << ", id = " << Part_id[index_true_photon_hadlevel]  */
/* 	 << ", prt = " << Part_prt[index_true_photon_hadlevel] << endl; */

/*     for(Int_t i=0; i<nppart; i++) { */
/*       cout << "parton #" << i << ": " << ppart[i][0] << " " << ppart[i][1] << " " << ppart[i][2] << " " << ppart[i][3] << ", id = " << idpart[i] << endl; */
/*     } */
 
    if(ppart[index_true_photon_partlevel][3] != Part_p[index_true_photon_hadlevel][3]) {
      cerr << "ppart[index_true_photon_partlevel][3] = " << ppart[index_true_photon_partlevel][3] << " != Part_p[index_true_photon_hadlevel][3] = " 
	   << Part_p[index_true_photon_hadlevel][3] << endl;
      here_is_true_prph = kFALSE;
      //      exit(-1);
    }
    //    if(index_true_electron_partlevel < 0) {
    //      cout << "warning: electron not found at parton level" << endl;
    //    }
  }
  vector<KtJet::KtLorentzVector> input_partons;
  vector<KtJet::KtLorentzVector> true_parton_jets;
  //
  // in order of time economy:
  //
  TLorentzVector tl_true_photon; tl_true_photon.SetPxPyPzE(ppart[index_true_photon_partlevel][0],
							   ppart[index_true_photon_partlevel][1],
							   ppart[index_true_photon_partlevel][2],
							   ppart[index_true_photon_partlevel][3]);
  if(tl_true_photon.Et() < 4. || tl_true_photon.Et() > 15. || tl_true_photon.Eta() < -0.7 || tl_true_photon.Eta() > 0.9)
    here_is_true_prph = kFALSE;

  if(take_pevent && here_is_true_prph) {
    //find parton level jet containing true photon
    KtJet::KtLorentzVector r;
    for (Int_t i=0; i<nppart; i++)
      {
	Double_t M2 = ppart[i][3]*ppart[i][3] - ppart[i][2]*ppart[i][2]-ppart[i][1]*ppart[i][1]-ppart[i][0]*ppart[i][0];
	// Skip DIS electron
	if(check_cuts && i==index_true_photon_partlevel) cout << "photon 4-momentum: " 
							      << ppart[i][0] << " "
							      << ppart[i][1] << " "
							      << ppart[i][2] << " "
							      << ppart[i][3] << ", M2 = "
							      << M2 << endl;
							   
	if (i==index_true_electron_partlevel) continue;
	if (TMath::Abs(ppart[i][3]) < 1.e-3) continue;
	if(ppart[i][3] <= 0) continue;
	//	if(M2 < 0 ) continue;
	//	if(TMath::Abs(ppart[i][2]) > ppart[i][3]) continue;
	if (i==index_true_photon_partlevel) index_photon_vector = input_partons.size();
	

	// create a KtJet with TRUE particles and put it onto
	// back of the input_particles vector
	if (TMath::Abs(ppart[i][2])<ppart[i][3])
	  r = KtJet::KtLorentzVector(ppart[i][0],ppart[i][1],ppart[i][2],ppart[i][3]);
	else r = KtJet::KtLorentzVector(ppart[i][0],ppart[i][1],ppart[i][2],
					TMath::Sqrt(ppart[i][0]*ppart[i][0] + ppart[i][1]*ppart[i][1] + ppart[i][2]*ppart[i][2]));
	//	cout << "i: " << r.px() << " " << r.py() << " " << r.pz() << " " << r.e() << " " << r.e() - TMath::Abs(r.pz()) << " " << M2 << endl;
	input_partons.push_back(r);
      }
    double rparameter = 1.0;
    //    cout << "get parton jets" << endl;
    KtJet::KtEvent ev(input_partons,3,2,1,rparameter);
    //
    // Get vector of jets
    //

    true_parton_jets = ev.getJetsEt();
    //    cout << "...OK" << endl;
    //    cout << "what is the problem?" << endl;
    if(check_cuts) {
      cout << "index_true_photon_partlevel = " << index_true_photon_partlevel << endl;
      cout << "index_photon_vector = " << index_photon_vector << endl; 
    }
    for(Int_t i=0; i<true_parton_jets.size(); i++) {
      if(check_cuts) cout << "ktjet #" << i << ": e = " << true_parton_jets[i].e() << ", et = " << true_parton_jets[i].et() << ", eta = " << true_parton_jets[i].eta() << endl;;
      if(true_parton_jets[i].contains(input_partons[index_photon_vector])) {
	index_photon_jet = i;	
	ephoton_over_ejet = input_partons[index_photon_vector].e() / true_parton_jets[i].e();
	if(check_cuts) {
	  cout << "found hadron level jet containing true photon: et = " << true_parton_jets[index_photon_jet].et()
	       << ", eta = " << true_parton_jets[index_photon_jet].eta() << " et_photon = " << input_partons[index_photon_vector].et()
	       << "e_photon = " << input_partons[index_photon_vector].e()
	       << "pz_photon = " << input_partons[index_photon_vector].pz() << endl;
	  cout << "eta_photon = " << input_partons[index_photon_vector].eta() << endl;
	  cout << "OK" << endl;
	}
	if(ephoton_over_ejet > 1 || ephoton_over_ejet < 0) {
	  cout << "strange in Eventnr = " << Eventnr << ": ephoton_over_ejet = " << ephoton_over_ejet << ", e_photon = " << input_partons[index_photon_vector].e() 
	       << ", e_jet = " << true_parton_jets[i].e() << endl;
	  cout << "et_photon = " << input_partons[index_photon_vector].et() 
	       << ", et_jet = " << true_parton_jets[i].et() << endl;
	  cout << "eta_photon = " << input_partons[index_photon_vector].eta() 
	       << ", eta_jet = " << true_parton_jets[i].et() << endl;
	}
	//break;
      }
    }
    //
    // cuts on parton level prph candidate
    //
    if(index_photon_vector >= 0) {

    Double_t true_photon_et  = input_partons[index_photon_vector].et();
    Double_t true_photon_eta = input_partons[index_photon_vector].eta();
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
    if(ephoton_over_ejet < 0.9) {
      is_true_prph_candidate = kFALSE;
      if(check_cuts)
	cout << "rejected by photon/jet energy ratio: " << ephoton_over_ejet << endl;
    }

    //
    // find accompanying hadronic jet
    //
    if(is_true_prph_candidate) {
      //      cout << "these must be equal: " << input_partons[index_photon_vector].e() << " = " << ppart[index_true_photon_partlevel][3] << endl;
      //      cout << "search jets without photon" << endl;
      Double_t max_et_jet = -999.;
      Int_t index_of_accomp_jet = -1;
      for(Int_t i=0; i<true_parton_jets.size(); i++) {
	take_true_jet = kTRUE;
	if(i == index_photon_jet) {
	  take_true_jet = kFALSE;
	}
	// phase space
	if (true_parton_jets[i].et() < ET_JET_CUT) {
	  take_true_jet = kFALSE;
	}  
	if (true_parton_jets[i].eta() < -1.5 || true_parton_jets[i].eta() > 1.8) {
	  take_true_jet = kFALSE;
	}
	if(take_true_jet)
	  {
	    here_is_true_jet = kTRUE;
	    count_nat++;
	    //	    cout << count_nat << " " << Runnr << " " << Eventnr << endl;
	    //	    cout << Siq2el[0] << " " << Mc_x << endl;
	    //	    cout << Mc_pfsl[3] << " " << input_partons[index_photon_vector].et() << " " << true_parton_jets[i].et() << endl;
	    //	    cout << endl;

	    if(true_parton_jets[i].et() > max_et_jet) {
	      max_et_jet = true_parton_jets[i].et();
	      index_of_accomp_jet = i;
	      v_true_parton_acc_jet->SetPxPyPzE(true_parton_jets[i].px(), true_parton_jets[i].py(), true_parton_jets[i].pz(), true_parton_jets[i].e());
	    }
	  }
      } // loop over paroton level jets to find accompanying to photon jet
      if(index_of_accomp_jet >=0)
	here_is_true_jet = kTRUE;
      if(here_is_true_jet) {
	accomp_jet_et = v_true_parton_acc_jet->Et();
	accomp_jet_eta = v_true_parton_acc_jet->Eta();
	
	hist.part_prph_e->Fill(input_partons[index_photon_vector].et(), wtx);
	hist.part_jet_e->Fill(accomp_jet_et, wtx);
	hist.part_Q2->Fill(Mc_q2, wtx);
	hist.part_x->Fill(Mc_x, wtx);
	
	hist.part_cross_et->Fill(input_partons[index_photon_vector].et(), wtx);
	hist.part_cross_eta->Fill(input_partons[index_photon_vector].eta(), wtx);
	hist.part_cross_Q2->Fill(Mc_q2, wtx);
	hist.part_cross_x->Fill(Mc_x, wtx);
	hist.part_cross_et_jet->Fill(accomp_jet_et, wtx);
	hist.part_cross_et_jet2->Fill(accomp_jet_et, wtx);
	hist.part_cross_eta_jet->Fill(accomp_jet_eta, wtx);

	part_et = input_partons[index_photon_vector].et();
	part_eta = input_partons[index_photon_vector].eta();
	part_et_jet = accomp_jet_et;
	part_eta_jet = accomp_jet_eta;
	part_x = Mc_x;
	part_Q2 = Mc_q2;
      }// if(here_is_true_jet
    } //    if(is_true_prph_candidate)
    } //    if(index_photon_vector >= 0)
  } //  if(take_pevent && here_is_true_prph) 
  if((take_pevent && here_is_true_prph && here_is_true_jet && is_true_prph_candidate)) {
    cout << Eventnr << " number of parton level jets: " << true_parton_jets.size() << endl;
    for(Int_t i=0; i<true_parton_jets.size(); i++) {
      cout << "jet #" << i << ": et = " << true_parton_jets[i].et() << ", eta = " << true_parton_jets[i].eta() << endl;
    }    
    cout << "found parton level jet containing true photon: et = " << true_parton_jets[index_photon_jet].et()
	 << ", eta = " << true_parton_jets[index_photon_jet].eta() << " et_photon = " << input_partons[index_photon_vector].et()
	 << ", eta_photon = " << input_partons[index_photon_vector].eta() << endl;
    cout << "parton level Eventr = " << Eventnr << ": q2 = " << Mc_q2 << ", et_photon = " << input_partons[index_photon_vector].et()
	 << ", eta_photon = " << input_partons[index_photon_vector].eta() 
	 << ", et_accomp_jet = " << accomp_jet_et 
	 << ", eta_accomp_jet = " << accomp_jet_eta 
	 << ", ephoton_over_ejet = " << ephoton_over_ejet << endl;

    cout << "================================PARTON LEVEL EVENT " << Eventnr  << " ACCEPTED============================================" << endl;
  }
  check_cuts = kFALSE;
  //  cout << "going out part. selection" << endl;
  return (take_pevent && here_is_true_prph && here_is_true_jet && is_true_prph_candidate);
}
