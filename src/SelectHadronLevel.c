#include <iostream>
#include <fstream>
#include <vector>
#include "selector.h"

using namespace std;

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"
using KtJet::KtLorentzVector;
using KtJet::KtEvent;

Bool_t selector::SelectHadronLevel(Bool_t take_det_event)
{
  if (check_cuts)cout << "going into hadr. selection" << endl;
  check_cuts = kFALSE;
  //param init
    Bool_t take_hevent = kTRUE;
    Bool_t here_is_true_prph = kTRUE;
    Bool_t is_true_prph_candidate = kTRUE;
    Bool_t here_is_true_jet = kFALSE;
    Bool_t take_true_jet = kTRUE;
    check_cuts = kTRUE;

    Double_t accomp_jet_et = 0.;
    Double_t accomp_jet_eta = -999.;
    Double_t accomp_jet_phi = -999.;
    hadron_level_jet_cont_photon_index = -1;
    hadron_level_ephoton_over_ejet = -1;
    cout << "param init" << endl;

  hist.mccorel_q2_y->Fill(Mc_q2, Mc_y, wtx);
  hist.mccorel_q2_x->Fill(Mc_q2, Mc_x, wtx);
  hist.mccorel_q2_y_noweight->Fill(Mc_q2, Mc_y);
  hist.mccorel_q2_x_noweight->Fill(Mc_q2, Mc_x);
  if (check_cuts) cout << "check cut on true q2 = " << Mc_q2 << endl;
    if (Mc_q2 < 10.) 
    {
      take_hevent = kFALSE;
      if (check_cuts) cout << "rejected by cut on true electron energy = " << Mc_pfsl[3] << endl;
    }

    //electron sel
      v_true_scattered_electron->SetPxPyPzE(Mc_pfsl[0], Mc_pfsl[1], Mc_pfsl[2], Mc_pfsl[3]); // need to compare with parton level
      TVector3 v_true_electron(Mc_pfsl[0], Mc_pfsl[1], Mc_pfsl[2]);
      if ((v_true_electron.Theta() * 180.0 / TMath::Pi() < 139.8 )||
          (v_true_electron.Theta() * 180.0 / TMath::Pi() > 180.0 ))  
      {
        take_hevent = kFALSE;
        if (check_cuts) cout << "rejected by cut on true electron theta = " << v_true_electron.Theta() * 180.0 / TMath::Pi() << endl;
      }
      if (Mc_pfsl[3] < 10.)
      {
        take_hevent = kFALSE;
        if (check_cuts) cout << "rejected by cut on true electron energy = " << Mc_pfsl[3] << endl;
      }

  // Find true photon on generated level --> here_is_true_prph, index_true_photon, index_jet
    Int_t index_true_photon = -1;
    Int_t index_jet = -1; 
    for(Int_t i = 0; i < Npart; i++)
    {
      cout << "Fmckin instance: " << i << " :: Part_id =" << Part_id[i]  << ", Part_prt = "<< Part_prt[i] << endl;
      if (Part_prt[i] == 29)
      {
        hist.dis_Q2_true->Fill(Mc_q2, wtx);
        hist.dis_electron_e_true->Fill(Mc_pfsl[3], wtx);
        hist.dis_electron_theta_true->Fill(v_true_electron.Theta() * 180.0 / TMath::Pi(), wtx);
        hist.dis_x_true->Fill(Mc_x, wtx);
        hist.dis_y_true->Fill(Mc_y, wtx);
      }
      if (take_hevent) cout << "\t\tpassed electron cuts" << endl;

  //find true photon on gen lev --> here_is_true_prph, index_true_photon, index_jet
     Int_t index_true_photon = -1;
     

    Double_t px_had_sum(0), py_had_sum(0), pz_had_sum(0), E_had_sum(0);
    for(Int_t i = 0; i < Npart; i++)
     {
       px_had_sum += Part_p[i][0];
       py_had_sum += Part_p[i][1];
       pz_had_sum += Part_p[i][2];
       E_had_sum += Part_p[i][3];
       if (check_cuts) cout << "Fmckin instance: " << i << " :: Part_id =" << Part_id[i]  << ", Part_prt = "<< Part_prt[i] << endl;
       if (Part_prt[i] == 29 && index_true_photon == -1)
       {
         index_true_photon = i;
         if (check_cuts) cout << "true photon found: " << index_true_photon << endl;
       }
     }
    if (check_cuts) 
    {
      cout << "Hadron E/M-conservation: (" << px_had_sum << ", " << py_had_sum << ", " << pz_had_sum << ", " << E_had_sum << ")" << endl;
      cout << "HAD: abs(E_had_sum - 947.5) " << abs(E_had_sum - 947.5) << endl;
    }
    if (check_en_mom_conservation && ( abs(E_had_sum  - E_cons)  > 1 ||
             abs(pz_had_sum - pz_cons) > 1 || 
             abs(py_had_sum - py_cons) > 0.1 || 
             abs(px_had_sum - px_cons) > 0.1 ) )
    {
      en_mom_conservation = false;
      return false;
    }
    // TODO: using FMCKin2 check that photon is a daughter of  photon who has no mother 

    if (index_true_photon < 0)
    {
      here_is_true_prph = kFALSE;
      if (check_cuts) cout << "there is no photon on hadron level by selecting by ID" << endl;
    }
    is_true_prph_candidate = here_is_true_prph;

  // Photon cuts
    if (take_hevent && here_is_true_prph) 
    {
      TVector3 v_true_photon(Part_p[index_true_photon][0], Part_p[index_true_photon][1], Part_p[index_true_photon][2]);
               v_true_prompt_photon->SetPxPyPzE(Part_p[index_true_photon][0], Part_p[index_true_photon][1], Part_p[index_true_photon][2], Part_p[index_true_photon][3]);
          
      Double_t true_photon_et  = Part_p[index_true_photon][3] * TMath::Sin(v_true_photon.Theta());
      Double_t true_photon_eta = v_true_photon.Eta();
      if (true_photon_et < 4. || true_photon_et > 15.) 
      {
        is_true_prph_candidate = kFALSE;
        if (check_cuts) cout << "rejected by tue_photon_et = " << true_photon_et << endl;
      }
      if (true_photon_eta < -0.7 || true_photon_eta > 0.9) 
      {
        is_true_prph_candidate = kFALSE;
        if (check_cuts) cout << "rejected by cut on true photon eta = " << true_photon_eta << endl;
      }
    } 
    if (take_hevent && is_true_prph_candidate && check_cuts) cout << "\t\tpassed photon cuts" << endl;
  
  if (take_hevent && is_true_prph_candidate) 
  {
    //parameters init
      vector<KtJet::KtLorentzVector> input_hadrons;
      vector<KtJet::KtLorentzVector> true_jets;
      Int_t index_true_electron = -1;
      Int_t index_photon_vector = -1;
      Int_t index_photon_jet = -1;
      Double_t ephoton_over_ejet = -1;

    //find electron on gen lev
      for(Int_t i = 0; i < Npart; i++)
        if ((Part_prt[i] == 23) || (Part_prt[i] == 24))
      	{
      	  index_true_electron = i;
      	  break;
      	}

    //jets constructions from particles on the hadron level
      KtJet::KtLorentzVector r;
      //construct input_hadrons made of Part_p --> input_hadrons
        for (Int_t i = 0; i < Npart; i++)
        {
          Double_t M2 = Part_p[i][3] * Part_p[i][3] 
                      - Part_p[i][2] * Part_p[i][2]
                      - Part_p[i][1] * Part_p[i][1]
                      - Part_p[i][0] * Part_p[i][0];
    	    // Skip DIS electron
          if (i == index_true_electron) continue;
          if (TMath::Abs(Part_p[i][3]) < 1.e-3) continue;
          if (Part_p[i][3] <= 0) continue;
          //	if (M2 < 0 ) continue;
          //  if (TMath::Abs(Part_p[i][2]) > Part_p[i][3]) continue;
          if (i == index_true_photon) index_photon_vector = input_hadrons.size();

          // create a KtJet with TRUE particles and put it to the end of the input_particles vector
          if (TMath::Abs(Part_p[i][2]) < Part_p[i][3])  r = KtJet::KtLorentzVector(Part_p[i][0], Part_p[i][1], Part_p[i][2], Part_p[i][3]);
                                                   else r = KtJet::KtLorentzVector(Part_p[i][0], Part_p[i][1], Part_p[i][2], 
    				                                                TMath::Sqrt(Part_p[i][0] * Part_p[i][0] 
                                                                      + Part_p[i][1]*Part_p[i][1] 
                                                                      + Part_p[i][2]*Part_p[i][2]));
          // TODO: убрать else
          input_hadrons.push_back(r);
        }

      double rparameter = 1.0;
      KtJet::KtEvent ev(input_hadrons, 3, 2, 1, rparameter);
      true_jets = ev.getJetsEt();
      hadron_level_jets = true_jets; // used in selector.c

    //find phtoton index in vector of jets --> index_photon_jet
      if ( (mc_type=="mc_prph" || mc_type=="mc_prph_rad" ) && index_true_photon < 0)
      {
        if (check_cuts)cout << "Warning in ReClusterize: index of TRUE photon in FMCKIN table < 0."<<endl;
        //	exit(-1);
        //  return jetsEtTRUE.size();
      }
      for (UInt_t i = 0; i < true_jets.size(); i++)
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

    //more cuts on photon
      if (index_photon_jet < 0)  
      { 
        is_true_prph_candidate = kFALSE;
        if (check_cuts) cout<<"Index of photon in jet < 0 - bad!!"<<endl; 
      }
      if (ephoton_over_ejet < 0.9) 
      { 
        is_true_prph_candidate = kFALSE;
        if (check_cuts) cout << "rejected, because ephoton_over_ejet at hadron level is " << ephoton_over_ejet << endl;
      }
    here_is_true_prph = is_true_prph_candidate;

    if (take_hevent && here_is_true_prph) 
    {    
      Double_t max_et_jet = -999.;
      Int_t index_of_accomp_jet = -1;

      // cuts on jet --> v_true_acc_jet, index_of_accomp_jet in true_jets
      for (UInt_t j = 0; j < true_jets.size(); j++)
      {
      	take_true_jet = kTRUE;
      	// cuts on jet
        	if (j == index_photon_jet) 
          {
        	  take_true_jet = kFALSE;
        	  if (check_cuts) cout << "jet " << j << " rejected, because it is a prompt photon jet" << endl;
        	}
        	if (true_jets[j].et() < ET_JET_CUT) 
          {
        	  take_true_jet = kFALSE;
        	  if (check_cuts) cout << "jet " << j << " rejected, because of its et: " << true_jets[j].et()<< endl;
        	}  
        	if (true_jets[j].eta() < -1.5 || true_jets[j].eta() > 1.8) 
          {
        	  take_true_jet = kFALSE;
        	  if (check_cuts) cout << "jet " << j << " rejected, because of its eta: " << true_jets[j].eta()<< endl;
        	}

        //max et jet selected
        	if (take_true_jet)
      	  {
      	    here_is_true_jet = kTRUE;
      	    count_nat++;

      	    if (true_jets[j].et() > max_et_jet) 
            {
      	      max_et_jet = true_jets[j].et();
      	      index_of_accomp_jet = j;
      	      v_true_acc_jet->SetPxPyPzE(true_jets[j].px(), true_jets[j].py(), true_jets[j].pz(), true_jets[j].e());
            }
          }
      }
    if (take_hevent && is_true_prph_candidate && here_is_true_jet) 
    {
      if (check_cuts)cout << "\t\tpassed jets cuts: " << index_of_accomp_jet << endl;
      //exit(-1);
    }

      //Fill histograms
        if (here_is_true_jet && !inclusive_prph_analysis) 
        {
          accomp_jet_et = true_jets[index_of_accomp_jet].et();
          accomp_jet_eta = true_jets[index_of_accomp_jet].eta();
          accomp_jet_phi = true_jets[index_of_accomp_jet].phi();
        }
        if (inclusive_prph_analysis) here_is_true_jet = kTRUE;
        if (take_hevent && here_is_true_prph && here_is_true_jet)
        {
          // Calculating some new variables
            Double_t x_pomeron = v_true_prompt_photon->E() + v_true_prompt_photon->Pz();
            Double_t x_gamma = v_true_prompt_photon->E() - v_true_prompt_photon->Pz();
            Double_t empz_particles_inside_jet = 0.;
            Double_t empz_particles_inside_jet_plus = 0.;

            for(Int_t zloop = 0; zloop < input_hadrons.size(); zloop++)
            { 
              if ( true_jets[index_of_accomp_jet].contains(input_hadrons[zloop]) )//if zufo in the jet of max_et
              {
                x_gamma += input_hadrons[zloop].e() - input_hadrons[zloop].pz();
                x_pomeron += input_hadrons[zloop].e() + input_hadrons[zloop].pz();
                
                empz_particles_inside_jet += input_hadrons[zloop].e() - input_hadrons[zloop].pz();
                empz_particles_inside_jet_plus += input_hadrons[zloop].e() + input_hadrons[zloop].pz();
              }
            }
            x_gamma /= (2. * E_e * Mc_y);
            x_pomeron /= (2. * E_p);
            if (x_gamma >= 1) x_gamma = 0.999;


            Double_t hardest_jet_et = true_jets[index_of_accomp_jet].et();
            Double_t hardest_jet_eta = true_jets[index_of_accomp_jet].eta();
            Double_t hardest_jet_phi = true_jets[index_of_accomp_jet].phi();
            Double_t temp_dphi = delta_phi(hardest_jet_phi, input_hadrons[index_photon_vector].phi()) * 180.0/TMath::Pi(); //? v_true_prompt_photon->Phi())
            Double_t temp_deta = hardest_jet_eta - input_hadrons[index_photon_vector].eta();
            Double_t temp_dphi_e_ph = delta_phi(v_true_electron.Phi(), input_hadrons[index_photon_vector].phi()) * 180.0/TMath::Pi();
            Double_t temp_deta_e_ph = -TMath::Log(TMath::Tan(v_true_electron.Theta()/2.)) - input_hadrons[index_photon_vector].eta();

        	hist.had_prph_e->Fill(input_hadrons[index_photon_vector].et(), wtx);
        	hist.had_jet_e->Fill(accomp_jet_et, wtx);
        	hist.had_Q2->Fill(Mc_q2, wtx);
        	hist.had_x->Fill(Mc_x, wtx);

          if (Mc_x >= 0.02) 
          {
            if (check_cuts)cout <<"==============>Mc_x exceeded upper limit\n";
            Mc_x = 0.015;
          }
          else if (Mc_x < 0.0002) 
          {
            if (check_cuts)cout <<"==============>Mc_x exceeded lower limit\n";
            Mc_x = 0.0005;
          }

                  if (temp_deta_e_ph >= -0.6) 
                  {
                    if (check_cuts)cout <<"==============>temp_deta_e_ph exceeded upper limit\n";
                    temp_deta_e_ph = -1.0;
                  }
                  else if (temp_deta_e_ph < -3.6) 
                  {
                    if (check_cuts)cout <<"==============>temp_deta_e_ph exceeded lower limit\n";
                    temp_deta_e_ph = -3.1;
                  }
          
        	hist.had_cross_et->Fill(input_hadrons[index_photon_vector].et(), wtx);
        	hist.had_cross_eta->Fill(input_hadrons[index_photon_vector].eta(), wtx);
        	hist.had_cross_Q2->Fill(Mc_q2, wtx);
        	hist.had_cross_x->Fill(Mc_x, wtx);
        	hist.had_cross_et_jet->Fill(accomp_jet_et, wtx);
        	hist.had_cross_et_jet2->Fill(accomp_jet_et, wtx);
        	hist.had_cross_eta_jet->Fill(accomp_jet_eta, wtx);
          hist.had_cross_xgamma->Fill(x_gamma, wtx);
          hist.had_cross_xp->Fill(x_pomeron, wtx);
          hist.had_cross_dphi->Fill(temp_dphi, wtx);
          hist.had_cross_deta->Fill(temp_deta, wtx);
          hist.had_cross_dphi_e_ph->Fill(temp_dphi_e_ph, wtx);
          hist.had_cross_deta_e_ph->Fill(temp_deta_e_ph, wtx);

        	had_et = input_hadrons[index_photon_vector].et();
        	had_eta = input_hadrons[index_photon_vector].eta();
        	had_et_jet = accomp_jet_et;
        	had_eta_jet = accomp_jet_eta;
        	had_x = Mc_x;
        	had_Q2 = Mc_q2;

          had_xgamma = (input_hadrons[index_photon_vector].e() - input_hadrons[index_photon_vector].pz() 
                        + true_jets[index_of_accomp_jet].e() - true_jets[index_of_accomp_jet].pz() ) / (2. * E_e * Mc_y);
          if (had_xgamma >=1) had_xgamma = 0.999;
          if (x_gamma != had_xgamma) exit(-1);

          had_xp = (input_hadrons[index_photon_vector].e() + input_hadrons[index_photon_vector].pz() 
                        + true_jets[index_of_accomp_jet].e() + true_jets[index_of_accomp_jet].pz() ) / (2. * E_p);
          had_dphi = delta_phi(true_jets[index_of_accomp_jet].phi(), input_hadrons[index_photon_vector].phi()) * 180.0/TMath::Pi();
          had_deta = true_jets[index_of_accomp_jet].eta() - had_eta;
          had_dphi_e_ph = delta_phi(v_true_electron.Phi(), input_hadrons[index_photon_vector].phi()) * 180.0/TMath::Pi();
          had_deta_e_ph = v_true_electron.Eta() - had_eta;

        	hist.prof_had_cross_et->Fill(input_hadrons[index_photon_vector].et(), input_hadrons[index_photon_vector].et(), wtx);
        	hist.prof_had_cross_eta->Fill(input_hadrons[index_photon_vector].eta(), input_hadrons[index_photon_vector].eta(), wtx);
        	hist.prof_had_cross_Q2->Fill(Mc_q2, Mc_q2, wtx);
        	hist.prof_had_cross_x->Fill(Mc_x, Mc_x, wtx);
        	hist.prof_had_cross_et_jet->Fill(accomp_jet_et, accomp_jet_et, wtx);
        	hist.prof_had_cross_eta_jet->Fill(accomp_jet_eta, accomp_jet_eta, wtx);
          hist.prof_had_cross_xgamma->Fill(x_gamma, x_gamma, wtx);
          hist.prof_had_cross_xp->Fill(x_pomeron, x_pomeron, wtx);
          hist.prof_had_cross_dphi->Fill(temp_dphi, temp_dphi, wtx);
          hist.prof_had_cross_deta->Fill(temp_deta, temp_deta, wtx);
          hist.prof_had_cross_dphi_e_ph->Fill(temp_dphi_e_ph, temp_dphi_e_ph, wtx);
          hist.prof_had_cross_deta_e_ph->Fill(temp_deta_e_ph, temp_deta_e_ph, wtx);

        	if (take_det_event) 
          {
            if(hist.had_cross_xgamma->FindBin(x_gamma) == det_cross_sec_bin_xgamma) //?
            {
              hist.hd_cross_xgamma->Fill(x_gamma, wtx);
              hist.prof_hd_cross_xgamma->Fill(x_gamma, x_gamma, wtx);
            }
            if(hist.had_cross_xp->FindBin(x_pomeron) == det_cross_sec_bin_xp) //?
            {
              hist.hd_cross_xp->Fill(x_pomeron, wtx);
              hist.prof_hd_cross_xp->Fill(x_pomeron, x_pomeron, wtx);
            }
            if(hist.had_cross_dphi->FindBin(temp_dphi) == det_cross_sec_bin_dphi) //?
            {
              hist.hd_cross_dphi->Fill(temp_dphi, wtx);
              hist.prof_hd_cross_dphi->Fill(temp_dphi, temp_dphi, wtx);
            }
            if(hist.had_cross_deta->FindBin(temp_deta) == det_cross_sec_bin_deta) //?
            {
              hist.hd_cross_deta->Fill(temp_deta, wtx);
              hist.prof_hd_cross_deta->Fill(temp_deta, temp_deta, wtx);
            }
            if(hist.had_cross_dphi_e_ph->FindBin(temp_dphi_e_ph) == det_cross_sec_bin_dphi_e_ph) //?
            {
              hist.hd_cross_dphi_e_ph->Fill(temp_dphi_e_ph, wtx);
              hist.prof_hd_cross_dphi_e_ph->Fill(temp_dphi_e_ph, temp_dphi_e_ph, wtx);
            }
            if(hist.had_cross_deta_e_ph->FindBin(temp_deta_e_ph) == det_cross_sec_bin_deta_e_ph) //?
            {
              hist.hd_cross_deta_e_ph->Fill(temp_deta_e_ph, wtx);
              hist.prof_hd_cross_deta_e_ph->Fill(temp_deta_e_ph, temp_deta_e_ph, wtx);
            }
        	  if (hist.had_cross_et->FindBin(input_hadrons[index_photon_vector].et()) == det_cross_sec_bin_et) 
            {
        	    hist.hd_cross_et->Fill(input_hadrons[index_photon_vector].et(), wtx);
        	    hist.prof_hd_cross_et->Fill(input_hadrons[index_photon_vector].et(), input_hadrons[index_photon_vector].et(), wtx);
        	  }
        	  if (hist.had_cross_eta->FindBin(input_hadrons[index_photon_vector].eta()) == det_cross_sec_bin_eta) 
            {
        	    hist.hd_cross_eta->Fill(input_hadrons[index_photon_vector].eta(), wtx);
        	    hist.prof_hd_cross_eta->Fill(input_hadrons[index_photon_vector].eta(), input_hadrons[index_photon_vector].eta(), wtx);
        	  }
        	  if (hist.had_cross_Q2->FindBin(Mc_q2) == det_cross_sec_bin_q2) 
            {
        	    hist.hd_cross_Q2->Fill(Mc_q2, wtx);
        	    hist.prof_hd_cross_Q2->Fill(Mc_q2, Mc_q2, wtx);
        	  }
        	  if (hist.had_cross_x->FindBin(Mc_x) == det_cross_sec_bin_x) 
            {
        	    hist.hd_cross_x->Fill(Mc_x, wtx);
        	    hist.prof_hd_cross_x->Fill(Mc_x, Mc_x, wtx);
        	  }
        	  if (hist.had_cross_et_jet->FindBin(accomp_jet_et) == det_cross_sec_bin_et_jet) 
            {
        	    hist.hd_cross_et_jet->Fill(accomp_jet_et, wtx);
        	    hist.prof_hd_cross_et_jet->Fill(accomp_jet_et, accomp_jet_et, wtx);
        	  }
        	  if (hist.had_cross_eta_jet->FindBin(accomp_jet_eta) == det_cross_sec_bin_eta_jet) 
            {
        	    hist.hd_cross_eta_jet->Fill(accomp_jet_eta, wtx);
        	    hist.prof_hd_cross_eta_jet->Fill(accomp_jet_eta, accomp_jet_eta, wtx);
            }
          }
        } 
    }

    //Additional output
      if ((take_hevent && here_is_true_prph && here_is_true_jet)) 
      {
        if (check_cuts)cout << "hadr level, Eventnr = " << Eventnr << ", N of hadrons = " << Npart << ", N of jets = " << true_jets.size() << ": " << endl;
        for(Int_t i = 0; i < true_jets.size(); i++) 
          if (check_cuts)cout << "ktjet #" << i << ": e = " << true_jets[i].e() << ", et = " << true_jets[i].et() << ", eta = " << true_jets[i].eta() << endl;;
        
        if (check_cuts)cout << "found hadron level jet containing true photon: et = " << true_jets[index_photon_jet].et()
          	 << ", eta = " << true_jets[index_photon_jet].eta() << " et_photon = " << input_hadrons[index_photon_vector].et()
          	 << "e_photon = " << input_hadrons[index_photon_vector].e()
          	 << "pz_photon = " << input_hadrons[index_photon_vector].pz() << endl;
              if (check_cuts)cout << "hadron level Eventnr = " << Eventnr << ": q2 = " << Mc_q2 << ", et_photon = " << input_hadrons[index_photon_vector].et()
          	 << ", eta_photon = " << input_hadrons[index_photon_vector].eta()
          	 << ", et_accomp_jet = " << accomp_jet_et 
          	 << ", eta_accomp_jet = " << accomp_jet_eta 
          	 << ", ephoton_over_ejet = " << ephoton_over_ejet << endl;
        if (check_cuts)cout << "================================HADRON LEVEL EVENT " << Eventnr << " ACCEPTED============================================" << endl;
      }
  }
  return (take_hevent && here_is_true_prph && here_is_true_jet);
}
