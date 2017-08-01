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

Double_t selector::GetXgamma(bool check_cuts = false)
{   
  Double_t w_part_xgamma = -1;
  //  initialisation of some starting parameters
    Bool_t take_pevent = kTRUE;
    Bool_t here_is_true_prph = kTRUE;
    Bool_t is_true_prph_candidate = kTRUE;
    Bool_t here_is_true_jet = kFALSE;
    Bool_t take_true_jet = kTRUE;
    Double_t ephoton_over_ejet = -1.;
    Double_t accomp_jet_et = -1;
    Double_t accomp_jet_eta = -1;

  cout << "GetXgamma::Mc_q2:: "<<  Mc_q2<< endl;
  // electron selection - the same as for hadron level
    if (Mc_q2 < 10.) 
    {
      take_pevent = kFALSE;
      if (check_cuts) cout << "GetXgamma::rejected by cut on true q2 = " << Mc_q2 << endl;
    }
    if (Mc_q2 < q2_cut_low || Mc_q2 > q2_cut_high) 
    {
      if (check_cuts && !nodebugmode) cout << "GetXgamma::PARTON LEVEL: special cut q2: "<< Mc_q2<< "< q2_cut_low ||  > q2_cut_high"<< endl;
      take_pevent = kFALSE;
    }
    TVector3 v_true_electron(Mc_pfsl[0], Mc_pfsl[1], Mc_pfsl[2]); // momentum of final state lepton
    if ((v_true_electron.Theta()*180.0 / TMath::Pi() < 140.0 ) ||
        (v_true_electron.Theta()*180.0 / TMath::Pi() > 180.0 ))  
    {
      take_pevent = kFALSE;
      if (check_cuts) cout << "GetXgamma::rejected by cut on true electron theta = " << v_true_electron.Theta()*180.0/TMath::Pi() << endl;
    }
    if (Mc_pfsl[3] < 10.)
    {
      take_pevent = kFALSE;
      if (check_cuts) cout << "GetXgamma::rejected by cut on true electron energy = " << Mc_pfsl[3] << endl;
    }

  //chek validity of parton level MC
    if (Nppart < 1)
    {
      take_pevent = kFALSE;
      if (check_cuts) cout << "GetXgamma::in this event Nppart < 1" << endl;
    }

  //Find true photon
  //on the HADRON level --> index_true_hadron_partlevel
    Int_t index_true_photon_hadlevel = -1;
    Int_t index_true_photon_hadlevel_fmckin2 = -1;
    int index_true_electron_partlevel_fmckin2 = -1;
    for(Int_t i = 0; i < Npart; i++)
    {
      if (check_cuts) cout << "GetXgamma::Fmckin instance: " << i << " :: Part_id =" << Part_id[i]  << ", Part_prt = "<< Part_prt[i] << endl;
      if (Part_prt[i] == 29)
      {
        index_true_photon_hadlevel = i;
        if (check_cuts) cout << "GetXgamma::true photon found: " << index_true_photon_hadlevel << endl;
        break;
      }
    }
    TLorentzVector v_PrPh;
    v_PrPh.SetPxPyPzE(Part_p[index_true_photon_hadlevel][0], Part_p[index_true_photon_hadlevel][1], Part_p[index_true_photon_hadlevel][2], Part_p[index_true_photon_hadlevel][3]);

    vector<TLorentzVector> parton_lev_prt;
    double part_px_sum_fmckin2 = 0;
    double part_py_sum_fmckin2 = 0;
    double part_pz_sum_fmckin2 = 0;
    double part_e_sum_fmckin2 = 0;
    if (part_lev_from_fmckin2)
    {
      bool passed_photon = false, passed_electron = false;
        for (int i = 0; i < Fmck_nstor; i++)
        {
          if (Fmck_isthep[i]%10000 == 1 && Fmck_prt[i] == 29)
          {
            index_true_photon_hadlevel_fmckin2 = i;
            passed_photon = true;
            v_PrPh.SetPxPyPzE(Fmck_px[i], Fmck_py[i], Fmck_pz[i], Fmck_e[i]);
            if (check_cuts) cout << i << " " << Fmck_id[i] << " " << Fmck_isthep[i] << " " << Fmck_daug[i] << " " << Fmck_prt[i] << " " << Fmck_px[i] << " " <<Fmck_py[i] << " " <<Fmck_pz[i] << " " << Fmck_e[i] <<  "PHOTON" << endl;
            continue;
          }
          if (Fmck_isthep[i]%10000 == 1 && (Fmck_prt[i] == 23 || Fmck_prt[i] == 24))
          {
            passed_electron = true; 
            index_true_electron_partlevel_fmckin2 = i;
            if (check_cuts) cout << i << " " << Fmck_id[i] << " " << Fmck_isthep[i] << " " << Fmck_daug[i] << " " << Fmck_prt[i] << " " << Fmck_px[i] << " " <<Fmck_py[i] << " " <<Fmck_pz[i] << " " << Fmck_e[i] <<  "ELECTRON" << endl;
            continue;
          }      
          if (!passed_photon) 
          {
            if (check_cuts) cout << i << " " << Fmck_id[i] << " " << Fmck_isthep[i] << " " << Fmck_daug[i] << " " << Fmck_prt[i] << " " << Fmck_px[i] << " " <<Fmck_py[i] << " " <<Fmck_pz[i] << " " << Fmck_e[i] <<  "BEFORE PHOTON" << endl;
            continue;
          };

          TLorentzVector v; 
          if ( abs(part_e_sum_fmckin2 + v_PrPh.E() + Mc_pfsl[3] + Fmck_e[i] - E_cons) > abs(part_e_sum_fmckin2 + v_PrPh.E() + Mc_pfsl[3] - E_cons) ) break;
          v.SetPxPyPzE(Fmck_px[i], Fmck_py[i], Fmck_pz[i], Fmck_e[i]);
          parton_lev_prt.push_back(v);
          part_px_sum_fmckin2 += Fmck_px[i];
          part_py_sum_fmckin2 += Fmck_py[i];
          part_pz_sum_fmckin2 += Fmck_pz[i];
          part_e_sum_fmckin2 += Fmck_e[i];
          if (check_cuts) cout << i << " " << Fmck_id[i] << " " << Fmck_isthep[i] << " " << Fmck_daug[i] << " " << Fmck_prt[i] << " " << Fmck_px[i] << " " <<Fmck_py[i] << " " <<Fmck_pz[i] << " " << Fmck_e[i] <<  "TAKEN" << endl;
          if (Fmck_prt[i] == 2092) break;
        }
    }
    else 
    {
      for (int i = 0; i < Nppart; i++)
      {
        TLorentzVector v; 
        v.SetPxPyPzE(Ppart[i][0], Ppart[i][1], Ppart[i][2], Ppart[i][3]);
        parton_lev_prt.push_back(v);
        part_px_sum_fmckin2 += Ppart[i][0];
        part_py_sum_fmckin2 += Ppart[i][1];
        part_pz_sum_fmckin2 += Ppart[i][2];
        part_e_sum_fmckin2 += Ppart[i][3];
      }
    }

  //Find true photon 
  //on the PARTON level --> index_true_photon_partlevel
    Int_t index_true_photon_partlevel = -1;

    //parameters init 
      is_true_prph_candidate = here_is_true_prph;// True if on had and part lev
      Int_t index_true_electron_hadlevel = -1;
      Int_t index_true_electron_partlevel = -1;
      Int_t index_photon_vector = -1;
      Int_t index_photon_jet = -1;

  // finding e and additional energy check for prph
    if (take_pevent && here_is_true_prph) 
    {
      //find the electron at hadron level
        for(Int_t i = 0; i < Npart; i++)
          if ((Part_prt[i] == 23) || (Part_prt[i] == 24))// e- e+
          {
            index_true_electron_hadlevel = i;
            break;
          }

      //find the electron at parton level
        for(Int_t i = 0; i < Nppart; i++)
          if ((Idpart[i]==23) || (Idpart[i]==24))// e- e+
          {
            index_true_electron_partlevel = i;
            break;
          }
        if (index_true_electron_partlevel != -1) cout << "GetXgamma::Strange index_true_electron_partlevel: " << index_true_electron_partlevel << endl;
    }

  // PrPh cuts on parton level(some):
    TLorentzVector tl_true_photon; //prph on the parton level
    tl_true_photon.SetPxPyPzE(v_PrPh.Px(),
                              v_PrPh.Py(),
                              v_PrPh.Pz(),
                              v_PrPh.E());
    if (tl_true_photon.Et() < 4. || tl_true_photon.Et() > 15. 
        || tl_true_photon.Eta() < -0.7 || tl_true_photon.Eta() > 0.9)
      here_is_true_prph = kFALSE;

  // finding jets on parton level, cuts on photon
    vector<KtJet::KtLorentzVector> input_partons;
    vector<KtJet::KtLorentzVector> true_parton_jets;
    if (take_pevent && here_is_true_prph) 
    {
        KtJet::KtLorentzVector r;
        for (int i = 0; i < parton_lev_prt.size(); i++)
        {
          double Energy_i = parton_lev_prt[i].E();
          double Pz_i = parton_lev_prt[i].Pz();
          double Py_i = parton_lev_prt[i].Py();
          double Px_i = parton_lev_prt[i].Px();
          
          if (TMath::Abs(Energy_i) < 1.e-3) continue; // Skip low energy
          if (Energy_i <= 0) continue;
      
          if (TMath::Abs(Pz_i) < Energy_i) r = KtJet::KtLorentzVector(Px_i, Py_i, Pz_i, Energy_i);
          else  r = KtJet::KtLorentzVector(Px_i, Py_i, Pz_i, TMath::Sqrt(Px_i * Px_i + Py_i * Py_i + Pz_i * Pz_i));
          input_partons.push_back(r);
        }


        //add the missing photon to the parton level instances
          {
            Double_t M2 = v_PrPh.E() * v_PrPh.E() 
                          - v_PrPh.Pz() * v_PrPh.Pz() 
                          - v_PrPh.Py() * v_PrPh.Py() 
                          - v_PrPh.Px() * v_PrPh.Px();// parton mass
             
             if (check_cuts) cout << "GetXgamma::photon 4-momentum: " 
                  << v_PrPh.Px() << " "
                  << v_PrPh.Py() << " "
                  << v_PrPh.Pz() << " "
                  << v_PrPh.E() << ", M2 = "
                  << M2 << endl;
                   
            if (index_true_photon_hadlevel != index_true_electron_partlevel
                && TMath::Abs(v_PrPh.E()) >= 1.e-3
                && v_PrPh.E() > 0) 
            {

              index_photon_vector = input_partons.size();
          
              if (TMath::Abs(v_PrPh.Pz()) < v_PrPh.E()) r = KtJet::KtLorentzVector(v_PrPh.Px(), v_PrPh.Py(), v_PrPh.Pz(), v_PrPh.E());
              else  r = KtJet::KtLorentzVector(v_PrPh.Px(), v_PrPh.Py(), v_PrPh.Pz(),
                        TMath::Sqrt(v_PrPh.Px() * v_PrPh.Px() 
                                  + v_PrPh.Py() * v_PrPh.Py() 
                                  + v_PrPh.Pz() * v_PrPh.Pz()));
              input_partons.push_back(r);
            }
          }
          if (check_cuts) cout << "GetXgamma::Now the index_photon_vector should be eq to  Nppart (from QCDPAR block) : " ;
          if (index_photon_vector == Nppart && (check_cuts)) cout << "GetXgamma::index_photon_vector == Nppart && (check_cuts): yes" << endl;
          else if (check_cuts) cout << "GetXgamma::index_photon_vector == Nppart && (check_cuts): no" << endl;
        
        // construct vector of jets --> true_parton_jets
          double rparameter = 1.0;
          KtJet::KtEvent ev(input_partons, 3, 2 ,1, rparameter);// pe deltaR e rparameter
          true_parton_jets = ev.getJetsEt();// jets on the parton lev
        
      // find # of jet containing photon - index_photon_jet
        for(Int_t i = 0; i < true_parton_jets.size(); i++) 
        {
          if (check_cuts) 
            if (check_cuts) cout << "GetXgamma::ktjet #" << i 
                  << ": e = " << true_parton_jets[i].e() 
                  << ", et = " << true_parton_jets[i].et() 
                  << ", eta = " << true_parton_jets[i].eta() 
                  << endl;

          if (true_parton_jets[i].contains(input_partons[index_photon_vector])) 
          {
            index_photon_jet = i; 
            ephoton_over_ejet = input_partons[index_photon_vector].e() / true_parton_jets[i].e();
          }
        }
      
      // cuts on photon on part level, find accomp jet on parton level
        if (index_photon_vector >= 0) 
        {
          Double_t true_photon_et  = input_partons[index_photon_vector].et();
          Double_t true_photon_eta = input_partons[index_photon_vector].eta();
          
          //cuts on photon on the parton level
            if (true_photon_et < 4. || true_photon_et > 15.) 
            {
              is_true_prph_candidate = kFALSE;
              if (check_cuts) cout << "GetXgamma::rejected by tue_photon_et = " << true_photon_et << endl;
            }
            if (true_photon_eta < -0.7 || true_photon_eta > 0.9) 
            {
              is_true_prph_candidate = kFALSE;
              if (check_cuts) cout << "GetXgamma::rejected by cut on true photon eta = " << true_photon_eta << endl;
            }
            if (ephoton_over_ejet < 0.9) 
            {
              is_true_prph_candidate = kFALSE;
              if (check_cuts) cout << "GetXgamma::rejected by photon/jet energy ratio: " << ephoton_over_ejet << endl;
            }
          
          // find accompanying jet --> index_of_accomp_jet and v_true_parton_acc_jet
            if (is_true_prph_candidate) 
            {
              // find accompanying hadronic jet - index_of_accomp_jet and v_true_parton_acc_jet
                Double_t max_et_jet = -999.;
                Int_t index_of_accomp_jet = -1;
                for(Int_t i = 0; i < true_parton_jets.size(); i++) 
                {
                  take_true_jet = kTRUE;
                  if (i == index_photon_jet) take_true_jet = kFALSE;
                  if (true_parton_jets[i].et() < ET_JET_CUT) take_true_jet = kFALSE;
                  if (true_parton_jets[i].eta() < -1.5 || true_parton_jets[i].eta() > 1.8) take_true_jet = kFALSE;
                  
                  if (take_true_jet)
                  {
                    here_is_true_jet = kTRUE;
                    count_nat++;


                    if (true_parton_jets[i].et() > max_et_jet) 
                    {
                      max_et_jet = true_parton_jets[i].et();
                      index_of_accomp_jet = i;
                      v_true_parton_acc_jet->SetPxPyPzE(true_parton_jets[i].px(), true_parton_jets[i].py(), true_parton_jets[i].pz(), true_parton_jets[i].e());
                    }
                  }
                } 

                w_part_xgamma = (input_partons[index_photon_vector].e() - input_partons[index_photon_vector].pz() 
                              + v_true_parton_acc_jet->E() - v_true_parton_acc_jet->Pz() ) / (2. * E_e * Mc_y);
                cout << "GetXgamma::Check intermetiate:" << input_partons[index_photon_vector].e() << " - " << input_partons[index_photon_vector].pz()  <<
                              " + " << v_true_parton_acc_jet->E() << " - " << v_true_parton_acc_jet->Pz() << " / 2 * " << E_e << " * " <<Mc_y << "\n";
                if (w_part_xgamma >= 1) w_part_xgamma = 0.999;
            }
        }
    }
  return w_part_xgamma;                
}