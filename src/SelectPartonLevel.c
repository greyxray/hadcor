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

Bool_t selector::SelectPartonLevel(Bool_t take_det_event, Bool_t take_had_event = kFALSE)
{
  cout << "...going into part. selection" << endl;
  if (take_had_event) check_cuts = kFALSE;
  // General check of mom-energy conservation on the Parton level:
    Double_t part_px_sum = 0;// = 892.48 * TMath::Sin(2.7e-4);
    Double_t part_py_sum = 0;// = 892.48 * TMath::Sin(1.3e-4);
    Double_t part_pz_sum = 0;// = - 892.48; 
    Double_t part_e_sum  = 0;// = - 27.52 - TMath::Sqrt( 920 * 920 + 0.9382 * 0.9382); //-947.5204783799
    Double_t part_et_sum = 0;// = 0;
    if (check_cuts) cout << "number of partons is " << Nppart << endl;
    if (true || check_cuts) cout << " Index :: Idpart :: Ppart_x :: Ppart_y :: Ppart_z :: Ppart_E " << endl;
    for(Int_t i = 0; i < Nppart; i++)
    {
      part_px_sum += Ppart[i][0];
      part_py_sum += Ppart[i][1];
      part_pz_sum += Ppart[i][2];
      part_e_sum  += Ppart[i][3];
      TLorentzVector v; 
      v.SetPxPyPzE(Ppart[i][0], Ppart[i][1], Ppart[i][2], Ppart[i][3]);
      part_et_sum += v.Et();
      if (true || check_cuts) cout << i << " " << Idpart[i] << " " << Ppart[i][0] << " " <<Ppart[i][1] << " " << Ppart[i][2] << " " << Ppart[i][3]<< endl;
    }
    
  // Additional check
    Double_t m_q2_true = -1. * (bosene * bosene - bospx * bospx - bospy*bospy - bospz*bospz); // parton level
    Double_t m_x_true = m_q2_true / (2. * m_Proton_Energy * (bosene-bospz)); // parton level
    Double_t m_y_true = m_q2_true / (4. * m_Lepton_Energy * m_Proton_Energy * m_x_true); // parton level

    hist.part_q2->Fill(m_q2_true, wtx);
    hist.part_y->Fill(m_y_true, wtx);
    //hist.h2d_etevent_part_had->Fill(part_et_sum, had_et_sum, wtx);
    hist.had_q2->Fill(Mc_q2, wtx);
    hist.had_y->Fill(Mc_y, wtx);

  //  initialisation of some starting parameters
    Bool_t take_pevent = kTRUE;
    Bool_t here_is_true_prph = kTRUE;
    Bool_t is_true_prph_candidate = kTRUE;
    Bool_t here_is_true_jet = kFALSE;
    Bool_t take_true_jet = kTRUE;
    Double_t ephoton_over_ejet = -1.;
    Double_t accomp_jet_et = -1;
    Double_t accomp_jet_eta = -1;

  // electron selection - the same as for hadron level
    if (Mc_q2 < 10.) 
    {
      take_pevent = kFALSE;
      if (check_cuts) cout << "rejected by cut on true q2 = " << Mc_q2 << endl;
    }
    TVector3 v_true_electron(Mc_pfsl[0], Mc_pfsl[1], Mc_pfsl[2]); // momentum of final state lepton
    if ((v_true_electron.Theta()*180.0 / TMath::Pi() < 140.0 ) ||
        (v_true_electron.Theta()*180.0 / TMath::Pi() > 180.0 ))  
    {
      take_pevent = kFALSE;
      if (check_cuts) cout << "rejected by cut on true electron theta = " << v_true_electron.Theta()*180.0/TMath::Pi() << endl;
    }
    if (Mc_pfsl[3] < 10.)
    {
      take_pevent = kFALSE;
      if (check_cuts) cout << "rejected by cut on true electron energy = " << Mc_pfsl[3] << endl;
    }

  //chek validity of parton level MC
    if (Nppart < 1)
    {
      take_pevent = kFALSE;
      if (check_cuts) cout << "in this event Nppart < 1" << endl;
    }

  //Find true photon
  //on the HADRON level --> index_true_hadron_partlevel
    Int_t index_true_photon_hadlevel = -1;
    Int_t index_true_photon_hadlevel_fmckin2 = -1;
    int index_true_electron_partlevel_fmckin2 = -1;
    for(Int_t i = 0; i < Npart; i++)
    {
      if (check_cuts) cout << "Fmckin instance: " << i << " :: Part_id =" << Part_id[i]  << ", Part_prt = "<< Part_prt[i] << endl;
      if (Part_prt[i] == 29)
      {
        index_true_photon_hadlevel = i;
        if (check_cuts) cout << "true photon found: " << index_true_photon_hadlevel << endl;
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
      if (check_cuts) cout << " Index :: Fmck_id :: Fmck_isthep :: Fmck_daug :: Fmck_prt :: Fmck_px :: Fmck_py :: Fmck_pz :: Fmck_e " << endl;
      bool passed_photon = false, passed_electron = false;
        for (int i = 0; i < Fmck_nstor; i++)
        {
          if (Fmck_isthep[i]%10000 == 1 && Fmck_prt[i] == 29) //(Fmck_id[i] == Part_id[index_true_photon_hadlevel]) 
          {
            index_true_photon_hadlevel_fmckin2 = i;
            passed_photon = true;
            v_PrPh.SetPxPyPzE(Fmck_px[i], Fmck_py[i], Fmck_pz[i], Fmck_e[i]);
            if (check_cuts) cout << i << " " << Fmck_id[i] << " " << Fmck_isthep[i] << " " << Fmck_daug[i] << " " << Fmck_prt[i] << " " << Fmck_px[i] << " " <<Fmck_py[i] << " " <<Fmck_pz[i] << " " << Fmck_e[i] <<  "PHOTON" << endl;
            continue;
          }
          // if (!passed_electron && abs(Fmck_px[i] - Mc_pfsl[0]) < 0.01 && 
          //                        abs(Fmck_py[i] - Mc_pfsl[1]) < 0.01 && 
          //                        abs(Fmck_pz[i] - Mc_pfsl[2]) < 0.01 && 
          //                        abs(Fmck_e[i] - Mc_pfsl[3]) < 0.01 ) 
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
          if ( abs(part_e_sum_fmckin2 + v_PrPh.E() + Mc_pfsl[3] + Fmck_e[i] - E_cons) > abs(part_e_sum_fmckin2 + v_PrPh.E() + Mc_pfsl[3] - E_cons) ) 
          {
            break;
          }
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

    if (check_cuts) 
    {
      cout << "Electron: (" << Mc_pfsl[0] << ", " << Mc_pfsl[1] << ", " << Mc_pfsl[2] << ", " << Mc_pfsl[3] << ")" << endl;
      cout << "Photon: (" << v_PrPh.Px() << ", " << v_PrPh.Py() << ", " << v_PrPh.Pz() << ", " << v_PrPh.E() << ")" << endl;
      cout << "Parton E/M-conservation with photon and electron : (" << part_px_sum_fmckin2 + v_PrPh.Px() + Mc_pfsl[0] << ", " << 
                                                        part_py_sum_fmckin2 + v_PrPh.Py() + Mc_pfsl[1] << ", " << 
                                                        part_pz_sum_fmckin2 + v_PrPh.Pz() + Mc_pfsl[2] << ", " << 
                                                        part_e_sum_fmckin2  + v_PrPh.E() + Mc_pfsl[3] << ")" 
                                                        << " and Et = " << part_et_sum << "\n";
      if (part_lev_from_fmckin2) cout << "FMCKIN2 Parton VS QCDPART :  px: " << part_px_sum_fmckin2 << " " << part_px_sum << 
                                     ", py: " << part_py_sum_fmckin2 << " " << part_py_sum << 
                                     ", pz: " << part_pz_sum_fmckin2 << " " << part_pz_sum <<
                                     ", en: " << part_e_sum_fmckin2  << " " << part_e_sum << endl;
    }

    // Check Momentum conservation
    if (check_en_mom_conservation_on_parton_level &&
        check_en_mom_conservation && ( 
             abs(part_e_sum_fmckin2  + v_PrPh.E() + Mc_pfsl[3] - E_cons)  > 1
             || abs(part_pz_sum_fmckin2 + v_PrPh.Pz() + Mc_pfsl[2] - pz_cons) > 1
             //|| abs(part_py_sum_fmckin2 + v_PrPh.Py() + Mc_pfsl[1] - py_cons) > 0.1
             //|| abs(part_px_sum_fmckin2 + v_PrPh.Px() + Mc_pfsl[0] - px_cons) > 0.1 
             ) 
        )
    {
      if (abs(part_e_sum_fmckin2  + v_PrPh.E() + Mc_pfsl[3] - E_cons)  > 1) cout << "E_part broken"<< endl;
      if (abs(part_pz_sum_fmckin2 + v_PrPh.Pz() + Mc_pfsl[2] - pz_cons) > 1 ) cout << "Pz_part broken"<< endl;
      if (abs(part_py_sum_fmckin2 + v_PrPh.Py() + Mc_pfsl[1] - py_cons) > 0.1 ) cout << "Py_part broken"<< endl;
      if (abs(part_px_sum_fmckin2 + v_PrPh.Px() + Mc_pfsl[0] - px_cons) > 0.1 ) cout << "Px_part broken"<< endl;
      en_mom_conservation = false;
      return false;
    }
    //else return true;

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
        if (index_true_electron_partlevel != -1) cout << "Strange index_true_electron_partlevel: " << index_true_electron_partlevel << endl;

      // Additional output
        if (check_cuts) 
        { 
          //v_true_scattered_electron - is filled in SelectHadronLevel
          cout << "true electron: " << v_true_scattered_electron->Px() << " " 
               << v_true_scattered_electron->Py() << " " 
               << v_true_scattered_electron->Pz() << " " 
               << v_true_scattered_electron->E()  
          << ", id = " << Part_id[index_true_electron_hadlevel] 
          << ", prt = " << Part_prt[index_true_electron_hadlevel]  
          << ", E from part_p = " << Part_p[index_true_electron_hadlevel][3] << endl; 
          
          cout << "true photon: " << v_true_prompt_photon->Px() << " " 
               << v_true_prompt_photon->Py() << " " 
               << v_true_prompt_photon->Pz() << " " 
               << v_true_prompt_photon->E() << ", id = " << Part_id[index_true_photon_hadlevel]  
          << ", prt = " << Part_prt[index_true_photon_hadlevel] << endl; 

          // for(Int_t i = 0; i<Nppart; i++)
          //  cout << "parton #" << i << ": " << Ppart[i][0] << " " << Ppart[i][1] << " " << Ppart[i][2] << " " << Ppart[i][3] << ", id = " << Idpart[i] << endl; 
        
        }
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
      // construct vector of r's --> input_partons
      // and find # of parton in input_partons that is true photon --> index_photon_vector
        KtJet::KtLorentzVector r;
        //if (check_cuts) cout << "Nppart " << Nppart<< endl;
        //for (Int_t i = 0; i < Nppart; i++)//Npart
        for (int i = 0; i < parton_lev_prt.size(); i++)
        {
          // Double_t M2 = Ppart[i][3] * Ppart[i][3] 
          //               - Ppart[i][2] * Ppart[i][2] 
          //               - Ppart[i][1] * Ppart[i][1] 
          //               - Ppart[i][0] * Ppart[i][0];// parton mass; 
          double Energy_i = parton_lev_prt[i].E();// Ppart[i][3];
          double Pz_i = parton_lev_prt[i].Pz();// Ppart[i][2]
          double Py_i = parton_lev_prt[i].Py();// Ppart[i][1]
          double Px_i = parton_lev_prt[i].Px();// Ppart[i][0]
          //if (i == index_true_electron_partlevel) continue; // Skip DIS electron
          if (TMath::Abs(Energy_i) < 1.e-3) continue; // Skip low energy
          if (Energy_i <= 0) continue;
          //  if (M2 < 0 ) continue;
          //  if (TMath::Abs(Ppart[i][2]) > Energy_i) continue;
      
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
             
             if (check_cuts) cout << "photon 4-momentum: " 
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
          if (check_cuts) cout << "Now the index_photon_vector should be eq to  Nppart (from QCDPAR block) : " ;
          if (index_photon_vector == Nppart && (check_cuts)) cout << "index_photon_vector == Nppart && (check_cuts): yes" << endl;
          else if (check_cuts) cout << "index_photon_vector == Nppart && (check_cuts): no" << endl;
      // construct vector of jets --> true_parton_jets
        double rparameter = 1.0;
        KtJet::KtEvent ev(input_partons, 3, 2 ,1, rparameter);// pe deltaR e rparameter
        true_parton_jets = ev.getJetsEt();// jets on the parton lev
        
      
      // find # of jet containing photon - index_photon_jet
        for(Int_t i = 0; i < true_parton_jets.size(); i++) 
        {
          if (check_cuts) 
            if (check_cuts) cout << "ktjet #" << i 
                  << ": e = " << true_parton_jets[i].e() 
                  << ", et = " << true_parton_jets[i].et() 
                  << ", eta = " << true_parton_jets[i].eta() 
                  << endl;

          if (true_parton_jets[i].contains(input_partons[index_photon_vector])) 
          {
            index_photon_jet = i; 
            ephoton_over_ejet = input_partons[index_photon_vector].e() / true_parton_jets[i].e();
            //Additional output
              {
                if (check_cuts) 
                  if (check_cuts) cout << "found hadron level jet containing true photon: Et_jet_with_gamma = " << true_parton_jets[index_photon_jet].et()
                        << ", eta_jet_with_gamma = " << true_parton_jets[index_photon_jet].eta() 
                        << " et_photon = " << input_partons[index_photon_vector].et()
                        << "e_photon = " << input_partons[index_photon_vector].e()
                        << "pz_photon = " << input_partons[index_photon_vector].pz()
                        << "eta_photon = " << input_partons[index_photon_vector].eta()
                        << endl;
                if (ephoton_over_ejet > 1 || ephoton_over_ejet < 0) 
                  if (check_cuts) cout << "strange in Eventnr = " << Eventnr 
                      << ": ephoton_over_ejet = " << ephoton_over_ejet 
                      << ", e_photon = " << input_partons[index_photon_vector].e() 
                      << ", e_jet = " << true_parton_jets[i].e() 
                      << "et_photon = " << input_partons[index_photon_vector].et() 
                      << ", et_jet = " << true_parton_jets[i].et() 
                      << "eta_photon = " << input_partons[index_photon_vector].eta() 
                      << ", eta_jet = " << true_parton_jets[i].et() << endl; 
              }
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
              if (check_cuts) cout << "rejected by tue_photon_et = " << true_photon_et << endl;
            }
            if (true_photon_eta < -0.7 || true_photon_eta > 0.9) 
            {
              is_true_prph_candidate = kFALSE;
              if (check_cuts) cout << "rejected by cut on true photon eta = " << true_photon_eta << endl;
            }
            if (ephoton_over_ejet < 0.9) 
            {
              is_true_prph_candidate = kFALSE;
              if (check_cuts) cout << "rejected by photon/jet energy ratio: " << ephoton_over_ejet << endl;
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
                    // Turn on for the check
                    //      cout << count_nat << " " << Runnr << " " << Eventnr << endl;
                    //      cout << Siq2el[0] << " " << Mc_x << endl;
                    //      cout << Mc_pfsl[3] << " " << input_partons[index_photon_vector].et() << " " << true_parton_jets[i].et() << endl;
                    //      cout << endl;

                    if (true_parton_jets[i].et() > max_et_jet) 
                    {
                      max_et_jet = true_parton_jets[i].et();
                      index_of_accomp_jet = i;
                      v_true_parton_acc_jet->SetPxPyPzE(true_parton_jets[i].px(), true_parton_jets[i].py(), true_parton_jets[i].pz(), true_parton_jets[i].e());
                    }
                  }
                } 
              
              //Fill histograms
                if (index_of_accomp_jet >= 0) here_is_true_jet = kTRUE;
                if (here_is_true_jet) 
                {
                  accomp_jet_et = v_true_parton_acc_jet->Et();
                  accomp_jet_eta = v_true_parton_acc_jet->Eta();


                  part_et = input_partons[index_photon_vector].et();
                  part_eta = input_partons[index_photon_vector].eta();
                  part_et_jet = accomp_jet_et;
                  part_eta_jet = accomp_jet_eta;
                  part_x = Mc_x;
                  part_Q2 = Mc_q2;
                  part_xgamma = (input_partons[index_photon_vector].e() - input_partons[index_photon_vector].pz() 
                                + v_true_parton_acc_jet->E() - v_true_parton_acc_jet->Pz() ) / (2. * E_e * Mc_y);
                  if (part_xgamma >= 1) part_xgamma = 0.999;
                  part_xp = (input_partons[index_photon_vector].e() + input_partons[index_photon_vector].pz() 
                                + v_true_parton_acc_jet->E() + v_true_parton_acc_jet->Pz() ) / (2. * E_p);
                  part_dphi = delta_phi(v_true_parton_acc_jet->Phi(), input_partons[index_photon_vector].phi()) * 180.0/TMath::Pi();
                  part_deta = accomp_jet_eta - part_eta;
                  part_dphi_e_ph = delta_phi(v_true_electron.Phi(), input_partons[index_photon_vector].phi()) * 180.0/TMath::Pi();
                  part_deta_e_ph = v_true_electron.Eta() - part_eta;

                  hist.part_prph_e->Fill(part_et, wtx);
                  hist.part_jet_e->Fill(accomp_jet_et, wtx);
                  hist.part_Q2->Fill(part_Q2, wtx);
                  hist.part_x->Fill(part_x, wtx);


                  if (part_x >= 2.e-2) 
                  {
                    if (check_cuts) cout <<"==============>parx_x exceeded upper limit\n";
                    part_x = 18.e-3;
                  }
                  else if (part_x < 2.e-4) 
                  {
                    if (check_cuts) cout <<"==============>parx_x exceeded lower limit\n";
                    part_x = 2.e-4;
                  }

                  if (part_deta >= 2) 
                  {
                    if (check_cuts) cout <<"==============>part_deta exceeded upper limit\n";
                    part_deta = 1.5;
                  }
                  else if (part_deta < -2.2) 
                  {
                    if (check_cuts) cout <<"==============>part_deta exceeded lower limit\n";
                    part_deta = -2.0;
                  }

                  if (part_deta_e_ph >= -0.6) 
                  {
                    if (check_cuts) cout <<"==============>part_deta_e_ph exceeded upper limit\n";
                    part_deta_e_ph = -1.0;
                  }
                  else if (part_deta_e_ph < -3.6) 
                  {
                    if (check_cuts) cout <<"==============>part_deta_e_ph exceeded lower limit\n";
                    part_deta_e_ph = -3.1;
                  }
                  hist.part_cross_et->Fill(part_et, wtx);
                  hist.part_cross_eta->Fill(part_eta, wtx);
                  hist.part_cross_Q2->Fill(part_Q2, wtx);
                  hist.part_cross_x->Fill(part_x, wtx);
                  hist.part_cross_et_jet->Fill(accomp_jet_et, wtx);
                  hist.part_cross_et_jet2->Fill(accomp_jet_et, wtx);
                  hist.part_cross_eta_jet->Fill(accomp_jet_eta, wtx);
                    hist.part_cross_xgamma->Fill(part_xgamma, wtx);
                    hist.part_cross_xp->Fill(part_xp, wtx);
                    hist.part_cross_dphi->Fill(part_dphi, wtx);
                    hist.part_cross_deta->Fill(part_deta, wtx);
                    if (check_cuts) cout << "===============> part_deta=============>" << part_deta << endl;
                    if (check_cuts) cout << "===============> jet=============>" << accomp_jet_eta << endl;
                    if (check_cuts) cout << "===============> prph=============>" << part_eta << endl;
                    hist.part_cross_dphi_e_ph->Fill(part_dphi_e_ph, wtx);
                    hist.part_cross_deta_e_ph->Fill(part_deta_e_ph, wtx);


                }// if (here_is_true_jet
            } //    if (is_true_prph_candidate)
        } 
    } 

  //Additional outup
    if ((take_pevent && here_is_true_prph && here_is_true_jet && is_true_prph_candidate)) 
    {
      if (check_cuts) cout << Eventnr << " number of parton level jets: " 
           << true_parton_jets.size() 
           << endl;

      for(Int_t i = 0; i < true_parton_jets.size(); i++)
        if (check_cuts) cout << "jet #" << i 
            << ": et = " << true_parton_jets[i].et() 
            << ", eta = " << true_parton_jets[i].eta() 
            << endl;  

      if (check_cuts) cout << "found parton level jet containing true photon: et = " << true_parton_jets[index_photon_jet].et()
           << ", eta = " << true_parton_jets[index_photon_jet].eta() 
           << " et_photon = " << input_partons[index_photon_vector].et()
           << ", eta_photon = " << input_partons[index_photon_vector].eta() 
           << endl;

      if (check_cuts) cout << "parton level Eventr = " << Eventnr 
           << ": q2 = " << Mc_q2 
           << ", et_photon = " 
           << input_partons[index_photon_vector].et()
           << ", eta_photon = " << input_partons[index_photon_vector].eta() 
           << ", et_accomp_jet = " << accomp_jet_et 
           << ", eta_accomp_jet = " << accomp_jet_eta 
           << ", ephoton_over_ejet = " << ephoton_over_ejet << endl;

      if (check_cuts) cout << "================================PARTON LEVEL EVENT " << Eventnr  << " ACCEPTED============================================" << endl;
    }
  
  check_cuts = kFALSE;//once is enought
  return (take_pevent && here_is_true_prph && here_is_true_jet && is_true_prph_candidate);
}
