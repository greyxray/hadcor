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

Bool_t selector::SelectPartonLevel(Bool_t take_det_event, Bool_t take_had_event = kFALSE)
{
  cout << "...going into part. selection" << endl;
  if (take_had_event) check_cuts = kTRUE;

  // General check of mom-energy conservation on the hadron level:
    Double_t had_px_sum = 892.48 * TMath::Sin(2.7e-4);
    Double_t had_py_sum = 892.48 * TMath::Sin(1.3e-4);
    Double_t had_pz_sum = - 892.48;
    Double_t had_e_sum = - 27.52 - TMath::Sqrt( 920 * 920 + 0.9382 * 0.9382);
    Double_t had_et_sum = 0;
    //  cout << "in SelectPartons: number of hadrons is " << Npart << endl;
    for(Int_t i = 0; i < Npart; i++)
    {
      had_px_sum += Part_p[i][0];
      had_py_sum += Part_p[i][1];
      had_pz_sum += Part_p[i][2];
      had_e_sum += Part_p[i][3];
      TLorentzVector v; 
      v.SetPxPyPzE(Part_p[i][0], Part_p[i][1], Part_p[i][2], Part_p[i][3]);
      had_et_sum += v.Et();
    }
    /*  
      hist.had_sum_e->Fill(had_e_sum, wtx);
      hist.had_sum_px->Fill(had_px_sum, wtx);
      hist.had_sum_py->Fill(had_py_sum, wtx);
      hist.had_sum_pz->Fill(had_pz_sum, wtx);
      hist.had_sum_et->Fill(had_et_sum, wtx);
      cout << "event #" << Eventnr << ": " << had_e_sum << " " << had_px_sum << " " << had_py_sum << " " << had_pz_sum << " " << had_et_sum << endl;
    */

  // General check of mom-energy conservation on the parton level:
    Double_t part_px_sum = 892.48 * TMath::Sin(2.7e-4);
    Double_t part_py_sum = 892.48 * TMath::Sin(1.3e-4);
    Double_t part_pz_sum = - 892.48;
    Double_t part_e_sum = - 27.52 - TMath::Sqrt( 920 * 920 + 0.9382 * 0.9382);
    Double_t part_et_sum = 0;
    //  cout << "number of partons is " << Nppart << endl;
    for(Int_t i = 0; i < Nppart; i++)
    {
      // hist.part_Idpart->Fill(Idpart[i], wtx);
      part_px_sum += Ppart[i][0];
      part_py_sum += Ppart[i][1];
      part_pz_sum += Ppart[i][2];
      part_e_sum += Ppart[i][3];
      TLorentzVector v; 
      v.SetPxPyPzE(Ppart[i][0], Ppart[i][1], Ppart[i][2], Ppart[i][3]);
      part_et_sum += v.Et();
    }
    /*  
      hist.part_sum_e->Fill(part_e_sum, wtx);
      hist.part_sum_px->Fill(part_px_sum, wtx);
      hist.part_sum_py->Fill(part_py_sum, wtx);
      hist.part_sum_pz->Fill(part_pz_sum, wtx);
      hist.part_sum_et->Fill(part_et_sum, wtx);
      hist.part_Nppart->Fill(Nppart, wtx);
    */
    //This is a check - not crusial
      Double_t m_q2_true = -1. * (bosene * bosene - bospx * bospx - bospy*bospy - bospz*bospz); // parton level
      Double_t m_x_true = m_q2_true / (2. * m_Proton_Energy * (bosene-bospz)); // parton level
      Double_t m_y_true = m_q2_true / (4. * m_Lepton_Energy * m_Proton_Energy * m_x_true); // parton level

      hist.part_q2->Fill(m_q2_true, wtx);
      hist.part_y->Fill(m_y_true, wtx);
      hist.h2d_etevent_part_had->Fill(part_et_sum, had_et_sum, wtx);
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
    TVector3 v_true_electron(Mc_pfsl[0], Mc_pfsl[1], Mc_pfsl[2]); // Four-momentum of final state lepton
    if ((v_true_electron.Theta()*180.0/TMath::Pi() < 139.8 ) ||
        (v_true_electron.Theta()*180.0/TMath::Pi() > 180.0 ))  
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

  /* additional filling before
    if (take_pevent) 
    {
      hist.dis_Q2_true_part->Fill(Mc_q2, wtx); 
      hist.dis_electron_e_true_part->Fill(Mc_pfsl[3], wtx); 
      hist.dis_electron_theta_true_part->Fill(v_true_electron.Theta()*180.0/TMath::Pi(), wtx); 
      hist.dis_x_true_part->Fill(Mc_x, wtx); 
      hist.dis_y_true_part->Fill(Mc_y, wtx); 
    }
  */

  //find true photon
  //on the HADRON level --> index_true_hadron_partlevel
    Int_t index_true_photon_hadlevel = -1;
    Int_t index_jet = -1;//this variable is to compare with Sanja, Nazar, Ian, Natasha
    
    // if (!Data)
    // {
    //   if (mc_type == "mc_bg_rad")/* || mc_type == "mc_bg_norad")*/
    //     for (Int_t i = 0; i < Npart; i++)
    //       if (Part_id[i] == 5)
    //       {
    //         index_jet = Part_jetid[i] - 1;
    //         index_true_photon_hadlevel = i;
    //         break;
    //       }
    //   else if (mc_type == "mc_prph")
    //     for (Int_t i = 0; i < Npart; i++)
    //       if (Part_id[i] == 12)
    //       {
    //         index_jet = Part_jetid[i] - 1;
    //         index_true_photon_hadlevel = i;
    //         break;

    //       }
    // }
    // if (index_true_photon_hadlevel < 0)
    // {
    //   here_is_true_prph = kFALSE;
    //   if (check_cuts) cout << "there is no photon in hadron level" << endl;
    // }
     //New
    for(Int_t i = 0; i < Npart; i++)
     {
       cout << "Fmckin instance: " << i << " :: Part_id =" << Part_id[i]  << ", Part_prt = "<< Part_prt[i] << endl;
       if (Part_prt[i] == 29)
       {
         index_true_photon_hadlevel = i;
         cout << "true photon found: " << index_true_photon_hadlevel << endl;
         break;
       }
     }

  //find true photon 
  //on the PARTON level --> index_true_photon_partlevel
    Int_t index_true_photon_partlevel = -1;

    //find the photon among partons
     // No such instance in the MC - will force to put gamma from hadron level
      // for(Int_t i = 0; i < Nppart; i++) 
      // {
      //   if (Idpart[i] == 29) //fmcprt -> gamma
      //   {
      //     index_true_photon_partlevel = i;
      //     break;
      //   }
      // } 
      // if (index_true_photon_partlevel < 0)
      // {
      //   here_is_true_prph = kFALSE;
      //   if (check_cuts) cout << "there is no photon in parton level" << endl;
      // }

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

      /*     
        v_true_scattered_electron - is filled in SelectHadronLevel
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

        for(Int_t i=0; i<Nppart; i++)
         cout << "parton #" << i << ": " << Ppart[i][0] << " " << Ppart[i][1] << " " << Ppart[i][2] << " " << Ppart[i][3] << ", id = " << Idpart[i] << endl; 
      */

      // ENERGY OF PHOTON_PARTON  != ENERGY OF PHOTON_HADRON => NO TRUE PRPH
       // No such instance in the MC - will force to put gamma from hadron level
        // if (Ppart[index_true_photon_partlevel][3] != Part_p[index_true_photon_hadlevel][3]) 
        // {
        //   cerr << "Ppart[index_true_photon_partlevel][3] = " << Ppart[index_true_photon_partlevel][3] 
        //        << " != Part_p[index_true_photon_hadlevel][3] = " << Part_p[index_true_photon_hadlevel][3] 
        //        << endl;
        //   here_is_true_prph = kFALSE;
        //   // exit(-1);
        // }
    }

  // PrPh cuts on parton level(some):
    
    TLorentzVector tl_true_photon; //prph on the parton level
    // No such instance in the MC - will force to put gamma from hadron level
    // tl_true_photon.SetPxPyPzE(Ppart[index_true_photon_partlevel][0],
  		// 					   Ppart[index_true_photon_partlevel][1],
  		// 					   Ppart[index_true_photon_partlevel][2],
  		// 					   Ppart[index_true_photon_partlevel][3]);
    tl_true_photon.SetPxPyPzE(Part_p[index_true_photon_hadlevel][0],
                              Part_p[index_true_photon_hadlevel][1],
                              Part_p[index_true_photon_hadlevel][2],
                              Part_p[index_true_photon_hadlevel][3]);
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
        cout << "Nppart " << Nppart<< endl;
        for (Int_t i = 0; i < Npart; i++)
        {
          Double_t M2 = Ppart[i][3] * Ppart[i][3] 
                        - Ppart[i][2] * Ppart[i][2] 
                        - Ppart[i][1] * Ppart[i][1] 
                        - Ppart[i][0] * Ppart[i][0];// parton mass
    					   
    	    if (i == index_true_electron_partlevel) continue; // Skip DIS electron
    	    if (TMath::Abs(Ppart[i][3]) < 1.e-3) continue; // Skip low energy
    	    if (Ppart[i][3] <= 0) continue;
    	    //	if (M2 < 0 ) continue;
    	    //	if (TMath::Abs(Ppart[i][2]) > Ppart[i][3]) continue;

    	    if (i == index_true_photon_partlevel) index_photon_vector = input_partons.size();
    	
        	if (TMath::Abs(Ppart[i][2]) < Ppart[i][3]) r = KtJet::KtLorentzVector(Ppart[i][0], Ppart[i][1], Ppart[i][2], Ppart[i][3]);
          else  r = KtJet::KtLorentzVector(Ppart[i][0], Ppart[i][1], Ppart[i][2],
    					     TMath::Sqrt(Ppart[i][0] * Ppart[i][0] 
                              + Ppart[i][1] * Ppart[i][1] 
                              + Ppart[i][2] * Ppart[i][2]));
          input_partons.push_back(r);
        }


        //add the missing photon to the parton level instances
          {
            Double_t M2 = Part_p[index_true_photon_hadlevel][3] * Part_p[index_true_photon_hadlevel][3] 
                          - Part_p[index_true_photon_hadlevel][2] * Part_p[index_true_photon_hadlevel][2] 
                          - Part_p[index_true_photon_hadlevel][1] * Part_p[index_true_photon_hadlevel][1] 
                          - Part_p[index_true_photon_hadlevel][0] * Part_p[index_true_photon_hadlevel][0];// parton mass
             
                     cout << "photon 4-momentum: " 
                          << Part_p[index_true_photon_hadlevel][0] << " "
                          << Part_p[index_true_photon_hadlevel][1] << " "
                          << Part_p[index_true_photon_hadlevel][2] << " "
                          << Part_p[index_true_photon_hadlevel][3] << ", M2 = "
                          << M2 << endl;
                   
            if (index_true_photon_hadlevel != index_true_electron_partlevel
                && TMath::Abs(Part_p[index_true_photon_hadlevel][3]) >= 1.e-3
                && Part_p[index_true_photon_hadlevel][3] > 0) 
            {

              index_photon_vector = input_partons.size();
          
              if (TMath::Abs(Part_p[index_true_photon_hadlevel][2]) < Part_p[index_true_photon_hadlevel][3]) r = KtJet::KtLorentzVector(Part_p[index_true_photon_hadlevel][0], Part_p[index_true_photon_hadlevel][1], Part_p[index_true_photon_hadlevel][2], Part_p[index_true_photon_hadlevel][3]);
              else  r = KtJet::KtLorentzVector(Part_p[index_true_photon_hadlevel][0], Part_p[index_true_photon_hadlevel][1], Part_p[index_true_photon_hadlevel][2],
                       TMath::Sqrt(Part_p[index_true_photon_hadlevel][0] * Part_p[index_true_photon_hadlevel][0] 
                                  + Part_p[index_true_photon_hadlevel][1] * Part_p[index_true_photon_hadlevel][1] 
                                  + Part_p[index_true_photon_hadlevel][2] * Part_p[index_true_photon_hadlevel][2]));
              input_partons.push_back(r);
            }
          }
          cout << "Now the index_photon_vector should be eq to  Nppart (from QCDPAR block) : " ;
          if (index_photon_vector == Nppart) cout << " yes" << endl;
            else cout << "no" << endl;
      // construct vector of jets --> true_parton_jets
        double rparameter = 1.0;
        KtJet::KtEvent ev(input_partons, 3, 2 ,1, rparameter);// pe deltaR e rparameter
        true_parton_jets = ev.getJetsEt();// jets on the parton lev
        //if (check_cuts)  cout << "index_true_photon_partlevel = " << index_true_photon_partlevel << endl << "index_photon_vector = " << index_photon_vector << endl; 
      
      // find # of jet containing photon - index_photon_jet
        for(Int_t i = 0; i < true_parton_jets.size(); i++) 
        {
          if (check_cuts) 
            cout << "ktjet #" << i 
                  << ": e = " << true_parton_jets[i].e() 
                  << ", et = " << true_parton_jets[i].et() 
                  << ", eta = " << true_parton_jets[i].eta() 
                  << endl;

          if (true_parton_jets[i].contains(input_partons[index_photon_vector])) 
          {
            index_photon_jet = i;	
            ephoton_over_ejet = input_partons[index_photon_vector].e() / true_parton_jets[i].e();
            //some output
              {
                if (check_cuts) 
                  cout << "found hadron level jet containing true photon: Et_jet_with_gamma = " << true_parton_jets[index_photon_jet].et()
                        << ", eta_jet_with_gamma = " << true_parton_jets[index_photon_jet].eta() 
                        << " et_photon = " << input_partons[index_photon_vector].et()
                        << "e_photon = " << input_partons[index_photon_vector].e()
                        << "pz_photon = " << input_partons[index_photon_vector].pz()
                        << "eta_photon = " << input_partons[index_photon_vector].eta()
                        << endl;
                if (ephoton_over_ejet > 1 || ephoton_over_ejet < 0) 
                  cout << "strange in Eventnr = " << Eventnr 
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
            //      cout << "these must be equal: " << input_partons[index_photon_vector].e() << " = " << Ppart[index_true_photon_partlevel][3] << endl;
          
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
                    //	    cout << count_nat << " " << Runnr << " " << Eventnr << endl;
                    //	    cout << Siq2el[0] << " " << Mc_x << endl;
                    //	    cout << Mc_pfsl[3] << " " << input_partons[index_photon_vector].et() << " " << true_parton_jets[i].et() << endl;
                    //	    cout << endl;

                    if (true_parton_jets[i].et() > max_et_jet) 
                    {
                      max_et_jet = true_parton_jets[i].et();
                      index_of_accomp_jet = i;
                      v_true_parton_acc_jet->SetPxPyPzE(true_parton_jets[i].px(), true_parton_jets[i].py(), true_parton_jets[i].pz(), true_parton_jets[i].e());
                    }
                  }
                } 
              
              //fill some histograms
                if (index_of_accomp_jet >= 0) here_is_true_jet = kTRUE;
                if (here_is_true_jet) 
                {
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
                }// if (here_is_true_jet
            } //    if (is_true_prph_candidate)
        } 
    } 

  //some outup
    if ((take_pevent && here_is_true_prph && here_is_true_jet && is_true_prph_candidate)) 
    {
      cout << Eventnr << " number of parton level jets: " 
           << true_parton_jets.size() 
           << endl;

      for(Int_t i = 0; i < true_parton_jets.size(); i++)
        cout << "jet #" << i 
            << ": et = " << true_parton_jets[i].et() 
            << ", eta = " << true_parton_jets[i].eta() 
            << endl;  

      cout << "found parton level jet containing true photon: et = " << true_parton_jets[index_photon_jet].et()
  	       << ", eta = " << true_parton_jets[index_photon_jet].eta() 
           << " et_photon = " << input_partons[index_photon_vector].et()
  	       << ", eta_photon = " << input_partons[index_photon_vector].eta() 
           << endl;

      cout << "parton level Eventr = " << Eventnr 
           << ": q2 = " << Mc_q2 
           << ", et_photon = " 
           << input_partons[index_photon_vector].et()
  	       << ", eta_photon = " << input_partons[index_photon_vector].eta() 
        	 << ", et_accomp_jet = " << accomp_jet_et 
        	 << ", eta_accomp_jet = " << accomp_jet_eta 
        	 << ", ephoton_over_ejet = " << ephoton_over_ejet << endl;

      cout << "================================PARTON LEVEL EVENT " << Eventnr  << " ACCEPTED============================================" << endl;
    }
  
  check_cuts = kFALSE;//once is enought
  return (take_pevent && here_is_true_prph && here_is_true_jet && is_true_prph_candidate);
}
