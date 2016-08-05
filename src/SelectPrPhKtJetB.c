#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#include "selector.h"

extern "C"
{
  float pcele_(int &Imod,float &Emes,float &Emip,float &The,float &Phi);
}

Bool_t selector::SelectPrPhKtJetB(Int_t i_true_rad_photon,  Int_t electron_in_zufo)
{
  //for MC we find the radiated_candidate_number as the closest to i_true_rad_photon - found by ID
  cout << "\033[1;31mSelectPrPhKtJetB\033[0m\n" << endl;

//FMCKin[0] -- electron
//FMCKin[1] -- photon

    cout << "\t\tzufo island containing the DIS electron: " << Fmc_ezisl << endl;
    cout << "\t\tzufo island containing the radiative photon: " << Fmc_gzisl << endl;
    cout << "\t\tzufo containing the DIS electron: " << Fmc_ezufo << endl;
    cout << "\t\tzufo containing the radiative photon: " << Fmc_gzufo << endl;
    cout << "\t\tzufo island number for electron from FMCKin: " << Fmc_zisl[0] << endl;
    cout << "\t\tzufo number for photon from FMCKin: " << Fmc_zufo[1] << endl;
    if(!Data && i_true_rad_photon > -1) 
    {
      if (check_cuts) cout << "photon : " << Part_p[i_true_rad_photon][0] << " " << Part_p[i_true_rad_photon][1] << " " << Part_p[i_true_rad_photon][2] << endl;
      TVector3 v_true_rad_photon(Part_p[i_true_rad_photon][0], Part_p[i_true_rad_photon][1], Part_p[i_true_rad_photon][2]);
      Double_t dr, dr_min = 999.;

      cout<< "Nzufos = "<< Nzufos<< endl;
      cout<< electron_in_zufo<< "electron_in_zufo" << endl;
      for(Int_t zloop = 0; zloop < Nzufos; zloop++)
      {
        if (zloop == electron_in_zufo) 
        {
          cout << "zloop" << zloop << " excluded as electron"<< endl;
          continue;
        }
        Double_t zufo_px = Zufo[zloop][0], 
                 zufo_py = Zufo[zloop][1], 
                 zufo_pz = Zufo[zloop][2];
        if (TMath::Sqrt( (zufo_px*zufo_px + zufo_py*zufo_py) / (zufo_pz*zufo_pz) ) < 1.e-10 
                        || zufo_px!=zufo_px || zufo_py!=zufo_py || zufo_pz!=zufo_pz 
                        || zufo_px*zufo_px+zufo_py*zufo_py+zufo_pz*zufo_pz<1.e-200)
        {
          cout << "zloop" << zloop << " excluded by position"<< endl;
          continue;
        }

        TVector3 v_zufo(zufo_px, zufo_py, zufo_pz);

        dr = v_zufo.DeltaR(v_true_rad_photon);
        cout << zloop << "# " << dr << endl; 
        if(dr < dr_min && dr < 0.2) 
        {
          dr_min = dr;
          radiated_candidate_number = zloop;
        }
      }

      // OLd not working and commented
      /*
        for(Int_t zloop = 0; zloop < Knzufos; zloop++)
        {
          Double_t zufo_px = Kzufos[zloop][0], zufo_py = Kzufos[zloop][1], zufo_pz = Kzufos[zloop][2];
          if (TMath::Sqrt( (zufo_px*zufo_px + zufo_py*zufo_py) / (zufo_pz*zufo_pz) ) < 1.e-10 
                          || zufo_px!=zufo_px || zufo_py!=zufo_py || zufo_pz!=zufo_pz 
                          || zufo_px*zufo_px+zufo_py*zufo_py+zufo_pz*zufo_pz<1.e-200)
              continue;

          TVector3 v_zufo(zufo_px, zufo_py, zufo_pz);

          dr = v_zufo.DeltaR(v_true_rad_photon);
          cout << zloop << "# " << dr << endl; 
          if(dr < dr_min && dr < 0.2) 
          {
            dr_min = dr;
            radiated_candidate_number = zloop;
          }
        }
      */

      cout << "\033[1m\033[33mEnd of Nzufos loop. radiated_candidate_number = " << radiated_candidate_number << ", dr_min = " << dr_min << "\033[0m\n"  << endl;
        if (radiated_candidate_number == electron_in_zufo)
        {
          cout << "radiated_candidate_number == electron_in_zufo" << endl;
          return kFALSE;
        }
        else 
        {
          cout << "radiated_candidate_number != electron_in_zufo" << endl;
        }
    }

  //
  cout << "Find zufo of prph"<< endl;
  Bool_t take_prph = kTRUE;

  //constructing vector on/of zufos without e and photon - input_vec, input_vec_to_zufos
          vector<KtLorentzVector> input_vec;
          vector<Int_t> input_vec_to_zufos;
          for(Int_t zloop = 0; zloop < Nzufos; zloop++)
          {
            if (electron_in_zufo > -1 && zloop == electron_in_zufo)
            {
              if (check_cuts) cout << "rejecting zufo number " << electron_in_zufo << " with e = " << Zufo[zloop][3] << " because it is most probably electron " << endl;
              continue;
            }
            
            KtLorentzVector p(Zufo[zloop][0], Zufo[zloop][1], Zufo[zloop][2], Zufo[zloop][3]);
            input_vec.push_back(p);
            input_vec_to_zufos.push_back(zloop);// check if this is 1 2 3 4 5 ...

          }
          if (check_cuts) cout << Nzufos << " " << input_vec.size() << " hadrons added" << endl;
        
        //constructing jets
        int type  = 3; // pe
        int angle = 2; // deltaR
        int recom = 1; // !pt
        double rparameter = 1.0;

        // Construct the KtEvent object 
        KtEvent ev(input_vec, type, angle, recom, rparameter);//Missed compare to full code
        vector<KtLorentzVector> jets = ev.getJetsEt();
  cout << "Nzufos " << Nzufos << endl;
  for (Int_t zloop = 0; zloop < Nzufos; zloop++)
  {
    if (zloop == electron_in_zufo) 
    {
      cout << "zloop" << zloop << " excluded as electron"<< endl;
      continue;
    }
    take_prph = kTRUE;
    cout << "zloop " << zloop << ": Nzufos " << Nzufos<< endl;

    //if (!Data && mc_type == "mc_bg_rad") cout << Runnr << " " << Eventnr << " " << zloop << " " << radiated_candidate_number << endl;
    if (!Data && mc_type == "mc_bg_rad"   && radiated_candidate_number != zloop)    take_prph = kFALSE;
    if (!Data && mc_type == "mc_bg_norad" && radiated_candidate_number == zloop)    take_prph = kFALSE;
    // sets global cuts
      //parameters initialising
        TLorentzVector z; 
                       z.SetPxPyPzE(Zufo[zloop][0], Zufo[zloop][1], Zufo[zloop][2], Zufo[zloop][3]);
        cout <<"4-vector: " << Zufo[zloop][0] <<" "<< Zufo[zloop][1] <<" "<< Zufo[zloop][2] <<" "<< Zufo[zloop][3] << endl;
        Float_t  zu_phi      = z.Phi();//Kzufophi[zloop];
        Double_t zu_eta      = z.Eta();//Kzufoeta[zloop];
        Double_t zu_pt       = TMath::Sqrt(Zufo[zloop][0] * Zufo[zloop][0] + Zufo[zloop][1] * Zufo[zloop][1]);
        Float_t  zu_theta    = TMath::ATan2(zu_pt, Zufo[zloop][2]);
        TVector3 v_photon; 
                 v_photon.SetXYZ(Zufo[zloop][0], Zufo[zloop][1], Zufo[zloop][2]);
        Float_t photon_e = Zufo[zloop][3];
        Float_t bpres_mips = (Float_t)BPRES_mips(v_photon, bpres_conerad);

        //cout << "calling correction " << pcele_mode << endl;
        //cout<< " " << photon_e << " " << bpres_mips << " " << zu_theta << " " << zu_phi << endl;
        // cout << "before pcele: " << 0. << " " << pcele_mode << " " << photon_e << " " << bpres_mips
        //       	   << " " << zu_theta << " " << zu_phi  << " " << photon_e * TMath::Sin(zu_theta) << endl;
        Double_t zu_e_corr = pcele_(pcele_mode, photon_e, bpres_mips, zu_theta, zu_phi);
        // cout << "after pcele: " << zu_e_corr  << " " << pcele_mode << " " << photon_e << " " << bpres_mips 
        //       	   << " " << zu_theta << " " << zu_phi << " " << zu_e_corr * TMath::Sin(zu_theta) << endl;
        Double_t zu_et_corr  = zu_e_corr * TMath::Sin(zu_theta);
              // cout << "zu_et_corr = " << zu_et_corr << endl;
        Double_t zu_et_uncorr = Zufo[zloop][3] * TMath::Sin(zu_theta);
        // OLd not working and commented Double_t zu_fmax     = Kzufofmax[zloop];//Calibration(Kzufofmax[zloop],"zu","fmax",zu_eta,zu_et_corr);
        // OLd not working and commentedif (zu_fmax > 0.999) zu_fmax=0.999;
        // OLd not working and commentedDouble_t zu_deltaz   = Kzufodeltaz[zloop]/5.45;//Calibration(Kzufodeltaz[zloop]/5.45,"zu","deltaz",zu_eta,zu_et_corr);
        //      myKzufodeltaz[zloop] = Kzufodeltaz[zloop]/5.45;
        if (check_cuts) cout << "prph candidate " << zloop  << " of " << Nzufos << ", eta = " << zu_eta  << ", et_corr = " << zu_et_corr << endl;
        // Cleaning cuts
        //      if (zu_fmax  <=0.05)  {
        //	take_prph = kFALSE;
        //	if(check_cuts)
        //	  cout << "candidate " << zloop << " rejected by cut on zu_fmax = " << zu_fmax << endl;
        //      }
        //      if (zu_deltaz <0.0)   {
        //     take_prph = kFALSE;
        //    }

      //cuts
        // OLd not working and commented if (zu_fmax    <= 0.05)        take_prph = kFALSE;
        // OLd not working and commented if (zu_deltaz  <  0.0)         take_prph = kFALSE; 
        if (zu_eta     <  -0.7)   
        {
          
          take_prph = kFALSE;
        }// Phase space cuts
        if (zu_eta     >  0.9)    
        {
          take_prph = kFALSE;
        }
        if (zu_et_corr <  4.)     
        {
          take_prph = kFALSE;
        }
        if (zu_et_corr >  15.)    
        {
          take_prph = kFALSE;
        }
        if (Zufoeemc[zloop] / Zufoecal[zloop]<  0.9)
        {
          take_prph = kFALSE;
        }//  Photon EMC Energy / Photon Total Energy > 0.90
        if (Tufo[zloop][0] != 31) 
        {
          take_prph = kFALSE;
        }// Only trackless zufos
    
        // track isolation Int_t ntrack=Tracks(zu_phi, zu_eta, 2, min_dist);  
          Int_t ntrack = 0;
          cout << "Trk_ntracks:" << Trk_ntracks << endl;
          for (Int_t i = 0; i < Trk_ntracks; i++)
          {
            TVector3 tr;	      
            tr.SetXYZ(Trk_px[i], Trk_py[i], Trk_pz[i]);

            if (tr.Mag()<0.25) continue;

            // Calculating phi, eta, distance of track
              Double_t ddr = (tr.PseudoRapidity()-zu_eta)*(tr.PseudoRapidity()-zu_eta);
              Double_t delta_phi = tr.Phi() - zu_phi;
              if(delta_phi > TMath::Pi())  delta_phi -= 2.*TMath::Pi();
              if(delta_phi < -TMath::Pi()) delta_phi += 2.*TMath::Pi();
              ddr = TMath::Sqrt(ddr + delta_phi * delta_phi);
              if (ddr < 0.2)      ntrack++;

            if(check_cuts) cout << "track # " << i << " " << tr.Eta() << " " << tr.PseudoRapidity() <<  " " << zu_eta << ", phi: " << tr.Phi() << " " << zu_phi << " " << ddr << " " << ntrack << endl;
          }
          if (ntrack > 0)      take_prph = kFALSE;

        // Fraction of zufo energy to jet energy
        // OLd not working and commented   Double_t zu_jetEnergy = Kpjets[ Kzufoidjet[zloop] - 1][3];
          Int_t candidate_jet_num = -1;  
          cout << jets[0] << endl;
          for(int jloop = 0; jloop < jets.size(); jloop++)
          {
            int pos = find(input_vec_to_zufos.begin(), input_vec_to_zufos.end(), zloop) - input_vec_to_zufos.begin();
            if(pos >= input_vec_to_zufos.size()) {
                cout << "SOMETHING WENT WEONG" << endl;
            }
            if( jets[jloop].contains(input_vec[input_vec_to_zufos[pos]]) )
            {
              candidate_jet_num = jloop; 
              cout << "\t\tcandidate_jet_num" << candidate_jet_num << endl;
              break;            
            }
          }
          cout << "candidate_jet_num " << candidate_jet_num << endl;
          Double_t zu_jetEnergy = jets[candidate_jet_num].et();
          cout << Zufo[zloop][3] << "/" << zu_jetEnergy << endl;
          if (Zufo[zloop][3] / zu_jetEnergy < 0.9)
          {
            take_prph = kFALSE;
            cout << "rejected by Zufo[zloop][3]/jets" << endl;
          }
                       
        //dR cut:
          TVector3 jet_contained_prph(jets[candidate_jet_num].px(), 
                                      jets[candidate_jet_num].py(),
                                      jets[candidate_jet_num].pz());

          Double_t  dR_jet = (jet_contained_prph.Eta() - zu_eta) * (jet_contained_prph.Eta() - zu_eta);
          Double_t dphi_jet = jet_contained_prph.Phi() - zu_phi;
          if (dphi_jet > TMath::Pi())  dphi_jet = dphi_jet - 2.*TMath::Pi();
          if (dphi_jet < -TMath::Pi()) dphi_jet = dphi_jet + 2.*TMath::Pi();
          dR_jet = TMath::Sqrt(dR_jet + dphi_jet*dphi_jet);
          if (check_cuts) cout << "cont jet: " << jet_contained_prph.Phi() << " " << jet_contained_prph.Eta() << ", photon: " << zu_phi<< " " << zu_eta  << " dR = " << dR_jet << endl;
          if(dR_jet > 0.2)
          {
            take_prph = kFALSE;
            if(check_cuts) cout << "candidate with eta " << zu_eta << " rejected by dR_jet = " << dR_jet << endl;
          }

    //if this is prph assighn index and properties
      if(take_prph)
      {
        here_is_prph = kTRUE;
        if (zu_et_corr > max_et_candidate)
        {
          jet_energy_frac         = Zufo[zloop][3] / zu_jetEnergy;
          cell_energy_frac        = Kzufoemcfrac[zloop];
          max_et_candidate        = zu_et_corr;
          max_et_candidate_uncorr = zu_et_uncorr;
          max_et_candidate_number = zloop;
          //candidate_jet_number = jloop;
          v_uncorr_prompt_photon->SetPxPyPzE(Zufo[zloop][0], 
                                             Zufo[zloop][1], 
                                             Zufo[zloop][2], 
                                             Zufo[zloop][3]);
          Double_t corfac = zu_e_corr / Zufo[zloop][3];
          v_corr_prompt_photon->SetPxPyPzE(corfac * Zufo[zloop][0],
                                           corfac * Zufo[zloop][1],
                                           corfac * Zufo[zloop][2],
                                           corfac * Zufo[zloop][3]);
          v_prompt_photon_jet->SetPxPyPzE(Kpjets[Kzufoidjet[zloop]-1][0],
                                          Kpjets[Kzufoidjet[zloop]-1][1],
                                          Kpjets[Kzufoidjet[zloop]-1][2],
                                          Kpjets[Kzufoidjet[zloop]-1][3]);
          glob_fmax   = Kzufofmax[zloop];
          glob_deltaz = myKzufodeltaz[zloop];
          glob_fmax   = stretch_calib(Kzufofmax[zloop], "zu", "fmax", v_corr_prompt_photon->Eta(), v_corr_prompt_photon->Et());
          glob_deltaz = stretch_calib(myKzufodeltaz[zloop], "zu", "deltaz", v_corr_prompt_photon->Eta(), v_corr_prompt_photon->Et());
        }
        if(check_cuts)
        cout << "====\n HERE IS PRPH CANDIDATE WITH ETA = " << zu_eta << ", it is zufo # " << max_et_candidate_number << "\n====" << endl;
        phot_count++;
      }
    

    // To select the most energetic zufo of this event (in case if more than 1 zufo)
    /*  
      if (zu_et_corr<zu_et_event){
      take_prph = kFALSE;
      }
      zu_et_event=zu_et_corr;
      index_RECO_phot = i;

      // Filling gamma candidate variables
      ZuCand.InitNull();
      ZuCand.SetE(Kzufos[zloop][3]);
      ZuCand.momentum.SetXYZ(Kzufos[zloop][0],Kzufos[zloop][1],Kzufos[zloop][2]);
      //    cout<<ZuCand.momentum.Perp()<<"  "<<zu_pt<<endl;
      ZuCand.SetECorr(zu_e_corr);
      ZuCand.SetEtCorr(zu_et_corr);
      ZuCand.SetPt(zu_pt);
      ZuCand.SetEta(zu_eta);
      ZuCand.SetTheta(zu_theta);
      ZuCand.SetFmax(zu_fmax);
      ZuCand.SetDeltaz(zu_deltaz);
      ZuCand.SetPhi(zu_phi);
      ZuCand.SetBPREmips(emip);
      ZuCand.SetEMCfrac(Kzufoemcfrac[zloop]);
      ZuCand.SetZufotype(Kzufotype[zloop]);
      ZuCand.SetJetEnergy(zu_jetEnergy);
      ZuCand.SetJetFrac(ZuCand.e/ZuCand.jetEnergy);
      ZuCand.SetTracks(ntrack);
      ZuCand.SetMinDist(min_dist);
    */    
    phot_count++; // no sence - will be equal to nzufos
  }

  return here_is_prph;
}
