#define selectorTermninate_c

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#include "constants.h"
#include "selector.h"
#include "runinfo.h"

void selector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  delete v_uncorr_prompt_photon;
  delete v_corr_prompt_photon;
  cout << "generating tree..." << endl;
  fChain->SetEventList(event_list);
  TTree *tree = fChain->CopyTree("");
  tree->Write();
  cout << "done" << endl;
  f->Write();
  {
    cout << "DIS03: " << cut_DIS03 << endl;
  }
  cout << "gamma had cut: " << cut_gammahad << endl;
  cout << "number of good tracks: " << cut_ngoodtrks << endl;
  cout << "electron is in outer parto of RCAL: " << cut_RCAL << endl;
  cout << "poltake: " << cut_poltake << endl;
  cout << "QED Compton rejecting: " << cut_QEDC << endl;
  cout << "y*(1-x)^2 > 0.004: " << cut_yx2 << endl;
  //		cout << "Sincand: " << cut_Sincand << endl;
  //		cout << "Siecorr && Emecorr: " << cut_Siecorr << endl;
  cout << "electron candidate: " << cutElCand << endl;
  cout << "electron energy: " << cutEC_1 << endl;
  cout << "electron track in ctd acc.: " << cutEC_6_1 << endl;
  cout << "electron track momentum: " << cut_track_momentum << endl;
  cout << "electron track dca: " << cut_electron_track_dca << endl;
  cout << "electron distance to the module edge: " << cut_electron_distance_module_edge << endl;
  cout << "electron distance to beam pipe: " << cutEC_4 << endl;
  cout << "enin: " << cutEC_5 << endl;
  cout << "zvtx: " << cutC_1 << endl;
  cout << "e-pz: " << cutC_2 << endl;
  cout << "q2: " << cutPS_1 << endl;
  cout << "y: " << cutPS_2 << endl;
  cout << "pt/et: " << cutC_3 << endl;
  cout << "chimney: " << cutEC_2 << endl;
  cout << "super crack: " << cutEC_3 << endl;
  
  //		cout << "Rejected cells: " << float(badCells/numberOfEntries) <<
  //				" per event" << endl;
  cout << "Total events selected: " << selectedEventsCount << endl;
  cout << "If in ctd acc: " << cutCTDacc << endl;
  cout << "Selected jets in this file: " << localSelectedJetsCount << endl;
  cout << "selected jets total: " << totalSelectedJetsCount << endl;
  cout << "eta_lab: "<< cutJPS_1 << endl;
  cout << "dijets assimetric cut (is already not assymetric): " << cut_assymetric << endl;
  cout << "Mjj > 20 GeV cut: " << cut_Mjj << endl;
  cout << "distance betw. el. and jet: "<< cutJC_1 << endl;
  cout << "E Breit > 5 && eta lab < -2: " << cutJC_2 << endl;
  cout << "EtLab > 3: " << cutJC_3 << endl;
  cout << "dijet events: " << detector_dijets << endl;
  
  if(!Data)
    {
      cout << "hadron q2: " << cut_q2_bos << endl;
      cout << "hadron y: " << cut_y_bos << endl;
      cout << "\n-----parton level-----" << endl;
      cout << "Et_jet: " << cut_part_Et_jet << endl;
      cout << "M_jj cut: " << cut_p_Mjj << endl;
      cout << "eta_jet: " << cut_part_eta << endl;
      cout << "parton dijet events: " << partonDijetEvents << endl;
      cout << "\n-----hadron level-----" << endl;
      cout << "Et_jet: " << cut_hadr_Et_jet << endl;
      cout << "M_jj cut: " << cut_h_Mjj << endl;
      cout << "eta_jet: " << cut_hadr_eta << endl;
      cout << "hadron dijet events: " << hadronDijetEvents << endl;
    }
  cout << "period: " << period << endl;
  cout << "use_corr: " << use_corr << endl;
  cout << endl;
  //		} //if % 10
  
  
  cout << "Total processed: " << processedEventsCount << " events" << endl;
  cout << take_det_dijet << endl;

  time_end.Set();
  cout << "\nbegin time:\t" << time_begin.AsString()  << endl;
  cout << "finish time:\t" << time_end.AsString() << endl << endl;
  Int_t d_days = time_end.GetDay() - time_begin.GetDay();
  Int_t d_hours = time_end.GetHour() - time_begin.GetHour();
  Int_t d_minutes = time_end.GetMinute() - time_begin.GetMinute();
  Int_t d_sec = time_end.GetSecond() - time_begin.GetSecond();

  if (d_sec < 0)
    {
      d_minutes--;
      d_sec = 60 + d_sec;
    }
  if (d_minutes < 0)
    {
      d_hours--;
      d_minutes = 60 + d_minutes;
    }
  if (d_hours < 0)
    {
      d_days--;
      d_hours = 24 + d_hours;
    }
  cout << "time of running: " << d_days << " days " << d_hours << " hours "
       << d_minutes << " minutes " << d_sec <<" seconds" << endl;

}
