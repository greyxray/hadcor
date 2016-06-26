#include <stdlib.h>
#include "TH1D.h"
#include <iostream>
#include "runinfo.h"
using std::cout;
using std::endl;

Bool_t poltake(Int_t Runnr, Int_t Eventnr, Int_t Lpoltake, Int_t Tpoltake, RunInfo runinfo[], const Int_t runinfo_runs_number)
{
  if(Runnr > 1 && Runnr > 37715)
    {
      Bool_t foundpoltakerun = kFALSE;
      Bool_t PolTake = kFALSE;
      Bool_t zerolumis = kFALSE;
      for(size_t i=0; i<runinfo_runs_number;i++)
        {
          if(Runnr == runinfo[i].GetRunnr())
            {
              foundpoltakerun = kTRUE;
              if(runinfo[i].GetTpol_lumi() > 0 || runinfo[i].GetLpol_lumi() > 0)
                {
                  if(runinfo[i].GetTpol_lumi() > runinfo[i].GetLpol_lumi())
                    {
                      if(Tpoltake>0)
                        {
                          PolTake = kTRUE;
                        }
                    }
                  else
                    {
                      if(Lpoltake>0)
                        {
                          PolTake = kTRUE;
                        }
                    }
                }
              else
                {
                  zerolumis = kTRUE;
                }
              break;
            }
        }
      if(!foundpoltakerun)
        {
          cout << "the run " << Runnr << " was not found in the pol lumi array" << endl;
          exit(-1);
        }
      if(!PolTake)
        {
          return kFALSE;
        }
      if(PolTake && zerolumis)
        {
          cout << "the run " << Runnr << " has no tpol and lpol lumi and poltake is 1 in event "<< Eventnr << "!" << endl;
          exit(-1);
        }
    }
  else
    return kTRUE;
}
