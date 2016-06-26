#include <iostream>
using namespace std;

#include "TH1D.h"
#include "TMath.h"
#include "selector.h"

void selector::fill_trigger_bits(TH1D* h[], Int_t trigger[])
{
  //  cout << "part1 = " << trigger[0] << " part2 = " << trigger[1] << endl;
  for(Int_t i=0; i<128; i++)
    {
      UInt_t bit_position;
      Int_t bit_number = i+1;
      Bool_t dstw_bit = kTRUE;
      bit_position = (UInt_t) pow((float)2, (float)(i%32));
      if((trigger[i/32] & bit_position) == 0)
	{
	  dstw_bit = kFALSE;
	}
      if(dstw_bit){
	if(i<64)
	  h[0]->Fill(bit_number-0.5);
	else
	  h[1]->Fill(bit_number-0.5);
      }
    }
}
