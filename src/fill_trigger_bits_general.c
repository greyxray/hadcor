#include <iostream>
using namespace std;

#include "TH1D.h"
#include "TMath.h"
#include "selector.h"

void selector::fill_trigger_bits_general(TH1D* h[], Int_t trigger[], Int_t n)
{
  //  cout << "part1 = " << trigger[0] << " part2 = " << trigger[1] << endl;
  for(Int_t i=0; i<32*n; i++)
    {
      UInt_t bit_position;
      Int_t bit_number = i+1;
      Bool_t bit = kTRUE;
      //      bit_position = (UInt_t) pow(2, i%32);
      bit_position = (UInt_t) pow((float)2, (float)(i%32));
      if((trigger[i/32] & bit_position) == 0)
	{
	  bit = kFALSE;
	}
      if(bit){
	h[i/32]->Fill(bit_number-0.5);
      }
    }
}
