void test()
{
  //  TFile *f = new TFile("/data/zenith234b/kuprash/prompt_photons/orange/oran_kt_parton/run/mc_bg_04-05e_1.root");
  TChain *ch = new TChain("h1");
  ch->Add("/data/zenith234b/kuprash/prompt_photons/orange/oran_kt_parton/run/mc_bg_04-05e_1.root");
  //  ch->Add("/data/zenith234b/zhmak/Topic/orange/nat_out/mc_bg/more_MC_2010aug/mc_bg_04-05e/mc_bg_04-05e_1.root");
  Float_t Zvtx;
  Int_t nppart;
  Int_t Npart;
  ch->SetBranchAddress("Zvtx", &Zvtx);
  ch->SetBranchAddress("npart", &Npart);
  ch->SetBranchAddress("nppart", &nppart);

  Int_t nentries = ch->GetEntries();
  for(Int_t i=0; i<10; i++) {
    ch->GetEntry(i);
    cout << "event " << i << ", zvtx = " << Zvtx << ", nppart = " << nppart << ", Npart = " << Npart << endl;
  }
}
