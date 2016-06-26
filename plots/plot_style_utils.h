void sign_window(TVirtualPad *p, TH2D* h, TString xaxis, TString yaxis, TString title)
{
  Double_t height = 1-p->GetTopMargin()-p->GetBottomMargin();
  Double_t width = 1-p->GetLeftMargin()-p->GetRightMargin();
  //  cout << "pad width is: " << width << endl;
  //  cout << "pad height is: " << height << endl;
  if((width == 0.) || (height == 0.)) {
    cout << "In sign_window(). Width or height of pad is 0. Cannot set text size." << endl;
    exit(-1);
  }
  Double_t fontsizex = 0.08/width; //0.06
  Double_t fontsizey = 0.08/height;

  h->SetStats(kFALSE);
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xaxis);
  h->GetYaxis()->SetTitle(yaxis);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(fontsizex);
  h->GetXaxis()->SetTitleFont(42);
  //  h->GetYaxis()->SetTitleSize(0.5*fontsizey);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelSize(fontsizex);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetLabelSize(fontsizey);
  h->GetYaxis()->SetLabelFont(42);
  if(p->GetLogy())
        h->GetYaxis()->SetTitleOffset(2.5);
  else
        h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetNdivisions(507, kTRUE);
  h->GetYaxis()->SetNdivisions(507, kTRUE);
}


void sign_window(TVirtualPad *p, TH2D* h, TString xaxis, TString yaxis, TString title, TString title_size)
{
  Double_t height = 1-p->GetTopMargin()-p->GetBottomMargin();
  Double_t width = 1-p->GetLeftMargin()-p->GetRightMargin();
  //  cout << "pad width is: " << width << endl;
  //  cout << "pad height is: " << height << endl;
  if((width == 0.) || (height == 0.)) {
    cout << "In sign_window(). Width or height of pad is 0. Cannot set text size." << endl;
    exit(-1);
  }
  Double_t fontsizex = 0.04/width; //0.06
  Double_t fontsizey = 0.04/height;

  if(title_size == "small")
    {
      fontsizex *= 0.7;
      fontsizey *= 0.7;
    }
  if(title_size == "middle")
    {
      fontsizex *= 1.4;
      fontsizey *= 1.4;
    }
  if(title_size == "large")
    {
      fontsizex *= 3.5;
      fontsizey *= 3.5;
    }

  h->SetStats(kFALSE);
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xaxis);
  h->GetYaxis()->SetTitle(yaxis);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(fontsizex);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(fontsizey);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelSize(fontsizex);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetLabelSize(fontsizey);
  h->GetYaxis()->SetLabelFont(42);
  if(p->GetLogy())
    h->GetYaxis()->SetTitleOffset(2.5);
  else
    h->GetYaxis()->SetTitleOffset(1.1);
  h->GetXaxis()->SetNdivisions(507, kTRUE);
  h->GetYaxis()->SetNdivisions(507, kTRUE);
  if(title_size == "large")
    {
      h->GetXaxis()->SetTitleOffset(0.8);
      h->GetYaxis()->SetTitleOffset(0.5);
      h->GetYaxis()->SetNdivisions(406, kTRUE);
      h->GetXaxis()->SetTickLength(0.07);
    }
}


void make_clean_canv(TCanvas *c)
{
  c->SetBorderMode(0);
  c->SetFillColor(kWhite);
  c->SetFillStyle(0);

}
void make_clean_pads(TCanvas *c, Int_t n, Int_t gridx, Int_t gridy)
{

  for(Int_t i=0; i<n; i++)
    {
      c->GetPad(i+1)->SetBottomMargin(0.18);
      c->GetPad(i+1)->SetTopMargin(0.05);
      c->GetPad(i+1)->SetLeftMargin(0.15);
      c->GetPad(i+1)->SetRightMargin(0.05);
      c->GetPad(i+1)->SetFrameBorderMode(0);
      c->GetPad(i+1)->SetBorderMode(0);
      c->GetPad(i+1)->SetFillColor(kWhite);
      c->GetPad(i+1)->SetGrid(gridx, gridy);
    }
}

void make_clean_pads(TVirtualPad *c, Int_t n, Int_t gridx, Int_t gridy)
{

  for(Int_t i=0; i<n; i++)
    {
      //          c->GetPad(i+1)->SetBottomMargin(0.17);
      //            c->GetPad(i+1)->SetTopMargin(0.05);
      //            c->GetPad(i+1)->SetLeftMargin(0.17);
      //            c->GetPad(i+1)->SetRightMargin(0.05);
      c->GetPad(i+1)->SetTicks(1,1);
      c->GetPad(i+1)->SetFrameBorderMode(0);
      c->GetPad(i+1)->SetBorderMode(0);
      c->GetPad(i+1)->SetFillColor(kWhite);
      c->GetPad(i+1)->SetGrid(gridx, gridy);
    }
}

Double_t hist_maximum(TH1D** h, Int_t n)
{
  vector<Double_t> maximums;
  for(Int_t i=0; i<n; i++)
    {
      Double_t max = h[i]->GetMaximum();
      maximums.push_back(max);
      cout << i << " added " << max << endl;
    }
  sort(maximums.begin(), maximums.end());
  cout << "max is " << maximums.back() << endl;
  return maximums.back();
}

Double_t hist_maximum(TH1D** h, Int_t start, Int_t end)
{
  vector<Double_t> maximums;
  for(Int_t i=start-1; i<end; i++)
    {
      Double_t max = h[i]->GetMaximum();
      maximums.push_back(max);
      cout << i << " added " << max << endl;
    }
  sort(maximums.begin(), maximums.end());
  cout << "max is " << maximums.back() << endl;
  return maximums.back();
}

void prepare_canvas_for_cross_sec(TCanvas *c)
{
  c->Divide(1, 2, 0., 0.);
  //  c->GetPad(1)->SetBottomMargin(0.18);
  //  c->GetPad(1)->SetTopMargin(0.05);
  c->GetPad(1)->SetLeftMargin(0.13);
  c->GetPad(1)->SetRightMargin(0.05);
  c->GetPad(1)->SetFrameBorderMode(0);
  c->GetPad(1)->SetBorderMode(0);
  c->GetPad(1)->SetFillColor(kWhite);
  c->GetPad(1)->SetTicky();  
  c->GetPad(1)->SetTickx();  
  c->GetPad(1)->SetPad(0.1, 0.4, 0.9, 0.9);
  c->GetPad(1)->SetFrameBorderMode(0);
  c->GetPad(1)->SetBorderMode(0);
  //  c->GetPad(1)->SetTopMargin(0.25);
  //  c->GetPad(1)->cd();

  //  c->GetPad(2)->SetBottomMargin(0.18);
  //  c->GetPad(2)->SetTopMargin(0.05);
  c->GetPad(2)->SetLeftMargin(0.13);
  c->GetPad(2)->SetRightMargin(0.05);
  c->GetPad(2)->SetFrameBorderMode(0);
  c->GetPad(2)->SetBorderMode(0);
  c->GetPad(2)->SetFillColor(kWhite);
  c->GetPad(2)->SetTicky();  
  c->GetPad(2)->SetTickx();  
  c->GetPad(2)->SetPad(0.1, 0.2, 0.9, 0.4);
  c->GetPad(2)->SetFrameBorderMode(0);
  c->GetPad(2)->SetBorderMode(0);
  c->GetPad(2)->SetBottomMargin(0.25);
  
}

void prepare_canvas_for_cross_sec(TVirtualPad *c)
{
  c->Divide(1, 2, 0., 0.);
  //  c->GetPad(1)->SetBottomMargin(0.18);
  //  c->GetPad(1)->SetTopMargin(0.05);
  c->GetPad(1)->SetLeftMargin(0.15);
  c->GetPad(1)->SetRightMargin(0.05);
  c->GetPad(1)->SetFrameBorderMode(0);
  c->GetPad(1)->SetBorderMode(0);
  c->GetPad(1)->SetFillColor(kWhite);
  c->GetPad(1)->SetTicky();  
  c->GetPad(1)->SetTickx();  
  c->GetPad(1)->SetPad(0., 0.4, 0.9, 1.);
  c->GetPad(1)->SetFrameBorderMode(0);
  c->GetPad(1)->SetBorderMode(0);
  //  c->GetPad(1)->SetTopMargin(0.25);
  //  c->GetPad(1)->cd();

  //  c->GetPad(2)->SetBottomMargin(0.18);
  //  c->GetPad(2)->SetTopMargin(0.05);
  c->GetPad(2)->SetLeftMargin(0.15);
  c->GetPad(2)->SetRightMargin(0.05);
  c->GetPad(2)->SetFrameBorderMode(0);
  c->GetPad(2)->SetBorderMode(0);
  c->GetPad(2)->SetFillColor(kWhite);
  c->GetPad(2)->SetTicky();  
  c->GetPad(2)->SetTickx();  
  c->GetPad(2)->SetPad(0., 0., 0.9, 0.4);
  c->GetPad(2)->SetFrameBorderMode(0);
  c->GetPad(2)->SetBorderMode(0);
  c->GetPad(2)->SetBottomMargin(0.35);
  
}
