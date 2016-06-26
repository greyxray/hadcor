{
//=========Macro generated from canvas: c1/c1
//=========  (Wed May 18 01:33:57 2011) by ROOT version5.26/00c
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,1400,600);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: c1_1
   TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.51,0.99,0.99);
   c1_1->Draw();
   c1_1->cd();
   c1_1->Range(0.9221061,-3.28502,2.875838,9.409534);
   c1_1->SetFillColor(0);
   c1_1->SetBorderMode(0);
   c1_1->SetBorderSize(2);
   c1_1->SetLogx();
   c1_1->SetLogy();
   c1_1->SetLeftMargin(0.13);
   c1_1->SetRightMargin(0.05);
   c1_1->SetTopMargin(0.05);
   c1_1->SetBottomMargin(0.18);
   c1_1->SetFrameBorderMode(0);
   c1_1->SetFrameBorderMode(0);
   
   TH2D *h_window0 = new TH2D("h_window0","",10,15,600,10,0.1,5.953969e+08);
   h_window0->SetStats(0);
   h_window0->GetXaxis()->SetTitle("Idfmckin");
   h_window0->GetXaxis()->CenterTitle(true);
   h_window0->GetXaxis()->SetNdivisions(507);
   h_window0->GetXaxis()->SetLabelFont(42);
   h_window0->GetXaxis()->SetLabelSize(0.05073171);
   h_window0->GetXaxis()->SetTitleSize(0.06341463);
   h_window0->GetXaxis()->SetTitleFont(42);
   h_window0->GetYaxis()->SetTitle("Particles");
   h_window0->GetYaxis()->CenterTitle(true);
   h_window0->GetYaxis()->SetNdivisions(507);
   h_window0->GetYaxis()->SetLabelFont(42);
   h_window0->GetYaxis()->SetLabelSize(0.05402597);
   h_window0->GetYaxis()->SetTitleSize(0.06753247);
   h_window0->GetYaxis()->SetTitleFont(42);
   h_window0->Draw("");
   
   TH1I *h_had_id1 = new TH1I("h_had_id1","Fill id for every hadron level object",2335,0,2335);
   h_had_id1->SetBinContent(18,77361);
   h_had_id1->SetBinContent(19,85314);
   h_had_id1->SetBinContent(20,77510);
   h_had_id1->SetBinContent(21,85492);
   h_had_id1->SetBinContent(22,2613);
   h_had_id1->SetBinContent(23,2613);
   h_had_id1->SetBinContent(24,435902);
   h_had_id1->SetBinContent(25,427949);
   h_had_id1->SetBinContent(26,85778);
   h_had_id1->SetBinContent(27,77796);
   h_had_id1->SetBinContent(30,5.953969e+07);
   h_had_id1->SetBinContent(55,2.507536e+07);
   h_had_id1->SetBinContent(56,2.353263e+07);
   h_had_id1->SetBinContent(58,8119);
   h_had_id1->SetBinContent(59,3289770);
   h_had_id1->SetBinContent(60,2978194);
   h_had_id1->SetBinContent(63,2958762);
   h_had_id1->SetBinContent(64,2957788);
   h_had_id1->SetBinContent(191,4011466);
   h_had_id1->SetBinContent(192,1345140);
   h_had_id1->SetBinContent(193,2623900);
   h_had_id1->SetBinContent(194,1250681);
   h_had_id1->SetBinContent(195,846145);
   h_had_id1->SetBinContent(196,427778);
   h_had_id1->SetBinContent(197,170753);
   h_had_id1->SetBinContent(198,102344);
   h_had_id1->SetBinContent(201,105755);
   h_had_id1->SetBinContent(202,89455);
   h_had_id1->SetBinContent(203,40617);
   h_had_id1->SetBinContent(204,33979);
   h_had_id1->SetBinContent(205,36840);
   h_had_id1->SetBinContent(206,33006);
   h_had_id1->SetBinContent(557,725);
   h_had_id1->SetBinContent(558,787);
   h_had_id1->SetEntries(1.32818e+08);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000000");
   h_had_id1->SetLineColor(ci);
   h_had_id1->GetXaxis()->SetMoreLogLabels();
   h_had_id1->Draw("E X0 HIST SAME");
   
   TLegend *leg = new TLegend(0.7,0.7,0.85,0.89,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("h_had_id1","Orange Ntuples","l");

   ci = TColor::GetColor("#000000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   c1_1->Modified();
   c1->cd();
  
// ------------>Primitives in pad: c1_2
   c1_2 = new TPad("c1_2", "c1_2",0.01,0.01,0.99,0.49);
   c1_2->Draw();
   c1_2->cd();
   c1_2->Range(-7.926829,-3.285102,53.04878,9.409907);
   c1_2->SetFillColor(0);
   c1_2->SetBorderMode(0);
   c1_2->SetBorderSize(2);
   c1_2->SetLogy();
   c1_2->SetLeftMargin(0.13);
   c1_2->SetRightMargin(0.05);
   c1_2->SetTopMargin(0.05);
   c1_2->SetBottomMargin(0.18);
   c1_2->SetFrameBorderMode(0);
   c1_2->SetFrameBorderMode(0);
   
   TH2D *h_window1 = new TH2D("h_window1","",10,0,50,10,0.1,5.958769e+08);
   h_window1->SetStats(0);
   h_window1->GetXaxis()->SetTitle("idpart");
   h_window1->GetXaxis()->CenterTitle(true);
   h_window1->GetXaxis()->SetNdivisions(507);
   h_window1->GetXaxis()->SetLabelFont(42);
   h_window1->GetXaxis()->SetLabelSize(0.05073171);
   h_window1->GetXaxis()->SetTitleSize(0.06341463);
   h_window1->GetXaxis()->SetTitleFont(42);
   h_window1->GetYaxis()->SetTitle("Particles");
   h_window1->GetYaxis()->CenterTitle(true);
   h_window1->GetYaxis()->SetNdivisions(507);
   h_window1->GetYaxis()->SetLabelFont(42);
   h_window1->GetYaxis()->SetLabelSize(0.05402597);
   h_window1->GetYaxis()->SetTitleSize(0.06753247);
   h_window1->GetYaxis()->SetTitleFont(42);
   h_window1->Draw("");
   
   TH1I *h_part_id1 = new TH1I("h_part_id1","Fill part id for every parton level object",50,0,50);
   h_part_id1->SetBinContent(1,2583539);
   h_part_id1->SetBinContent(2,8741648);
   h_part_id1->SetBinContent(3,2068836);
   h_part_id1->SetBinContent(4,2.268905e+07);
   h_part_id1->SetBinContent(5,4242053);
   h_part_id1->SetBinContent(6,1789950);
   h_part_id1->SetBinContent(7,2833490);
   h_part_id1->SetBinContent(8,2542024);
   h_part_id1->SetBinContent(9,2517813);
   h_part_id1->SetBinContent(10,36297);
   h_part_id1->SetBinContent(11,29956);
   h_part_id1->SetBinContent(34,5.958769e+07);
   h_part_id1->SetBinContent(41,1.606212e+07);
   h_part_id1->SetBinContent(43,3516132);
   h_part_id1->SetBinContent(45,7350512);
   h_part_id1->SetEntries(1.365911e+08);

   ci = TColor::GetColor("#000000");
   h_part_id1->SetLineColor(ci);
   h_part_id1->Draw("E X0 HIST SAME");
   c1_2->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
