{
//=========Macro generated from canvas: c1/c1
//=========  (Thu May  5 12:13:40 2011) by ROOT version5.26/00c
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
   c1_1->Range(0.7581625,-3.377635,3.079989,9.831448);
   c1_1->SetFillColor(0);
   c1_1->SetBorderMode(0);
   c1_1->SetBorderSize(2);
   c1_1->SetLogx();
   c1_1->SetLogy();
   c1_1->SetLeftMargin(0.18);
   c1_1->SetRightMargin(0.13);
   c1_1->SetTopMargin(0.08);
   c1_1->SetBottomMargin(0.18);
   c1_1->SetFrameBorderMode(0);
   c1_1->SetFrameBorderMode(0);
   
   TH2D *h_window0 = new TH2D("h_window0","",10,15,600,10,0.1,5.952804e+08);
   h_window0->SetStats(0);
   h_window0->GetXaxis()->SetTitle("Idfmckin");
   h_window0->GetXaxis()->CenterTitle(true);
   h_window0->GetXaxis()->SetNdivisions(507);
   h_window0->GetXaxis()->SetLabelFont(42);
   h_window0->GetXaxis()->SetLabelSize(0.07536232);
   h_window0->GetXaxis()->SetTitleSize(0.07536232);
   h_window0->GetXaxis()->SetTitleFont(42);
   h_window0->GetYaxis()->SetTitle("Particles");
   h_window0->GetYaxis()->CenterTitle(true);
   h_window0->GetYaxis()->SetNdivisions(507);
   h_window0->GetYaxis()->SetLabelFont(42);
   h_window0->GetYaxis()->SetLabelSize(0.07027027);
   h_window0->GetYaxis()->SetTitleSize(0.07027027);
   h_window0->GetYaxis()->SetTitleOffset(1.2);
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
   
   TH1I *h_had_id2 = new TH1I("h_had_id2","Fill id for every hadron level object",2335,0,2335);
   h_had_id2->SetBinContent(18,77329);
   h_had_id2->SetBinContent(19,85556);
   h_had_id2->SetBinContent(20,77298);
   h_had_id2->SetBinContent(21,85623);
   h_had_id2->SetBinContent(22,2539);
   h_had_id2->SetBinContent(23,2539);
   h_had_id2->SetBinContent(24,436404);
   h_had_id2->SetBinContent(25,428176);
   h_had_id2->SetBinContent(26,85899);
   h_had_id2->SetBinContent(27,77573);
   h_had_id2->SetBinContent(30,5.952804e+07);
   h_had_id2->SetBinContent(55,2.507699e+07);
   h_had_id2->SetBinContent(56,2.353435e+07);
   h_had_id2->SetBinContent(58,8105);
   h_had_id2->SetBinContent(59,3292365);
   h_had_id2->SetBinContent(60,2978747);
   h_had_id2->SetBinContent(63,2962676);
   h_had_id2->SetBinContent(64,2959730);
   h_had_id2->SetBinContent(191,4009292);
   h_had_id2->SetBinContent(192,1344564);
   h_had_id2->SetBinContent(193,2625222);
   h_had_id2->SetBinContent(194,1250720);
   h_had_id2->SetBinContent(195,846090);
   h_had_id2->SetBinContent(196,427600);
   h_had_id2->SetBinContent(197,170713);
   h_had_id2->SetBinContent(198,102258);
   h_had_id2->SetBinContent(201,106059);
   h_had_id2->SetBinContent(202,89750);
   h_had_id2->SetBinContent(203,40693);
   h_had_id2->SetBinContent(204,33881);
   h_had_id2->SetBinContent(205,36678);
   h_had_id2->SetBinContent(206,32981);
   h_had_id2->SetBinContent(557,747);
   h_had_id2->SetBinContent(558,801);
   h_had_id2->SetBinError(18,88.15973);
   h_had_id2->SetBinError(19,92.731);
   h_had_id2->SetBinError(20,88.14191);
   h_had_id2->SetBinError(21,92.76744);
   h_had_id2->SetBinError(22,15.97745);
   h_had_id2->SetBinError(23,15.97745);
   h_had_id2->SetBinError(24,209.4318);
   h_had_id2->SetBinError(25,207.4482);
   h_had_id2->SetBinError(26,92.91672);
   h_had_id2->SetBinError(27,88.29901);
   h_had_id2->SetBinError(30,2446.014);
   h_had_id2->SetBinError(55,1587.581);
   h_had_id2->SetBinError(56,1537.975);
   h_had_id2->SetBinError(58,28.54274);
   h_had_id2->SetBinError(59,575.244);
   h_had_id2->SetBinError(60,547.1608);
   h_had_id2->SetBinError(63,545.6827);
   h_had_id2->SetBinError(64,545.4113);
   h_had_id2->SetBinError(191,634.7928);
   h_had_id2->SetBinError(192,367.6114);
   h_had_id2->SetBinError(193,513.6664);
   h_had_id2->SetBinError(194,354.5506);
   h_had_id2->SetBinError(195,291.6127);
   h_had_id2->SetBinError(196,207.3085);
   h_had_id2->SetBinError(197,130.9881);
   h_had_id2->SetBinError(198,101.3792);
   h_had_id2->SetBinError(201,103.2459);
   h_had_id2->SetBinError(202,94.97663);
   h_had_id2->SetBinError(203,63.95332);
   h_had_id2->SetBinError(204,58.35552);
   h_had_id2->SetBinError(205,60.71577);
   h_had_id2->SetBinError(206,57.5747);
   h_had_id2->SetBinError(557,8.668689);
   h_had_id2->SetBinError(558,8.977809);
   h_had_id2->SetEntries(1.321483e+09);
   h_had_id2->SetMarkerStyle(20);
   h_had_id2->SetMarkerSize(0.5);
   h_had_id2->GetXaxis()->SetMoreLogLabels();
   h_had_id2->Draw("E P X0 SAME");
   
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
   entry=leg->AddEntry("h_had_id2","Common Ntuples","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(0.5);
   leg->Draw();
   c1_1->Modified();
   c1->cd();
  
// ------------>Primitives in pad: c1_2
   c1_2 = new TPad("c1_2", "c1_2",0.01,0.01,0.99,0.49);
   c1_2->Draw();
   c1_2->cd();
   c1_2->Range(-13.04348,-3.195352,59.42029,9.001046);
   c1_2->SetFillColor(0);
   c1_2->SetBorderMode(0);
   c1_2->SetBorderSize(2);
   c1_2->SetLogy();
   c1_2->SetLeftMargin(0.18);
   c1_2->SetRightMargin(0.13);
   c1_2->SetTopMargin(0.08);
   c1_2->SetBottomMargin(0.18);
   c1_2->SetFrameBorderMode(0);
   c1_2->SetFrameBorderMode(0);
   
   TH2D *h_window1 = new TH2D("h_window1","",10,0,50,10,0.1,1.06007e+08);
   h_window1->SetStats(0);
   h_window1->GetXaxis()->SetTitle("idpart");
   h_window1->GetXaxis()->CenterTitle(true);
   h_window1->GetXaxis()->SetNdivisions(507);
   h_window1->GetXaxis()->SetLabelFont(42);
   h_window1->GetXaxis()->SetLabelSize(0.07536232);
   h_window1->GetXaxis()->SetTitleSize(0.07536232);
   h_window1->GetXaxis()->SetTitleFont(42);
   h_window1->GetYaxis()->SetTitle("Particles");
   h_window1->GetYaxis()->CenterTitle(true);
   h_window1->GetYaxis()->SetNdivisions(507);
   h_window1->GetYaxis()->SetLabelFont(42);
   h_window1->GetYaxis()->SetLabelSize(0.07027027);
   h_window1->GetYaxis()->SetTitleSize(0.07027027);
   h_window1->GetYaxis()->SetTitleOffset(1.2);
   h_window1->GetYaxis()->SetTitleFont(42);
   h_window1->Draw("");
   
   TH1I *h_part_id1 = new TH1I("h_part_id1","Fill part id for every parton level object",50,0,50);
   h_part_id1->SetBinContent(2,1220935);
   h_part_id1->SetBinContent(3,477149);
   h_part_id1->SetBinContent(4,4192238);
   h_part_id1->SetBinContent(5,1625689);
   h_part_id1->SetBinContent(6,375261);
   h_part_id1->SetBinContent(7,373192);
   h_part_id1->SetBinContent(8,718038);
   h_part_id1->SetBinContent(9,722302);
   h_part_id1->SetBinContent(10,8584);
   h_part_id1->SetBinContent(11,8593);
   h_part_id1->SetBinContent(34,1.060099e+07);
   h_part_id1->SetBinContent(41,1911383);
   h_part_id1->SetBinContent(43,530007);
   h_part_id1->SetBinContent(45,866741);
   h_part_id1->SetEntries(2.36311e+07);

   ci = TColor::GetColor("#000000");
   h_part_id1->SetLineColor(ci);
   h_part_id1->Draw("E X0 HIST SAME");
   
   TH1I *h_part_id2 = new TH1I("h_part_id2","Fill part id for every parton level object",50,0,50);
   h_part_id2->SetBinContent(2,1221767);
   h_part_id2->SetBinContent(3,478038);
   h_part_id2->SetBinContent(4,4190621);
   h_part_id2->SetBinContent(5,1624303);
   h_part_id2->SetBinContent(6,375879);
   h_part_id2->SetBinContent(7,372949);
   h_part_id2->SetBinContent(8,718321);
   h_part_id2->SetBinContent(9,723084);
   h_part_id2->SetBinContent(10,8609);
   h_part_id2->SetBinContent(11,8656);
   h_part_id2->SetBinContent(34,1.06007e+07);
   h_part_id2->SetBinContent(41,1911687);
   h_part_id2->SetBinContent(43,529115);
   h_part_id2->SetBinContent(45,867363);
   h_part_id2->SetBinError(2,350.3997);
   h_part_id2->SetBinError(3,219.1801);
   h_part_id2->SetBinError(4,648.946);
   h_part_id2->SetBinError(5,404.0201);
   h_part_id2->SetBinError(6,194.3541);
   h_part_id2->SetBinError(7,193.595);
   h_part_id2->SetBinError(8,268.676);
   h_part_id2->SetBinError(9,269.5652);
   h_part_id2->SetBinError(10,29.41393);
   h_part_id2->SetBinError(11,29.49502);
   h_part_id2->SetBinError(34,1032.135);
   h_part_id2->SetBinError(41,438.3064);
   h_part_id2->SetBinError(43,230.5921);
   h_part_id2->SetBinError(45,295.2365);
   h_part_id2->SetEntries(2.351502e+08);
   h_part_id2->SetMarkerStyle(20);
   h_part_id2->SetMarkerSize(0.5);
   h_part_id2->Draw("E P X0 SAME");
   c1_2->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
