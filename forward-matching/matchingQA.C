#include <MFTTracking/Constants.h>

TCanvas c("c", "c", 1200, 800);
c.SetGridx(true);
c.SetGridy(true);

constexpr int nPoints = 3;

TH1* GetTH1(TFile* f, TString histname)
{
  return (TH1*)f->Get(histname);

  // TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey* key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " << key << std::endl;
  if (!key)
    return NULL;
  return (TH1*)key->ReadObjectAny(TH1::Class());
}

TH2* GetTH2(TFile* f, TString histname)
{
  return (TH2*)f->Get(histname);

  // TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey* key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " << key << std::endl;
  if (!key)
    return NULL;
  return (TH2*)key->ReadObjectAny(TH2::Class());
}

TH1F* GetGoodRankingFraction(TH2* rankingHist)
{
  TH1F* hist = new TH1F((std::string(rankingHist->GetName()) + "-good-ranking-fraction").c_str(),
                        rankingHist->GetTitle(),
                        rankingHist->GetXaxis()->GetNbins(),
                        rankingHist->GetXaxis()->GetXmin(),
                        rankingHist->GetXaxis()->GetXmax());
  hist->SetTitle("Fraction of correctly ranked matches");
  hist->GetXaxis()->SetTitle(rankingHist->GetXaxis()->GetTitle());

  for (int i = 1; i <= rankingHist->GetXaxis()->GetNbins(); i++) {
    float total = 0;
    for (int j = 2; j <= rankingHist->GetYaxis()->GetNbins(); j++) {
      total += rankingHist->GetBinContent(i, j);
    }

    // ratio between matches with correct ranking and all ranked matches
    float fraction = (total > 0) ? rankingHist->GetBinContent(i, 2) / total : 0.f;
    hist->SetBinContent(i, fraction);
  }

  hist->SetMinimum(0.5);

  return hist;
}

void MatchRankingCompareMatchingMethods(TFile* rootFile,
                                        const std::vector<std::string>& matchingMethods,
                                        std::string histName)
{
  TH1* h1;
  TH2* h2;
  TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);

  // same particle in MFT and MCH
  legend->Clear();
  c.Clear();
  c.SetLogy(kTRUE);
  int index = 0;
  for (auto method : matchingMethods) {
    std::string path = std::string("qa-matching/matching/MC/") + method + "/";
    h1 = GetTH1(rootFile, (path + histName).c_str());
    if (!h1) {
      std::cout << "Histogram \"" << (path + histName) << "\" not found" << std::endl;
      continue;
    }
    h1->SetLineColor(index + 1);
    // Normalize to unit area
    h1->Scale(1.0 / h1->Integral());
    if (index == 0) {
      h1->Draw("HIST");
      h1->SetMaximum(1.0);
    } else
      h1->Draw("same");
    legend->AddEntry(h1, method.c_str(), "l");
    index += 1;
  }
  legend->Draw();
  c.SaveAs("matchingQA.pdf");

  legend->SetY1NDC(0.2);
  legend->SetY2NDC(0.5);
  std::vector<std::string> variables{
    "P", "Pt"};
  for (auto variable : variables) {
    legend->Clear();
    c.Clear();
    c.SetLogy(kFALSE);
    int index = 0;
    for (auto method : matchingMethods) {
      std::string path = std::string("qa-matching/matching/MC/") + method + "/";
      h2 = GetTH2(rootFile, (path + histName + "Vs" + variable).c_str());
      h1 = GetGoodRankingFraction(h2);
      h1->SetLineColor(index + 1);
      if (index == 0)
        h1->Draw();
      else
        h1->Draw("same");
      legend->AddEntry(h1, method.c_str(), "l");
      index += 1;
    }
    legend->Draw();
    c.SaveAs("matchingQA.pdf");
  }
}

void MatchingScoreCompareMatchingMethods(TFile* rootFile,
                                         const std::vector<std::string>& matchingMethods,
                                         std::string histName)
{
  TH1* h1;
  TH2* h2;
  TLegend* legend = new TLegend(0.15, 0.8, 0.85, 0.9);
  legend->SetNColumns(3);

  legend->Clear();
  c.Clear();
  c.SetLogy(kFALSE);
  float max = 0;
  // matching score distributions
  std::vector<TH1*> h1vec;
  // cumulative distributions
  std::vector<TH1*> ch1vec;
  for (auto method : matchingMethods) {
    std::string path = std::string("qa-matching/matching/MC/") + method + "/";
    h1 = GetTH1(rootFile, (path + histName).c_str());
    if (h1->GetMaximum() > max) {
      max = h1->GetMaximum();
    }
    h1vec.push_back(h1);

    TH1F* ch1 = (TH1F*)h1->Clone();
    float integral = h1->Integral();
    int binMax = h1->GetXaxis()->GetNbins();
    for (int bin = 1; bin <= h1->GetXaxis()->GetNbins(); bin++) {
      float cumulant = h1->Integral(bin, binMax);
      ch1->SetBinContent(bin, cumulant / integral);
    }
    ch1vec.push_back(ch1);
  }

  for (int i = 0; i < h1vec.size(); i++) {
    auto* hist = h1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(1.2f * max);
      hist->Draw();
    } else {
      hist->Draw("same");
    }
    legend->AddEntry(hist, matchingMethods[i].c_str(), "l");
  }
  legend->Draw();
  c.SaveAs("matchingQA.pdf");

  for (int i = 0; i < ch1vec.size(); i++) {
    auto* hist = ch1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(1.2f);
      hist->Draw();
    } else {
      hist->Draw("same");
    }
    // legend->AddEntry(hist, matchingMethods[i].c_str(), "l");
  }
  c.SaveAs("matchingQA.pdf");
}

// Identical to the above, except splits into variable ranges for each matching method
// Each matching method will have its own histogram to compare the score distributions vs variable ranges
// at some point make a generic function that can plot X vs Y while varying Z
void MatchingScoreCompareVsVariable(TFile* rootFile,
                                    std::string matchingMethod,
                                    std::string histName,
                                    std::string variable,
                                    const std::vector<std::pair<double, double>> ranges)
{
  TH1* h1;
  TH2* h2;
  TLegend* legend = new TLegend(0.15, 0.8, 0.85, 0.9);
  legend->SetNColumns(3);

  legend->Clear();
  c.Clear();
  c.SetLogy(kFALSE);
  float max = 0;
  // matching score distributions
  std::vector<TH1*> h1vec;
  // cumulative distributions
  std::vector<TH1*> ch1vec;
  std::string path = std::string("qa-matching/matching/MC/") + matchingMethod + "/";
  h2 = GetTH2(rootFile, (path + histName + "Vs" + variable).c_str());

  for (const auto& range : ranges) {
    h1 = h2->ProjectionY(
      TString::Format("%g<%s<%g", range.first,variable.c_str(), range.second).Data(),
      h2->GetXaxis()->FindBin(range.first),
      h2->GetXaxis()->FindBin(range.second) - 1);

    h1->Scale(1.0 / h1->Integral()); // noramalize to 1
    if (h1->GetMaximum() > max) {
      max = h1->GetMaximum();
    }
    h1vec.push_back(h1);

    TH1F* ch1 = (TH1F*)h1->Clone();
    float integral = h1->Integral();
    int binMax = h1->GetXaxis()->GetNbins();
    for (int bin = 1; bin <= h1->GetXaxis()->GetNbins(); bin++) {
      float cumulant = h1->Integral(bin, binMax);
      ch1->SetBinContent(bin, cumulant / integral);
    }
    ch1vec.push_back(ch1);
  }

  for (int i = 0; i < h1vec.size(); i++) {
    auto* hist = h1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(1.2f * max);
      hist->Draw("HIST");

    } else {
      hist->Draw("HIST SAME");
    }
    legend->AddEntry(hist, hist->GetName(), "l");
  }
  legend->Draw();
  c.SaveAs("matchingQA.pdf");

  legend->Clear();
  for (int i = 0; i < ch1vec.size(); i++) {
    auto* hist = ch1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(1.2f);
      hist->Draw("HIST");
    } else {
      hist->Draw("HIST SAME");
    }
    legend->AddEntry(hist, hist->GetName(), "l");
  }
  legend->Draw();
  c.SaveAs("matchingQA.pdf");
}

void PlotMatchRanking(TFile* rootFile)
{
  std::vector<std::string> matchingMethods{
    "Prod",
    "MatchXYPhiTanl",
    "MatchXYPhiTanlMom"};

  TH1* h1;
  TH2* h2;
  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6);
  TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  const std::vector<std::pair<double, double>> ranges = {
    {0.f, 10.f},
    {10.f, 15.f},
    {15.f, 20.f},
    {20.f, 30.f},
    {30.f, 50.f},
    {50.f, 100.f}};

  // true match ranking for good MCH tracks
  c.Clear();
  title->Clear();
  title->AddText("True match ranking");
  title->AddText("Good MCH tracks");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchRankingCompareMatchingMethods(rootFile, matchingMethods, "trueMatchRankingGoodMCH");

  c.Clear();
  title->Clear();
  title->AddText("True match ranking");
  title->AddText("Good MCH tracks, paired with MFT");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchRankingCompareMatchingMethods(rootFile, matchingMethods, "trueMatchRankingGoodPairedMCH");

  c.Clear();
  title->Clear();
  title->AddText("Matching score distribution");
  title->AddText("Good MCH tracks, all particles");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchingScoreCompareMatchingMethods(rootFile, matchingMethods, "trueMatchScore");

  c.Clear();
  title->Clear();
  title->AddText("Matching score distribution");
  title->AddText("Good MCH tracks, fake matches");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchingScoreCompareMatchingMethods(rootFile, matchingMethods, "fakeMatchScore");

  // beginning of my additions
  // opt for per method comparison of p ranges for simplicity, can of course swap toper p range comparison of methods if desired

  for (std::string entry : matchingMethods) {
    c.Clear();
    title->Clear();
    title->AddText(("Matching score distribution, for " + entry).c_str());
    title->AddText(("Good MCH tracks, all particles, for " + entry).c_str());
    title->Draw();
    c.SaveAs("matchingQA.pdf");
    MatchingScoreCompareVsVariable(rootFile, entry, "trueMatchScore", "P", ranges);
  }

  for (std::string entry : matchingMethods) {
    c.Clear();
    title->Clear();
    title->AddText(("Matching score distribution, P dependence, for " + entry).c_str());
    title->AddText(("Good MCH tracks, fake matches, P dependence for " + entry).c_str());
    title->Draw();
    c.SaveAs("matchingQA.pdf");
    MatchingScoreCompareVsVariable(rootFile, entry, "fakeMatchScore", "P", ranges);
    
  }
}

void matchingQA()
{
  TFile* fAnalysisResults;

  fAnalysisResults = new TFile("O-O_DQ_LHC25i4.root");

  gStyle->SetOptStat(0);
  // gStyle->SetOptStat(1111);
  // gStyle->SetOptFit(1111);

  c.SaveAs("matchingQA.pdf(");

  PlotMatchRanking(fAnalysisResults);

  c.Clear();
  c.SaveAs("matchingQA.pdf)");
}
