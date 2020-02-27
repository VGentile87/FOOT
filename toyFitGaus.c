/////
//
//
void toyFitGaus(){

  
  TFile * file1 = TFile::Open("histoVR01.root");
  TH1F * h1 = (TH1F*) file1->Get("hVR012");  
  
  TFile * file2 = TFile::Open("histoV01_inV012.root");
  TH1F * h2 = (TH1F*) file2->Get("hv012");

  TFile * file3 = TFile::Open("histoVR012.root");
  TH1F * h3 = (TH1F*) file3->Get("h");

  TFile * file4 = TFile::Open("Profile_2_1.root");
  TH1F * h4 = (TH1F*) file4->Get("profileX");
  

  TH1F *hdiv = new TH1F("hdiv","hdiv",50,-3,7);
  hdiv = (TH1F*)h2->Clone();
  hdiv->Divide(h1);
  
  float pesi[50]={};
  float shift[50]={};
  for(int i=0;i<50;i++){
    pesi[i]=hdiv->GetBinContent(i+1);
    if(pesi[i]==0)pesi[i]=1;
    pesi[0]=0.12;
    shift[i]=h4->GetBinContent(i+1);
    //cout << i << " " << pesi[i] << endl;
  }
  
  TF1 * g1 = new TF1("g1","gaus",-3,7);
  TF1 * g2 = new TF1("g2","gaus",-3,7);
  TF1 * g3 = new TF1("g3","gaus",-3,7);
  TF1 * gsum = new TF1("gsum","g1+g2+g3",-3,7);
  
  /*
  gsum->SetParameter(0,83.7);
  gsum->SetParameter(1,-1.1);
  gsum->SetParameter(2,1.012);
  gsum->SetParameter(3,12.44);
  gsum->SetParameter(4,1.778);
  gsum->SetParameter(5,1.049);
  gsum->SetParameter(6,5.478);
  gsum->SetParameter(7,4.624);
  gsum->SetParameter(8,0.619);
*/
  
  /*gsum->SetParameter(0,164.9);
  gsum->SetParameter(1,-1.1);
  gsum->SetParameter(2,1.084);
  gsum->SetParameter(3,24.04);
  gsum->SetParameter(4,1.913);
  gsum->SetParameter(5,0.8211);
  gsum->SetParameter(6,11.63);
  gsum->SetParameter(7,4.499);
  gsum->SetParameter(8,0.8993);
  */
  /*VR01/
  gsum->SetParameter(0,209.6);
  gsum->SetParameter(1,-0.4142);
  gsum->SetParameter(2,0.57);
  gsum->SetParameter(3,111.8);
  gsum->SetParameter(4,0.38);
  gsum->SetParameter(5,0.99);
  gsum->SetParameter(6,24.85);
  gsum->SetParameter(7,0.19);
  gsum->SetParameter(8,1.64);
  */

  // VR012
  gsum->SetParameter(0,255.4);
  gsum->SetParameter(1,-0.82);
  gsum->SetParameter(2,1.0);
  gsum->SetParameter(3,32.26);
  gsum->SetParameter(4,1.90);
  gsum->SetParameter(5,0.65);
  gsum->SetParameter(6,18.2);
  gsum->SetParameter(7,3.76);
  gsum->SetParameter(8,0.67);

  TRandom3 *rnd = new TRandom3();

  TH1F * heff = new TH1F("heff","heff",50,-3,7);
  TH1F * hnoeff = new TH1F("hnoeff","hnoeff",50,-3,7);
  
  for(int i=0;i<100000;i++){
    float val_f = gsum->GetRandom(-3,7);
    
    
    float start=-3;
    float step=0.2;
    float delta=0;

    for(int j=0;j<50;j++){
      if(val_f>start && val_f<(start+step))delta=shift[j];
      start+=step;
    }

    float corr = 0.54*val_f + 0.17;
    //float val_i = val_f -  rnd->Gaus(delta,0.66);
    float val_i = rnd->Gaus(corr,0.66);
    start=-3;
    for(int j=0;j<50;j++){
      //cout << pesi[j] << endl;
      if(val_i>start && val_i<(start+step)){
	heff->Fill(val_i,(pesi[j]));
	hnoeff->Fill(val_i);
      }
	start+=step;
    }
    //val_i -= rnd->Gaus(delta,0.66);
  }
  

  cout << heff->Integral() << endl;
  heff->Scale(h2->Integral()/heff->Integral());
  hnoeff->Scale(h1->Integral()/hnoeff->Integral());

  cout <<"Kolmogorov test no eff " << hnoeff->KolmogorovTest(h1) << endl;
  cout <<"Kolmogorov test eff " << heff->KolmogorovTest(h2) << endl;
  
  TCanvas *c1 = new TCanvas("c1","c1",1500,1000);
  c1->Divide(3,2);
  c1->cd(1);
  h1->Draw();
  h2->Draw("sames");
  c1->cd(2);
  hdiv->Draw();
  c1->cd(3);
  h4->Draw();
  h4->Fit("pol1","R","",-1.8,3);
  c1->cd(4);
  h3->Draw();
  gsum->Draw("sames");
  c1->cd(5);
  h1->Draw();
  h1->SetLineColor(kRed);
  hnoeff->Draw("sames");
  c1->cd(6);
  h2->Draw();
  h2->SetLineColor(kGreen);
  heff->Draw("sames");
  
}
