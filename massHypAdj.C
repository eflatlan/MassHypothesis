#include <stdio.h>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH2F.h>
#include <Riostream.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TVector2.h>
#include <TH1D.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <memory>
#include "SaveData.cpp"
#include "RandomValues.cpp"

#include <Math/Vector3D.h>
#include <Math/Vector2D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "Math/Vector3D.h"
#include "Math/Vector2D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"



static constexpr double arrWaveLenDefault[30] = {
  162, 164, 166, 168, 170, 172, 174, 176, 178, 180,
  182, 184, 186, 188, 190, 192, 194, 196, 198, 200,
  202, 204, 206, 208, 210, 212, 214, 216, 218, 220};
static constexpr double nm2eV = 1239.842609;



class ParticleInfo : public TObject {
public:
    // Your data members here
    double momentum;
    double mass;
    double energy;
    float refractiveIndex;
    double ckov;
    TH2F* map;

    ParticleInfo() {
        // Initialize your data members
    }

    virtual ~ParticleInfo() {}

    ClassDef(ParticleInfo,1) // For ROOT dictionary
};


/*
struct ParticleInfo {
    double momentum;
    double mass;
    double energy;
    float refractiveIndex;
    double ckov;
    TH2F* map;
}; */ 


double arrW[750]= {0.};


using namespace o2;
using namespace o2::hmpid;

#include <HMPIDBase/Param.h>
void setStyleInd(TH2* th1f, double ratio = 1.2);
void setStyleInd(TH1* th1f, double ratio = 1.2);



const double mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in
std::array<double, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
// get Ring-radius from Cherenkov-angle
double getRadiusFromCkov(double ckovAngle);
double calcRingGeom(double ckovAng, int level);
float GetFreonIndexOfRefraction(float x);
float GetQuartzIndexOfRefraction(float x);
double BackgroundFunc(double *x, double *par);

const double fDTheta = 0.001;  // increment
double kThetaMax = 0.75;
int nChannels = (int)(kThetaMax / fDTheta + 0.5);


void setStyle();


auto phots = new TH1D("Photon Candidates", "Photon Candidates;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
auto photsw = new TH1D("Photon Weights", "Photon Weights;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
auto resultw = new TH1D("Sum of Weights in Window at Bin", "Sum of Weights in Window at bin;angle [rad]; counts/1 mrad", nChannels, 0, kThetaMax);
TH1F *hTheta = new TH1F("Background","Background; angle [rad]; counts/1 mrad",750,0.,0.75); 

TH1F *hThetaCh = new TH1F("Cherenkov Photons","Cherenkov Photons; angle [rad]; counts/1 mrad",750,0.,0.75); 

double meanCherenkovAngle;


std::vector<double> photonCandidates;

double ckovTrackOut = 0;
double /*std::array<TH1D*, 3>*/ houghResponse(std::vector<double>& photonCandidates, double fWindowWidth);



const double defaultPhotonEnergy = 6.75; 
const double refIndexFreon = GetFreonIndexOfRefraction(defaultPhotonEnergy);
const double refIndexQuartz = GetQuartzIndexOfRefraction(defaultPhotonEnergy);
const double  refIndexCH4 = 1.00; 


TH2F* backgroundStudy(double ckovActual = 0.5, double occupancy = 0.03);

const double CH4GapWidth = 8;
const double  RadiatorWidth = 1.;
const double  QuartzWindowWidth = 0.5;
const double  EmissionLenght = RadiatorWidth/2;


TH1* getMaxInRange(TH1* th1, double& up, double mid, double width);
double getMaxInRange(TH1* th1, int start, int width);
double getMaxInRange(TH1* th1, double mid, double width);


// mass_Pion_sq mass_Kaon_sq mass_Proton_sq GeV/c^2
const double mass_Pion_sq = mass_Pion*mass_Pion, mass_Kaon_sq = mass_Kaon*mass_Kaon,  mass_Proton_sq = mass_Proton*mass_Proton;


TRandom2* rndInt = new TRandom2(1); 

double calcCkovFromMass(double p, double n, double m);
std::array<double, 3> calcCherenkovHyp(double p, double n);


double randomMass(); double randomEnergy();


double randomMomentum()
{
  return 1+4*rndInt->Gaus(0.5, 0.25);
}


// create a number of random particles 
/*double[] randomParticles(int numParticles)
{
  DataSaver dataSaver(Form("RandomParticles%3d.root",numParticles));

  double[] randomMomentum = 

}*/

TH2F* tHistMass = new TH2F("test", "test; Momentum (GeV/c); Cherenkov Angle, #theta_{ch} (rad)", 5000, 0., 5., 800, 0., 0.8);
TCanvas *tCkov = new TCanvas("ckov","ckov",800,800);  
void testHyp()
{  

//TH2F *hClusterMap = new TH2F("Cluster Map", "Cluster Map; x [cm]; y [cm]",1000,-10.,10.,1000,-10.,10.);
  rndInt->SetSeed(0);
  for(double p = 0.; p < 5; p+= 0.001)
  { 

    
    auto photonEnergy = randomEnergy();
    auto n = GetFreonIndexOfRefraction(photonEnergy); // refractive index
    Printf("P =  %f  || n = %f", p, n);
    auto ckovAngles = calcCherenkovHyp(p, n);
    for(auto& ckovAngle:ckovAngles){
      if(!TMath::IsNaN(ckovAngle)){tHistMass->Fill(p, ckovAngle);}
    }
  }
  tCkov->cd();
  tHistMass->Draw();
  //calcCherenkovHyp(1.5, 1.289);

}

std::array<double, 3> ckovMassHyp;
const double ckovConstraintWhidth = 0.015; // 15 mrad per side for ckovangle constraint 
					   // from mass hypotheses


double ckovActual = 0.;
double massActual = 0.;
double pActual = 0.;

void saveDataInst();
/*std::shared_ptr<TFile>*/void saveParticleInfoToROOT(const std::vector<ParticleInfo>& particleVector);

void testRandomMomentum(int numObjects = 10)
{  


// create random mass, energy (and from this ref-index) and momentum :


  
  std::vector<RandomValues> randomObjects(numObjects);
  rndInt->SetSeed(0);

  std::vector<ParticleInfo> particleVector;

  int i = 0;
  for(const auto& randomValue : randomObjects){



     // get cherenkov angle from mass momentum and refindex
     const auto& ckov = calcCkovFromMass(randomValue.momentum, randomValue.refractiveIndex, randomValue.mass); //  calcCkovFromMass(momentum, n, mass)

     // get the map with a given occupancy and ckov angle calculated 
     const auto& map = backgroundStudy(ckov, 0.01);
     Printf("CkovAngle %f Mass %f RefIndex %f Momentum %f", ckov, randomValue.mass, randomValue.refractiveIndex, randomValue.momentum);   

     // make sure the momentum is valid for the given particle (e.g., no - in the square-root in calcCkovFromMass and acos [-1..1])
    if (ckov == 0) {
      continue;
    }
     i++;
     ParticleInfo particle;
     particle.momentum = randomValue.momentum;
     particle.mass = randomValue.mass;
     particle.energy = randomValue.energy;
     particle.refractiveIndex = randomValue.refractiveIndex;
     particle.ckov = ckov;
     particle.map = map;  
     particleVector.emplace_back(particle);
     //map->SaveAs(Form("map%d.root", i));
  }

  // save object
  
  saveParticleInfoToROOT(particleVector);
}




TH2F* backgroundStudy(double ckovActual = 0.5, double occupancy = 0.03)   
{

  auto ckovAngle = ckovActual;

  Int_t NumberOfEvents = 1; Int_t NumberOfClusters = 13; double Hwidth = 15.;
  //testRandomMomentum();
  gStyle->SetOptStat("ei");

  TRandom2* rndP = new TRandom2(1); 
  rndP->SetSeed(0);


  // number of cherenkov photons in the cherenkov ring:
  const auto numberOfCkovPhotons = rndP->Poisson(13);

  photonCandidates.clear();
  double ThetaP=0,PhiP=0,PhiF=0,DegThetaP=0,DegPhiP=0;
  //float /*RadiatorWidth,*/ QuartzWindowWidth,CH4GapWidth,EmissionLenght;
  float FreonIndexOfRefraction,QuartzIndexOfRefraction,CH4IndexOfRefraction;
  double ThetaF1,ThetaF2,ThetaF=0,ThetaLimite;
  float Xpi=0,Ypi=0,Xf=0,Yf=0,Xf1=0,Yf1=0,Xf2=0,Yf2=0,Xp=0,Yp=0; 

  float PhotonEnergy = 6.75; 
  
  FreonIndexOfRefraction = GetFreonIndexOfRefraction(PhotonEnergy);
  QuartzIndexOfRefraction = GetQuartzIndexOfRefraction(PhotonEnergy);
  CH4IndexOfRefraction = 1.00;

  setStyle();

  TH2F *hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",1000,-25.,25.,1000,-25.,25.);
  double mapArray[40][40]{};
   
  float Deltax = (RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght)*TMath::Tan(ThetaP)*TMath::Cos(PhiP);
  float Deltay = (RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght)*TMath::Tan(ThetaP)*TMath::Sin(PhiP);
	  
  Xpi = Xp - Deltax;
  Ypi = Yp - Deltay;
	  
  float ThetaCherenkov[100000] = {0x0}, PhiCherenkov[100000] = {0x0}, DegPhiCherenkov[100000] = {0x0};
    

  
  for(Int_t iEvt = 0; iEvt<NumberOfEvents; iEvt++){
    
    //Printf("event number = %i",iEvt);
    
    float Xcen[100000],Ycen[100000];
     
    DegThetaP = 4.;//0.*(1 - 2*gRandom->Rndm(iEvt));
    gRandom->SetSeed(0);
    DegPhiP   = 360*gRandom->Rndm(iEvt);
        
    ThetaP = TMath::Pi()*DegThetaP/180;
    PhiP = TMath::Pi()*DegPhiP/180;  
     

    const int numBackGround = occupancy*120*120; //occupancy = 0.03
    NumberOfClusters = numBackGround;

    for(Int_t n1=0;n1<NumberOfClusters; n1++) {// clusters loop
      
    //  Printf("cluster = %i",n1);
      
      Xcen[n1] = 60*(1 - 2*gRandom->Rndm(n1));
      Ycen[n1] = 60*(1 - 2*gRandom->Rndm(n1));

      //noiseMap->Fill(Xcen[n1], Ycen[n1]);
      hSignalAndNoiseMap->Fill(Xcen[n1], Ycen[n1]);
      //mapArray[Xcen[n1]+20][Ycen[n1]+20] = 1;
      //hSignalAndNoiseMap2->Fill(Xcen[n1], Ycen[n1]);


     
      
   }
	  
  Int_t NphotTot = 0;
  float MeanTot = 0;
     
 }


 // rndm value in ranfe 0.4, 0.7?
 TRandom2* rnd = new TRandom2(1);
 rnd->SetSeed(0);



 for(Int_t i=0; i < numberOfCkovPhotons; i++) {
   

   double ckovAnglePhot = rnd->Gaus(ckovAngle, 0.012);		    // random CkovAngle, with 0.012 std-dev
   


   // hThetaCh->Fill(ckovAngle);



   double ringRadius = getRadiusFromCkov(ckovAnglePhot); // R in photon-map
   Printf("Cherenkov Photon : Angle = %f Radius = %f", ckovAnglePhot, ringRadius);	



   // populate the photon maps randomly radially
   double alpha = static_cast<double>((3.14159)*(1-2*gRandom->Rndm(1)));    // angle in photon-map (-pi to pi)

   // get x and y values of Photon-candidate:
   double x = static_cast<double>(TMath::Cos(alpha)*ringRadius);
   double y = static_cast<double>(TMath::Sin(alpha)*ringRadius);  



   // populating the pad
   hSignalAndNoiseMap->Fill(x,y);
   // 
   //mapArray[Xcen[n1]+20][Ycen[n1]+20] = 1;

   // add Actual Cherenkov Photon to Candidates
   photonCandidates.emplace_back(ckovAnglePhot);

  } 


	
 //Printf("Hough Window size = %f", Hwidth);
 
 /*auto ckovAnglePredicted = houghResponse(photonCandidates,  Hwidth); */

 return hSignalAndNoiseMap;
 
}
//**********************************************************************************************************************************************************************************************************
float GetFreonIndexOfRefraction(float photonEnergy)
{
  auto x = photonEnergy;
  float k = 1.177 + (0.0172)*x;
  return k;
}
//**********************************************************************************************************************************************************************************************************
float GetQuartzIndexOfRefraction(float x)
  
{
  float k = TMath::Sqrt(1 + 46.411/(113.763556 - x) + 228.71/(328.51563 - x));
  return k;
}

//*********************************************************************************************************************************************************************************************************
double BackgroundFunc(double *x, double *par)
{
 double xx = x[0];
  
 double f = par[0]*TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*(1+TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*TMath::Tan(TMath::ASin(1.2903*TMath::Sin(xx)))*1.2903*TMath::Cos(xx)/cos(asin(1.2903*TMath::Sin(xx))));
  
 return f;
}       
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// get Ring-radius from Cherenkov-angle
double getRadiusFromCkov(double ckovAngle)
{

  //// refIndexFreon refIndexQuartz refIndexCH4
  double sin_ckov = static_cast<double>(TMath::Sin(ckovAngle));
  double sin_qz = static_cast<double>(sin_ckov*(refIndexFreon/refIndexQuartz));
  double sin_theta0 = static_cast<double>(sin_qz*(refIndexQuartz/refIndexCH4));

  double R_ckov = sin_ckov*(RadiatorWidth - EmissionLenght);
  double R_qz = sin_qz * QuartzWindowWidth;
  double R_0 = sin_theta0*CH4GapWidth;
  //Printf("Radiuses  , R_ckov  % f + R_qz %f + R_0 %f", R_ckov, R_qz, R_0);
  double R = static_cast<double>(R_ckov + R_qz + R_0);
  return R;
} 


/*
*/

double /*std::array<TH1D*, 3>*/ houghResponse(std::vector<double>& photonCandidates, double fWindowWidth)
{

  // ckovMassHyp applies constraints!!
  
  auto l0 = ckovMassHyp[0] - ckovConstraintWhidth;
  auto u0 = ckovMassHyp[0] + ckovConstraintWhidth;
  auto l1 = ckovMassHyp[1] - ckovConstraintWhidth;
  auto u1 = ckovMassHyp[1] + ckovConstraintWhidth;
  auto l2 = ckovMassHyp[2] - ckovConstraintWhidth;
  auto u2 = ckovMassHyp[2] + ckovConstraintWhidth;
  Printf("Regions %f %f || %f %f || %f %f  ", l0, u0, l1, u1, l2, u2);

  
  
  /* ef : changed from this:
  // TH1D *resultw = new TH1D("resultw","resultw"       ,nChannels,0,kThetaMax);
  // TH1D *phots   = new TH1D("Rphot"  ,"phots"         ,nChannels,0,kThetaMax);
  // TH1D *photsw  = new TH1D("RphotWeighted" ,"photsw" ,nChannels,0,kThetaMax); */

  int nBin = (int)(kThetaMax / fDTheta);


  // ef : nCorrBand w the previous setting led to the bin1 = 1 bin2 = 750
  int nCorrBand = (int)(fWindowWidth / (2 /** fDTheta*/));

  Printf("nBin %d nCorrBand %d", nBin, nCorrBand);
  int binMax, sumMax = 0;
  std::vector<double> okAngles;
  okAngles.clear();
  for (const auto& angle : photonCandidates) { // photon cadidates loop

    if (angle < 0 || angle > kThetaMax)
      continue;

    // ef : check if angle in ckovMassHyp:
    if (angle < l2 || angle > u0)
      continue;


    phots->Fill(angle);

    int bin = (int)(0.5 + angle / (fDTheta));
    double weight = 1.;
    if (true) {
      double lowerlimit = ((double)bin) * fDTheta - 0.5 * fDTheta;
      double upperlimit = ((double)bin) * fDTheta + 0.5 * fDTheta;


      double rLow = getRadiusFromCkov(lowerlimit);
      double areaLow =  0.5*TMath::Pi()*TMath::Sq(rLow);// calcRingGeom(lowerlimit, 2);
 
      double rHigh = getRadiusFromCkov(upperlimit);
      double areaHigh =  0.5*TMath::Pi()*TMath::Sq(rHigh);// calcRingGeom(lowerlimit, 2);
      //Printf("Areas : areaLow %f, areaHigh %f ", areaLow, areaHigh);

      double diffArea = areaHigh - areaLow;
      
      if (diffArea > 0)
        weight = 1. / diffArea;
    }
    okAngles.emplace_back(angle);
    photsw->Fill(angle, weight);

    int nnn = static_cast<int>(angle*1000);
    arrW[nnn] += weight;
    //fPhotWei.emplace_back(weight); ef: do i need this?
  } // photon candidates loop

  for (int i = 1; i <= nBin; i++) {
    int bin1 = i - nCorrBand;
    int bin2 = i + nCorrBand;
    if (bin1 < 1)
      bin1 = 1;
    if (bin2 > nBin)
      bin2 = nBin;
    double sumPhots = phots->Integral(bin1, bin2);

    /*Printf("bin1 %d ; bin2 %d; sumPhots %f ", bin1, bin2, sumPhots);
    if (sumPhots < 3)
      continue; // if less then 3 photons don't trust to this ring*/
    double sumPhotsw = photsw->Integral(bin1, bin2);
    if ((double)((i /*+ 0.5*/) * fDTheta) > 0.7)
      continue;

    if (sumPhotsw > sumMax){
      binMax = i;
      sumMax = sumPhotsw;
      
    }
    resultw->Fill((double)((i /*+ 0.5*/) * fDTheta), sumPhotsw);
  }
  // evaluate the "BEST" theta ckov as the maximum value of histogramm

  double* pVec = resultw->GetArray();
  int locMax = TMath::LocMax(nBin, pVec);


  double smtest = 0; int ent=0;

  for(const auto& ok:okAngles){
    if (TMath::Abs(ok*1000-locMax) > nCorrBand)
      continue;     
    smtest+=ok; ent++;	
  }
  auto avgTest = smtest/ent;
  Printf("avgTest %f ent %d", avgTest, ent);
  



  Printf("pVec %f locMax %d", *pVec, locMax);
  Printf("sumMax %d, binMax %d", sumMax, binMax);
//photsw
  Printf("Entries : resultw %f photsw %f phots %f", resultw->GetEntries(), photsw->GetEntries(), phots->GetEntries());
  Printf("Max resultw %f photsw %f phots %f",  resultw->GetMaximum(10000.), photsw->GetMaximum(10000.), phots->GetMaximum(10000.));
  Printf("Min resultw %f photsw %f phots %f",  resultw->GetMinimum(-1.), photsw->GetMinimum(-1.), phots->GetMinimum(-1.));

  // ef: not this method, raw-pointers should not be used with new/delete-keywords
  //     smart-pointers are deleted when the fcuntion exits scope :
  // delete phots;delete photsw;delete resultw; // Reset and delete objects

  
  double ckovTrack = static_cast<double>(locMax * fDTheta + 0.5 * fDTheta); // final most probable track theta ckov
  ckovTrackOut = ckovTrack;


  double sumCkov = 0.;
  int entries = 0;
  double ckovSliding = phots->Integral(locMax-nCorrBand, locMax+nCorrBand);

  //binMax
  //locMax
  for(int i = locMax-nCorrBand; i < locMax+nCorrBand; i++)
  {

     auto val = phots->GetBinContent(i);
     auto w = photsw->GetBinContent(i);
     //Printf("Ckov Photons : Bin %d, Value %f; Weight %f", i, val,w);

     if(val > 0.){
       sumCkov += fDTheta*i;//+0.5*fDTheta;
       entries++;
     }
  }

  const double avgCkov = sumCkov/static_cast<double>(entries);
  ckovTrack = avgCkov;
  Printf("Avg Ckov = %f", avgCkov);

  Printf("sumCkov = %f || ckovSliding = %f", sumCkov,ckovSliding);

  TCanvas *houghCanvas = new TCanvas("Hough Canvas","Hough Canvas", 800, 800);
  houghCanvas->Divide(2,2);

  houghCanvas->cd(1);
  auto hThetaCl2 = static_cast<TH1*>(hTheta->Clone());
  setStyleInd(hThetaCl2);
  hThetaCl2->Draw();

  houghCanvas->cd(2);
  setStyleInd(phots);
  phots->Draw();
  TLatex lt2;
  lt2.DrawLatexNDC(.15, .85, Form("Window #eta_{c} :"));
  lt2.DrawLatexNDC(.16, .775, Form("Width = %.0f [mRad]", fWindowWidth));
  lt2.DrawLatexNDC(.16, .7, Form("Entries = %d", entries));

  auto up = getMaxInRange(phots, ckovTrack, fWindowWidth);
  TLine* l3 = new TLine(ckovTrack-(fWindowWidth/2)*0.001, 0, ckovTrack-(fWindowWidth/2)*0.001, up);
  TLine* l4 = new TLine(ckovTrack+(fWindowWidth/2)*0.001, 0, ckovTrack+(fWindowWidth/2)*0.001, up);
  TLine* l5 = new TLine(ckovTrack-(fWindowWidth/2)*0.001,up, ckovTrack+(fWindowWidth/2)*0.001, up);
  l3->SetLineColor(kGreen);  l4->SetLineColor(kGreen);l5->SetLineColor(kGreen);
  l3->SetLineWidth(2);    l4->SetLineWidth(2);l5->SetLineWidth(2);
  l3->Draw();  l4->Draw();l5->Draw();


  Printf("TBox Max %f ", photsw->GetBinContent(phots->GetMaximumBin()));
  Printf("TBox MaxBin %f ", phots->GetMaximumBin());


  int binBox = phots->GetYaxis()->GetLast();
  int yMax = phots->GetYaxis()->GetXmax();
  Printf("TBox binY %d ", binBox);
    Printf("TBox maxY %d ", yMax);

  TBox* box = new TBox(/*ckovTrack*/ckovTrack-nCorrBand*0.001,0.,/*ckovTrack*/ckovTrack+nCorrBand*0.001,phots->GetBinContent(ckovTrack*1000/*photsw->GetMaximumBin()*/));
  box->SetLineColor(kGreen);
  box->SetFillStyle(0);  
  box->SetLineWidth(2);  


  houghCanvas->cd(3);
  setStyleInd(photsw);
  double up_ = 0;
  auto th = getMaxInRange(photsw, up_, ckovTrack, fWindowWidth);
  th->Draw();




  auto l3c = static_cast<TLine*>(l3->Clone());  auto l4c = static_cast<TLine*>(l4->Clone());  auto l5c = static_cast<TLine*>(l5->Clone());
  l5c->SetY1(up_);
  l3c->SetY2(up_); l4c->SetY2(up_);l5c->SetY2(up_);
  l3c->Draw();  l4c->Draw();l5c->Draw();



  auto pad4 = static_cast<TPad*>(houghCanvas->cd(4));
  pad4->PaintText(.2, .9, Form("Track Cherenkov"));
  pad4->PaintText(.2, .8, Form("Actual Value = %.3f", meanCherenkovAngle));
  pad4->PaintText(.2, .7, Form("Predicted Value = %.3f", ckovTrack));
  pad4->PaintTextNDC(.2, .9, Form("Track Cherenkov"));
  pad4->PaintTextNDC(.2, .8, Form("Actual Value = %.3f", meanCherenkovAngle));
  pad4->PaintTextNDC(.2, .7, Form("Predicted Value = %.3f", ckovTrack));
  
  gStyle->SetPaintTextFormat("1.3f");
  setStyleInd(resultw);
  resultw->Draw();


  TLine* l = new TLine(/*ckovTrack*/ckovTrack, 0, /*ckovTrack*/ckovTrack, resultw->GetBinContent(locMax)*2);

  //TLine* l = new TLine(resultw->GetMaximumBin(), 0., resultw->GetMaximumBin(), resultw->GetBinContent(resultw->GetMaximumBin()));
  l->SetLineColor(kGreen);
  l->SetLineWidth(2);  
  l->Draw();

   
  Printf("ckovTrack = %f", ckovTrack);
  return ckovTrack;

}




void setStyleInd(TH1* th1f, double ratio = 1.2)
{
  th1f->SetTitleSize((th1f->GetTitleSize("x")*ratio), "xy");
  th1f->SetLabelSize((th1f->GetLabelSize("x")*ratio), "xy");
}


/*
void setPad(TPad*, double l, double, r, double t, double b)
{
  th1f->SetPadMaring((th1f->GetTitleSize("x")*ratio), "xy");
}*/


void setStyleInd(TH2* th1f, double ratio)
{
  th1f->SetTitleSize((th1f->GetTitleSize("x")*ratio), "xy");
  th1f->SetLabelSize((th1f->GetLabelSize("x")*ratio), "xy");
}





double getMaxInRange(TH1* th1, double mid, double width)
{


  //Printf("mid %f || width %f ", mid, width);
  double max = 1.0;

  const int startBin = static_cast<int>(mid*1000);
  const int nBin2 = static_cast<int>(width/2);
  //Printf("startBin %d || endBin %d ", startBin-nBin2, startBin+nBin2);


  int start = static_cast<int>(startBin-nBin2);
  int end = static_cast<int>(startBin+nBin2);
  for(int i = start; i < end; i++){
    auto binEnt = th1->GetBinContent(i);
    //if (binEnt > max) max = binEnt;
    //Printf("ent i %d || val %f ", i, binEnt);
  }
  return max;
}

TH1* getMaxInRange(TH1* th1, double& up, double mid, double width)
{
  TH1* thOut = static_cast<TH1*>(th1);

  //Printf("mid %f || width %f ", mid, width);
  double max = 1.0;

  const int startBin = static_cast<int>(mid*1000);
  const int nBin2 = static_cast<int>(width/2);
  //Printf("startBin %d || endBin %d ", startBin-nBin2, startBin+nBin2);


  int start = static_cast<int>(startBin-nBin2);
  int end = static_cast<int>(startBin+nBin2);

  double r = 0;
  for(const auto& i : arrW){
    if (i > r) 
      r = i;
  }
  
  for(int i = start; i < end; i++){
    auto binEnt = th1->GetBinContent(i);
    auto binent = arrW[i-1];
    thOut->SetBinContent(binent, i);
    //if (binEnt > max) max = binent;
    //Printf("ent i %d || val %f || val2 %f", i, binEnt, binent);

  }
  thOut->GetXaxis()->SetRangeUser(0., r);
  up = max;
  return thOut;
}

// mass_Pion_sq mass_Kaon_sq mass_Proton_sq
std::array<double, 3> calcCherenkovHyp(double p, double n)
{
  const double p_sq = p*p;
  const double cos_ckov_denom = p*refIndexFreon;
  const auto cos_ckov_Pion = TMath::Sqrt(p_sq + mass_Pion_sq)/(cos_ckov_denom); // n = refIndexFreon 1.289 later make it random?

  const auto cos_ckov_Kaon = TMath::Sqrt(p_sq + mass_Kaon_sq)/(cos_ckov_denom); 
  const auto cos_ckov_Proton = TMath::Sqrt(p_sq + mass_Proton_sq)/(cos_ckov_denom);
  

  const auto ckovAnglePion = TMath::ACos(cos_ckov_Pion); 
  const auto ckovAngleKaon = TMath::ACos(cos_ckov_Kaon); 
  const auto ckovAngleProton = TMath::ACos(cos_ckov_Proton); 

  Printf("Pion %.3f Kaon %.3f Proton %.3f", ckovAnglePion, ckovAngleKaon, ckovAngleProton);

  return {ckovAnglePion, ckovAngleKaon, ckovAngleProton};
}



double calcCkovFromMass(double p, double n, double m)
{
  const double p_sq = p*p;
  const double cos_ckov_denom = p*refIndexFreon;

  // sanity check ;)
  if(p_sq + m*m < 0){
    return 0;
  }

  const auto cos_ckov = TMath::Sqrt(p_sq + m*m)/(cos_ckov_denom);

  // sanity check ;)
  if(cos_ckov > 1 || cos_ckov < -1)
    return 0;

  const auto ckovAngle = TMath::ACos(cos_ckov);

  return ckovAngle;
}




double randomMass() 
{  
  auto index = static_cast<int>(rndInt->Integer(3));
  Printf("randomMass indes = %d", index);
  return masses[index];
}

double randomEnergy()
{
  // random energy from the arrWaveLenDefault
  auto index = static_cast<int>(rndInt->Integer(30));
  Printf("rEn indes = %d", index);
  double photonEnergy = static_cast<double>(nm2eV/arrWaveLenDefault[index]);
  return photonEnergy;
}

void setStyle()
{
  gStyle->SetOptStat("ei");
  phots->SetTitleSize(phots->GetTitleSize("x")*1.3, "xy");
  phots->SetLabelSize(phots->GetLabelSize("x")*1.3, "xy");

  photsw->SetTitleSize(phots->GetTitleSize("x"), "xy");
  photsw->SetLabelSize(phots->GetLabelSize("x"), "xy");

  resultw->SetTitleSize(phots->GetTitleSize("x"), "xy");
  resultw->SetLabelSize(phots->GetLabelSize("x"), "xy");

  hTheta->SetTitleSize(phots->GetTitleSize("x"), "xy");
  hTheta->SetLabelSize(phots->GetLabelSize("x"), "xy");
}


void saveDataVector()
{
  std::vector<DataInfo> dataVector(100);  // Assuming you have a vector of Data

    // Fill the vector with some data
    // In a real case, you would likely fill this from your actual data source
    for (int i = 0; i < 100; i++) {
        dataVector[i].momentum = i * 0.5;  // some dummy values
        dataVector[i].typeOfParticle = i % 3;
        dataVector[i].refractiveIndex = i * 0.1;
        for(int j=0; j<10; j++){
            for(int k=0; k<10; k++){
                dataVector[i].pads[j][k] = i+j+k; // some dummy values
            }
        }
    }

    DataSaver dataSaver("data.root");
    dataSaver.fillData(dataVector);
    dataSaver.save();
}


void saveDataInst()
{

  //std::vector<int> dataVector(100);  // Assuming you have a vector of Data	
  std::vector<DataInfo> dataVector(100);  // Assuming you have a vector of Data

    // Fill the vector with some data
    // In a real case, you would likely fill this from your actual data source
    for (int i = 0; i < 100; i++) {
        dataVector[i].momentum = i * 0.5;  // some dummy values
        dataVector[i].typeOfParticle = i % 3;
        dataVector[i].refractiveIndex = i * 0.1;
        for(int j=0; j<10; j++){
            for(int k=0; k<10; k++){
                dataVector[i].pads[j][k] = i+j+k; // some dummy values
            }
        }
    }

    DataSaver dataSaver("data.root");
    dataSaver.fillData(dataVector);
    dataSaver.save();
}


/*
// save TH2F* in own TTree since it caused segmentation fault when writing to same TTree as the other elements
std::shared_ptr<TFile> saveParticleInfoToROOT(const std::vector<ParticleInfo>& particleVector) {
    // Create a smart pointer for the TFile
    std::shared_ptr<TFile> outputFile(new TFile("outputFile.root", "RECREATE"));

    // Create a smart pointer for the TTree
    std::shared_ptr<TTree> tree(new TTree("tree", "ParticleInfo Tree"));

    // Create variables to hold the values of ParticleInfo properties
    double momentum;
    double mass;
    double energy;
    float refractiveIndex;
    double ckov;

    // Set the branch addresses for the TTree
    tree->Branch("momentum", &momentum);
    tree->Branch("mass", &mass);
    tree->Branch("energy", &energy);
    tree->Branch("refractiveIndex", &refractiveIndex);
    tree->Branch("ckov", &ckov);

    // Loop over the ParticleInfo objects and fill the TTree
    for (const auto& particle : particleVector) {
        momentum = particle.momentum;
        mass = particle.mass;
        energy = particle.energy;
        refractiveIndex = particle.refractiveIndex;
        ckov = particle.ckov;

        tree->Fill();

        // Ensure the map is not a nullptr before writing it to the file
        if (particle.map) {
            particle.map->Write(); // This will write the TH2F* to the file
        }
    }

    // Write the TTree to the TFile
    outputFile->cd();
    tree->Write();

    return outputFile;
}*/ 

void saveParticleInfoToROOT2(const std::vector<ParticleInfo>& particleVector) {
    // Create a new TFile
    TFile* outputFile = new TFile("outputFile2.root", "RECREATE");

    // Create a TTree
    TTree* tree = new TTree("tree", "ParticleInfo Tree");

    // Create a ParticleInfo pointer to hold the values
    ParticleInfo* particleInfo = nullptr;

    // Create a branch for ParticleInfo objects
    tree->Branch("ParticleInfo", "ParticleInfo", &particleInfo, 32000, 0);

    // Loop over the ParticleInfo objects and fill the TTree
    for (const auto& particle : particleVector) {
        // Delete the old object and create a new one with the current particle's values
        delete particleInfo;
        particleInfo = new ParticleInfo(particle);

        // Fill the tree
        tree->Fill();
    }

    // Write the TTree to the TFile
    tree->Write();

    // Close the TFile
    outputFile->Close();

    // Delete the last object created
    delete particleInfo;
}


void readParticleInfoFromROOT() {
    // Open the file
    TFile* inputFile = new TFile("outputFile.root", "READ");

    // Get the TTree
    TTree* tree = (TTree*)inputFile->Get("tree");

    // Variables to hold the values of ParticleInfo properties
    double momentum;
    double mass;
    double energy;
    float refractiveIndex;
    double ckov;

    // Set the branch addresses for the TTree
    tree->SetBranchAddress("momentum", &momentum);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("refractiveIndex", &refractiveIndex);
    tree->SetBranchAddress("ckov", &ckov);

    // Get number of entries in the TTree
    Long64_t nEntries = tree->GetEntries();

    // Loop over the entries in the TTree
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        // Load the entry
        tree->GetEntry(iEntry);

        // Print the ParticleInfo object's values
        std::cout << "ParticleInfo object " << iEntry << ":\n";
        std::cout << "  momentum: " << momentum << "\n";
        std::cout << "  mass: " << mass << "\n";
        std::cout << "  energy: " << energy << "\n";
        std::cout << "  refractiveIndex: " << refractiveIndex << "\n";
        std::cout << "  ckov: " << ckov << "\n";

        // Go to the "maps" directory
        TDirectory* mapsDir = (TDirectory*)inputFile->Get("maps");
        mapsDir->cd();

        // Get the associated map
        TString mapName = TString::Format("hist_%lld", iEntry);
        TH2F* map = (TH2F*)gDirectory->Get(mapName);

        if (map) {
            // Print the number of entries in the map
            std::cout << "  map entries: " << map->GetEntries() << "\n";
        } else {
            std::cout << "  map not found\n";
        }

        // Go back to the top directory
        inputFile->cd();
    }

    // Close the file
    inputFile->Close();
}

/*
void readParticleInfoFromROOT2() {
    // Open the file
    TFile* inputFile = new TFile("outputFile2.root", "READ");

    // Get the TTree
    TTree* tree = (TTree*)inputFile->Get("tree");

    // Create a ParticleInfo pointer
    ParticleInfo* particleInfo = nullptr;

    // Set the branch address
    Printf("Setting ParticleInfo adress");
    tree->SetBranchAddress("ParticleInfo", &particleInfo);

    // Get number of entries in the TTree
    Long64_t nEntries = tree->GetEntries();

    Printf("Entering Loop");

    // Loop over the entries in the TTree
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        // Load the entry


        if(particleInfo) {
            delete particleInfo;
            particleInfo = nullptr;
        }
        std::cout << "tree->GetEntry(iEntry); " << iEntry << ":\n";
        tree->GetEntry(iEntry);

        if(particleInfo) {
            std::cout << "ParticleInfo object " << iEntry << ":\n";
            std::cout << "  momentum: " << particleInfo->momentum << "\n";
            std::cout << "  mass: " << particleInfo->mass << "\n";
            std::cout << "  energy: " << particleInfo->energy << "\n";
            std::cout << "  refractiveIndex: " << particleInfo->refractiveIndex << "\n";
            std::cout << "  ckov: " << particleInfo->ckov << "\n";
        } else {
            std::cout << " ParticleInfo object not found\n";
        }

        // Get the associated histogram
        TString histName = TString::Format("hist_%lld", iEntry);
        TH2F* hist = (TH2F*)inputFile->Get(histName);
        if (hist) {
            // Print the number of points in the histogram
            std::cout << "  map points: " << hist->GetEntries() << "\n";
        } else {
            std::cout << "  map not found\n";
        }
    }

    // Close the file
    inputFile->Close();
}*/ 



void saveParticleInfoToROOT(const std::vector<ParticleInfo>& particleVector) {
    // Create a new TFile
    TFile* outputFile = new TFile("outputFile.root", "RECREATE");

    // Create a TTree
    TTree* tree = new TTree("tree", "ParticleInfo Tree");

    // Create variables to hold the values of ParticleInfo properties
    double momentum;
    double mass;
    double energy;
    float refractiveIndex;
    double ckov;

    // Set the branch addresses for the TTree
    tree->Branch("momentum", &momentum, "momentum/D");
    tree->Branch("mass", &mass, "mass/D");
    tree->Branch("energy", &energy, "energy/D");
    tree->Branch("refractiveIndex", &refractiveIndex, "refractiveIndex/F");
    tree->Branch("ckov", &ckov, "ckov/D");

    // Create a directory to store the maps
    TDirectory* mapsDir = outputFile->mkdir("maps");
    mapsDir->cd();

    // Counter for unique histogram names
    int histCounter = 0;

    // Loop over the ParticleInfo objects and fill the TTree
    for (const auto& particle : particleVector) {
        momentum = particle.momentum;
        mass = particle.mass;
        energy = particle.energy;
        refractiveIndex = particle.refractiveIndex;
        ckov = particle.ckov;
	//TCanvas* canvas = new TCanvas("canvas", "Map Canvas", 800, 800);
        // Write each histogram to the maps directory with a unique name
        TString histName = TString::Format("hist_%d", histCounter++);
        TH2F* histCopy = new TH2F(*particle.map);
        histCopy->SetName(histName);

	// Draw the map on the canvas with a 1:1 aspect ratio
	//canvas->SetCanvasSize(800, 800);
	//canvas->SetFixedAspectRatio();
        histCopy->Write();
        tree->Fill();
    }

    // Go back to the top directory
    outputFile->cd();

    // Write the TTree to the TFile
    tree->Write();

    // Close the TFile
    outputFile->Close();


    // save by class: some problems
    //saveParticleInfoToROOT2(particleVector);


    // works good
    Printf("\n\n Reading from file now...");
    readParticleInfoFromROOT();
}




