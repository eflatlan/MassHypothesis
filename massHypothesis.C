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


double arrW[750]= {0.};


using namespace o2;
using namespace o2::hmpid;

#include <HMPIDBase/Param.h>
void setStyleInd(TH2* th1f, double ratio = 1.2);
void setStyleInd(TH1* th1f, double ratio = 1.2);

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




const double CH4GapWidth = 8;
const double  RadiatorWidth = 1.;
const double  QuartzWindowWidth = 0.5;
const double  EmissionLenght = RadiatorWidth/2;


TH1* getMaxInRange(TH1* th1, double& up, double mid, double width);
double getMaxInRange(TH1* th1, int start, int width);
double getMaxInRange(TH1* th1, double mid, double width);


// mass_Pion_sq mass_Kaon_sq mass_Proton_sq
const double mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in GeV/c^2

const double mass_Pion_sq = mass_Pion*mass_Pion, mass_Kaon_sq = mass_Kaon*mass_Kaon,  mass_Proton_sq = mass_Proton*mass_Proton;


TRandom2* rndInt = new TRandom2(1); 

double calcCkovFromMass(double p, double n, double m);
std::array<double, 3> calcCherenkovHyp(double p, double n);

std::array<double, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
double randomMass(); double randomEnergy();


double randomMomentum()
{
  return 1+4*rndInt->Gaus(0.5, 0.25);
}

TH2F* tHistMass = new TH2F("test", "test; Momentum (GeV/c); Cherenkov Angle, #theta_{ch} (rad)", 5000, 0., 5., 800, 0., 0.8);
TCanvas *tCkov = new TCanvas("ckov","ckov",800,800);  
void testHyp()
{  

//TH2F *hClusterMap = new TH2F("Cluster Map", "Cluster Map; x [cm]; y [cm]",1000,-10.,10.,1000,-10.,10.);
  rndInt->SetSeed(0);
  for(double p = 0.; p < 5; p+= 0.001)
  { 

    
    auto photonEnergy = randomEnergy();
    auto n = GetFreonIndexOfRefraction(photonEnergy);
    Printf("P =  %f  || n = %f", p, n);
    auto ckovAngles = calcCherenkovHyp(p, n);
    for(auto& ckovAngle:ckovAngles){
      if(!TMath::IsNaN(tt)){tHistMass->Fill(p, ckovAngle);}
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

void testRandomMomentum()
{  

//TH2F *hClusterMap = new TH2F("Cluster Map", "Cluster Map; x [cm]; y [cm]",1000,-10.,10.,1000,-10.,10.);
  rndInt->SetSeed(0);


  // create random momentum (model has this information):
  momentum = randomMomentum();


  // the actual mass (model is blind to this):
  mass = randomMass();

  // "random" refractive index (known to model)
  auto photonEnergy = randomEnergy();
  auto n = GetFreonIndexOfRefraction(photonEnergy);
  
  // create the ckovAngles from the mass hypotheses:
  auto ckovAngles = calcCherenkovHyp(p, n);
  ckovMassHyp = ckovAngles;
  

  ckovActual = calcCkovFromMass(mass);
  //calcCherenkovHyp(1.5, 1.289);

}

void backgroundStudy(Int_t NumberOfEvents, Int_t NumberOfClusters, double Hwidth, double occupancy = 0.03)   
{
    gStyle->SetOptStat("ei");

  TRandom2* rndP = new TRandom2(1); 
  rndP->SetSeed(0);
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
  

  TH1F *ThetaHough = new TH1F("ThetaHough","ThetaHough",1000,0,1);
  TH1F *NPhoton    = new TH1F("NPhoton","NPhoton",20,0.,20.); 
  TH1F *NHough     = new TH1F("NHuogh","NHough",1000,0.,1.); 

  setStyle();


  TH1F *hTheta2    = new TH1F("Photon Cherenkov Angle Distribution","Photon Cherenkov Angle Distribution;angle [rad]; counts/1 mrad",750,0.,0.75);
  TH1F *hThetawg   = new TH1F("Track Cherenkov Angle Distribution","Track Cherenkov Angle Distribution;angle [rad]; counts/1 mrad",750,0.,0.75);
  TH1F *hThetaRing = new TH1F("hThetaRing","Ring Cherenkov; Cherenkov Angle [rad];",1000,0.,1.);

  TH2F *hClusterMap = new TH2F("Cluster Map", "Cluster Map; x [cm]; y [cm]",1000,-10.,10.,1000,-10.,10.);


  TH2F *hSignalAndNoiseMap = new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]",1000,-25.,25.,1000,-25.,25.);

TH2F *hSignalAndNoiseMap2 = new TH2F("Signal and Noise 2", "Signal and Noise 2; x [cm]; y [cm]",1000,-10.,10.,1000,-10.,10.);


  TH2F *hphotonMap = new TH2F("Photon Map","Photon Map; x [cm]; y [cm]",1000,-10.,10.,1000,-10.,10.);

  TH2F *hphotonMap2 = new TH2F("Photon Map2","Photon Map2; x [cm]; y [cm]",1000,-25.,25.,1000,-25.,25.);


  TH2F *signalMap = new TH2F("Hit Map","Hit Map; x [cm]; y [cm]",1000,-7.,7.,1000,-7.,8.);
  TH2F *noiseMap = new TH2F("noiseHitMap","noiseHitMap; x [cm]; y [cm]",1000,-7.,7.,1000,-7.,8.);
 
   
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
      hphotonMap->Fill(Xcen[n1], Ycen[n1]);
      hphotonMap2->Fill(Xcen[n1], Ycen[n1]);



      noiseMap->Fill(Xcen[n1], Ycen[n1]);
      hSignalAndNoiseMap->Fill(Xcen[n1], Ycen[n1]);
      hSignalAndNoiseMap2->Fill(Xcen[n1], Ycen[n1]);

//      Xcen[n1] = 130*gRandom->Rndm(n1); Ycen[n1] = 130*gRandom->Rndm(n1); 
      
      TVector3 v2(Xcen[n1]-Xpi-EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP),Ycen[n1]-Ypi-EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP),RadiatorWidth+QuartzWindowWidth+CH4GapWidth-EmissionLenght); 
      
      PhiF = v2.Phi();     



     // R_ckov = sin (theta_ckov)*(RadiatorWidth - EmissionLenght)
     // R_qz = sin(theta_qz) * QuartzWindowWidth
     // R_0 = sin(theta_0)*CH4GapWidth
     
	




     //theta_qz = TMath::Asin(sin_ckov*n_GAP/n_quartz);

     // sin_ckov * n_Radiator = sin_qz * n_quartz
     // sin_qz * n_quartz = sin_theta0 * n_GAP
       
      
      ThetaLimite = TMath::ASin(CH4IndexOfRefraction/QuartzIndexOfRefraction);
      
      double ThetaF0 = TMath::ASin(QuartzIndexOfRefraction/FreonIndexOfRefraction*TMath::Sin(ThetaLimite))-0.00001;


      
    //  Printf("ThetaF0 = %f",ThetaF0*TMath::RadToDeg());
	      
      double ThetaF01 = TMath::ASin((FreonIndexOfRefraction/QuartzIndexOfRefraction)*(TMath::Sin(ThetaF0)));      
      
      double ThetaF02 = TMath::ASin((QuartzIndexOfRefraction/CH4IndexOfRefraction)*(TMath::Sin(ThetaF01)));
	      
      float X01 = EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP);
	      
      float Y01 =  EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP);
	      
      float X02 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF0)*TMath::Cos(PhiF)+QuartzWindowWidth*TMath::Tan(ThetaF01)*TMath::Cos(PhiF)+CH4GapWidth*TMath::Tan(ThetaF02)*TMath::Cos(PhiF);
	      
      float Y02 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF0)*TMath::Sin(PhiF) + QuartzWindowWidth*TMath::Tan(ThetaF01)*TMath::Sin(PhiF) + CH4GapWidth*TMath::Tan(ThetaF02)*TMath::Sin(PhiF);  
	  
      float X0 = X01 + X02;
      float Y0 = Y01 + Y02;
	      
      double ThetaMin = 0;
//      double ThetaMax = 0.75+ThetaP;
      double ThetaMax = ThetaF0;
	      
      Xf = 999;
      Yf = 999;
      
      Int_t nWhile = 0;
      
      while(TMath::Sqrt((Xf-Xcen[n1]+Xpi)*(Xf-Xcen[n1]+Xpi)+(Yf-Ycen[n1]+Ypi)*(Yf-Ycen[n1]+Ypi))>0.0001)
		
	{ 
          nWhile++;
          
	  ThetaF = (double) (0.5*(ThetaMax - ThetaMin) + ThetaMin);
	  
	  ThetaF1 = TMath::ASin((FreonIndexOfRefraction/QuartzIndexOfRefraction)*(TMath::Sin(ThetaF)));     
	  
	  ThetaF2 = TMath::ASin((QuartzIndexOfRefraction/CH4IndexOfRefraction)*(TMath::Sin(ThetaF1)));
	  
	  Xf1 = EmissionLenght*TMath::Tan(ThetaP)*TMath::Cos(PhiP);
		  
	  Yf1 =  EmissionLenght*TMath::Tan(ThetaP)*TMath::Sin(PhiP);
	  
	  Xf2 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF)*TMath::Cos(PhiF)+QuartzWindowWidth*TMath::Tan(ThetaF1)*TMath::Cos(PhiF)+CH4GapWidth*TMath::Tan(ThetaF2)*TMath::Cos(PhiF);
	  
		  
	  Yf2 = (RadiatorWidth - EmissionLenght)*TMath::Tan(ThetaF)*TMath::Sin(PhiF) + QuartzWindowWidth*TMath::Tan(ThetaF1)*TMath::Sin(PhiF) + CH4GapWidth*TMath::Tan(ThetaF2)*TMath::Sin(PhiF);  
		  

	 

	  Xf = Xf1 + Xf2;
	  Yf = Yf1 + Yf2;
	  //Printf("Positions Xf1 %f Yf1 %f Xf2 %f Yf2 %f" , Xf1, Yf1, Xf2, Yf2);
		  


	  if(TMath::Sqrt((Xf-X0)*(Xf-X0)+(Yf-Y0)*(Yf-Y0))>TMath::Sqrt((Xcen[n1]-Xpi-X0)*(Xcen[n1]-Xpi-X0)+(Ycen[n1]-Ypi-Y0)*(Ycen[n1]-Ypi-Y0)))
		    
	    {
	      ThetaMin = ThetaF;
	    }
		  
	  else 
		    
	    {
	      ThetaMax = ThetaF;   
	    }  
            
            if(nWhile>30) break;
            
	} // while 
        
	      
      TVector3 vP((TMath::Sin(ThetaP))*(TMath::Cos(PhiP)),(TMath::Sin(ThetaP))*(TMath::Sin(PhiP)),(TMath::Cos(ThetaP)));
      TVector3 vz(0.,0.,1.);
	      
      TVector3 v1 = vP.Cross(vz);
	      
      TVector3 vF((TMath::Sin(ThetaF))*(TMath::Cos(PhiF)),(TMath::Sin(ThetaF))*(TMath::Sin(PhiF)),(TMath::Cos(ThetaF)));
      
      if(ThetaP==0)
	{	      
	  ThetaCherenkov[n1] = ThetaF;		  
	  PhiCherenkov[n1] = PhiF;	      
	}	  
      else		
	{
	  vF.Rotate(ThetaP,v1);
	  
	  ThetaCherenkov[n1] = vF.Theta();
	  PhiCherenkov[n1] = vF.Phi();
	}
      
      DegPhiCherenkov[n1] = 180*PhiCherenkov[n1]/(TMath::Pi());

      if(DegPhiCherenkov[n1]<0) DegPhiCherenkov[n1]+=360;	      

      TVector2 v5(Xcen[n1]-Xp,Ycen[n1]-Yp);
      
      //double Phi = v5.Phi();
      
      //float DegPhi = 180*Phi/TMath::Pi(); 
      
      hTheta->Fill(ThetaCherenkov[n1]);
      // add background Photon to Candidates
      photonCandidates.emplace_back(ThetaCherenkov[n1]);
	      	      
      //if(k1==2233) cout << " ThetaCherenkov " <<ThetaCherenkov[n1] <<"   PhiCherenkov = "<<DegPhiCherenkov[n1]<<"  event number = "<<k1<< endl;

      Int_t dstep = (Int_t)DegPhiCherenkov[n1]/20;
      if(dstep==18) dstep = 0;
      if(dstep<0 || dstep>18) continue;	      
      
   }//clusters loop 	  
  
  Int_t HCS[700] = {0x0};
	  
  for(Int_t k2=0;k2<700;k2++)
    {
      
      HCS[k2] = 0;             
      Int_t NphotHough = 0;
      
      for(Int_t p=0;p<NumberOfClusters;p++)
	{
		  
	  if(ThetaCherenkov[p]>(0.001*k2) && ThetaCherenkov[p]<(0.001*Hwidth+0.001*k2)) NphotHough++;
	}
	      
      HCS[k2] = NphotHough;
	      
      NHough->Fill(0.001*k2+0.001/2.,(float)NphotHough);
	      
   }
	  
  Int_t LocPos = TMath::LocMax(700,HCS);
	  
  Int_t NphotTot = 0;
  float MeanTot = 0;
  
  for(Int_t p=0;p<NumberOfClusters;p++)
    {
      if(ThetaCherenkov[p]>(0.001*LocPos) && ThetaCherenkov[p]<(Hwidth+0.001*LocPos)) 
	{
	  ThetaHough->Fill(ThetaCherenkov[p]);
		  
	  NphotTot++;
	  MeanTot+=ThetaCherenkov[p];		 
	  // 	  ThetaPhoton[NphotTot-1]=ThetaCherenkov[p];
	}
    }
	    
   float RingThetaCherenkov = MeanTot/(float)NphotTot;

   NPhoton->Fill(NphotTot);
   
   if(NphotTot<3) continue;
   
   hThetaRing->Fill(RingThetaCherenkov); // gjor denne på weighted m signal
        
 }// events loop

 // endre weight til fra Recon fra Giacomo sin branch

 TF1 *fBackGround = new TF1("fBackGround",BackgroundFunc,0,0.75,1);


 auto hThetaClone = static_cast<TH1*>(hTheta->Clone());
 hTheta->Fit(fBackGround,"RQ");
 //hTheta->Draw();



 

 
 Printf("par0 = %f",fBackGround->GetParameter(0));
 


 // rndm value in ranfe 0.4, 0.7?
 TRandom2* rnd = new TRandom2(1);
 rnd->SetSeed(0);

 const double ckovAngleMean = 0.4+0.25*rnd->Gaus(0.5, 0.25);
 meanCherenkovAngle = ckovAngleMean;
 Printf("ckovAngleMean Generated %f", ckovAngleMean);


 for(Int_t i=0; i < numberOfCkovPhotons; i++) {
   

   
   double ckovAngle = rnd->Gaus(ckovAngleMean, 0.012);		    // random CkovAngle
   


  hThetaCh->Fill(ckovAngle);



   double ringRadius = getRadiusFromCkov(ckovAngle); // R in photon-map
   Printf("Cherenkov Photon : Angle = %f Radius = %f", ckovAngle, ringRadius);	

   double alpha = static_cast<double>((3.14159)*(1-2*gRandom->Rndm(1)));    // angle in photon-map (-pi to pi)

   // get x and y values of Photon-candidate:
   double x = static_cast<double>(TMath::Cos(alpha)*ringRadius);
   double y = static_cast<double>(TMath::Sin(alpha)*ringRadius);  


   signalMap->Fill(x,y);//signalMap noiseMap hSignalAndNoiseMap
   hSignalAndNoiseMap->Fill(x,y);
   hSignalAndNoiseMap2->Fill(x,y);


 // hSignalAndNoiseMap, hSignalAndNoiseMap2 i s

   // add Actual Cherenkov Photon to Candidates
   photonCandidates.emplace_back(ckovAngle);

   hTheta2->Fill(ckovAngle);
   hThetawg->Fill(ckovAngle);
  } 
 
 for(Int_t n2=0;n2<NumberOfClusters; n2++) {
  
   if(ThetaCherenkov[n2]>0.75) continue;
   
   Int_t bin = (Int_t)(ThetaCherenkov[n2]/(0.001))+1;
   
  // Printf("bin = %i, hTheta->GetBinContent(bin)=%f",bin, hTheta->GetBinContent(bin));
   
  // if(hTheta->GetBinContent(bin)==0) bin=bin+1;
   
   double weight = 1 - fBackGround->Eval(ThetaCherenkov[n2])/hTheta->GetBinContent(bin);

   hTheta2->Fill(ThetaCherenkov[n2]);
    
   hThetawg->Fill(ThetaCherenkov[n2],weight);


 }

	
 Printf("Hough Window size = %f", Hwidth);
 
 auto ckovAnglePredicted = houghResponse(photonCandidates,  Hwidth);
 
 Printf("Total Number of photonCandidates = %zu", photonCandidates.size()); 
 Printf("Total Number of Background Clusters = %i", NumberOfClusters); 
 Printf("Number of  numberOfCkovPhotons = %i", numberOfCkovPhotons); 
 Printf("Cherenkov Photon Mean Angle = %f " , ckovAngleMean);
 Printf("Predicted Cherenkov Photon Angle = %f " , ckovAnglePredicted);
 


 TCanvas *signalNoiseCanvas = new TCanvas("Signal and Noise Hit Map","Signal and Noise Hit Map",800,800);  
 signalNoiseCanvas->Divide(2,1);
 signalNoiseCanvas->cd(1);
 signalMap->Draw();
 signalNoiseCanvas->cd(2);
 noiseMap->Draw();



 TCanvas *signalPNoiseCanvas = new TCanvas("Signal+Noise Hit Map","Signal+Noise Hit Map",800,800);
 signalPNoiseCanvas->Divide(2,1);
 signalPNoiseCanvas->cd(1);
 setStyleInd(noiseMap);
 noiseMap->SetMarkerStyle(kBlack);
 noiseMap->SetMarkerStyle(2);
 noiseMap->SetMarkerColor(kRed);
 noiseMap->SetTitle(Form("Photon  Hit Map"));


 signalMap->SetMarkerStyle(3);
 signalMap->SetMarkerColor(kBlue);



 for(double d = 0.1; d < 0.4; d+= 0.05){
   double rmin = getRadiusFromCkov(d+0.001*Hwidth/2); // R in photon-map
   double rmax = getRadiusFromCkov(d-0.001*Hwidth/2); // R in photon-map
   Printf("d%f Rmin%f, Rmax%f", d, rmin, rmax);
 }



 double rMax = getRadiusFromCkov(ckovTrackOut+0.001*Hwidth/2); // R in photon-map
 double rMin = getRadiusFromCkov(ckovTrackOut-0.001*Hwidth/2); // R in photon-map
 TEllipse* telMin = new TEllipse(0., 0., rMin);
 TEllipse* telMax = new TEllipse(0., 0., rMax);

 TArc* arcMin = new TArc(0., 0., rMin);
 TArc* arcMax = new TArc(0., 0., rMax);

 telMin->SetLineColor(kGreen); telMax->SetLineColor(kGreen);
 

 telMin->SetFillStyle(0);telMax->SetFillStyle(0);
 telMin->SetLineWidth(2);telMax->SetLineWidth(2);


  auto up1 = getMaxInRange(hThetaClone, ckovTrackOut, Hwidth);
  auto up2 = getMaxInRange(hThetaCh, ckovTrackOut, Hwidth);
  auto up = std::max(up1, up2);

 TBox* box2 = new TBox(ckovTrackOut-(Hwidth/2)*0.001,0.,ckovTrackOut+(Hwidth/2)*0.001,up);
 box2->SetLineColor(kGreen);
 box2->SetFillStyle(0);  
 box2->SetLineWidth(2);  
  



 noiseMap->GetXaxis()->SetRangeUser(-7, 7);
 noiseMap->GetYaxis()->SetRangeUser(-7, 8); 

 signalMap->GetXaxis()->SetRangeUser(-7, 7);
 signalMap->GetYaxis()->SetRangeUser(-7, 8);

 TLatex lt3;

 //lt3.DrawLatexNDC(.25, .8, Form("Background"));

 noiseMap->Draw();
 signalMap->Draw("Same"); 
 //lt3.DrawLatexNDC(.15, .85, Form("#color[4]{Cherenkov} Background"));
 lt3.DrawLatexNDC(.15, .85, Form("#color[4]{Cherenkov}"));
 lt3.DrawLatexNDC(.15, .8, Form("#color[2]{Background}"));




 telMin->Draw(); 
 telMax->Draw(); 
 auto pad = static_cast<TPad*>(signalPNoiseCanvas->cd(2));
 //hTheta2->SetLineColor(kRed);
 hThetaClone->SetLineColor(kRed);//hThetaCh
 
 hThetaCh->SetLineColor(kBlue);
 hThetaCh->SetTitle("Cherenkov Photons And Background");
 hThetaCh->Draw();

 hThetaClone->Draw("Same");
 hThetaClone->SetTitle("Histogram");
 hThetaClone->GetXaxis()->SetRangeUser(0, 0.65);



  TLine* l = new TLine(ckovTrackOut-(Hwidth/2)*0.001, 0, ckovTrackOut-(Hwidth/2)*0.001, up);

  TLine* l2 = new TLine(ckovTrackOut+(Hwidth/2)*0.001, 0, ckovTrackOut+(Hwidth/2)*0.001, up);
  TLine* ll = new TLine(ckovTrackOut-(Hwidth/2)*0.001,  up, ckovTrackOut+(Hwidth/2)*0.001, up);
  l->SetLineColor(kGreen);  l2->SetLineColor(kGreen);ll->SetLineColor(kGreen);
  l->SetLineWidth(2);    l2->SetLineWidth(2);  ll->SetLineWidth(2); 

 l->Draw();  l2->Draw();ll->Draw();

 //box2->Draw();
 //lt2.DrawLatexNDC(.15, .85, Form("#color[4]{Cherenkov Photon}"));
 //lt2.DrawLatexNDC(.15, .8, Form("Background Photon"));
 signalPNoiseCanvas->Show();



/*
 hSignalAndNoiseMap2->SetMarkerStyle(3);
 hSignalAndNoiseMap2->Draw();

*/



 /*
 TCanvas *c1 = new TCanvas("c1","c1",800,800);
 TCanvas *c2 = new TCanvas("c2","c2",800,800);
 TCanvas *c3 = new TCanvas("c3","c3",800,800);
 TCanvas *c4 = new TCanvas("c4","c4",600,400);

 c1->cd();

 c2->cd();
 hTheta2->Draw(); 
 
 c3->cd();
 hThetawg->Draw();
 
 c4->cd();
 NPhoton->Draw(); */ 

 TCanvas* ringCanvas = new TCanvas("Ring Cherenkov", "Ring Cherenkov", 800,800);
 ringCanvas->cd();
 hThetaRing->Draw();


  gStyle->SetOptStat("ei");
 /*
  gStyle->SetOptStat("erm");
 {
   TCanvas* testCanvas = new TCanvas("test" ,"test", 800,800);
   testCanvas->Divide(2,1);
   auto pad = static_cast<TPad *>(testCanvas->cd(1));
   //pad->SetLogy(1);
   hThetawg->Draw();

   testCanvas->cd(2);
   auto pad2 = static_cast<TPad *>(testCanvas->cd(2));
   hphotonMap->Draw();

 } */ 


/*
 {

    TCanvas *houghCanvas = new TCanvas("Hough Canvas","Hough Canvas", 800, 800);
  houghCanvas->Divide(2,2);

   houghCanvas->cd(1);
   hTheta->Draw();

   houghCanvas->cd(2);
   phots->Draw();

   houghCanvas->cd(3);
   photsw->Draw();

   houghCanvas->cd(4);
   resultw->Draw();

   houghCanvas->Show();
   houghCanvas->SaveAs("houghCanvas.png");
 }*/

  
 /*
 {
   TCanvas* testCanvas2 = new TCanvas("test2" ,"test2", 1600,800);
   testCanvas2->Divide(2,2);
   testCanvas2->cd(1);
   auto hTheta_cl = static_cast<TH1F*>(hTheta->Clone());
   hTheta_cl->Draw();
   hTheta_cl->SetTitle("a) Background Angle Distribution;angle [rad]; counts/1 mrad");
   hTheta_cl->SetTitleSize(hTheta_cl->GetTitleSize("x")*1.4, "xy");
   hTheta_cl->SetLabelSize(hTheta_cl->GetLabelSize("x")*1.25, "xy");

   auto pad2 = static_cast<TPad *>(testCanvas2->cd(2));
   pad2->SetLogy(1);
   auto hTheta2Log = static_cast<TH1F*>(hTheta2->Clone());
   hTheta2Log->SetTitle("b) LogY Photon Cherenkov Angle Distribution;angle [rad]; counts/1 mrad");
   hTheta2Log->SetTitleSize(hTheta_cl->GetTitleSize("x"), "xy");
   hTheta2Log->SetLabelSize(hTheta_cl->GetLabelSize("x"), "xy");
   hTheta2Log->Draw();


   auto pad3 = static_cast<TPad *>(testCanvas2->cd(3));
   //pad->SetLogy(1);
   auto hTheta2_cl = static_cast<TH1F*>(hTheta2->Clone());
   hTheta2_cl->SetTitleSize(hTheta_cl->GetTitleSize("x"), "xy");
   hTheta2_cl->SetLabelSize(hTheta_cl->GetLabelSize("x"), "xy");
   hTheta2_cl->Draw();
   hTheta2_cl->SetTitle("c) Photon Cherenkov Angle Distribution;angle [rad]; counts/1 mrad");
   auto pad4 = static_cast<TPad *>(testCanvas2->cd(4));
   //pad->SetLogy(1);
   auto hThetawg_cl = static_cast<TH1F*>(hThetawg->Clone());
   hThetawg_cl->SetTitleSize(hTheta_cl->GetTitleSize("x"), "xy");
   hThetawg_cl->SetLabelSize(hTheta_cl->GetLabelSize("x"), "xy");
   hThetawg_cl->SetTitle("d) (Weighted) Track Photon Cherenkov Angle Distribution;angle [rad]; counts/1 mrad");
   hThetawg_cl->Draw();
 } */

 
 TFile *outFile = TFile::Open("background.root","RECREATE");
 
 outFile->WriteObject(hThetaRing,"hThetaRing");
 
 outFile->Write();
 outFile->Close();
 
}
//**********************************************************************************************************************************************************************************************************
float GetFreonIndexOfRefraction(float x)
  
{
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




  //TLine* l = new TLine(ckovTrack-0.012, ckovTrack+0.012, 0.1, resultw->GetMaximum(10.));

  //0, kThetaMax
  TLatex lt;
  lt.DrawLatexNDC(.15, .85, Form("Track Cherenkov #Theta_{c} :"));
  lt.DrawLatexNDC(.16, .775, Form("Actual = %.3f", meanCherenkovAngle));
  lt.DrawLatexNDC(.16, .7, Form("Predicted = %.3f", ckovTrack));

  houghCanvas->Show();
  houghCanvas->SaveAs("houghCanvas.png");


  /*
  TCanvas *photsCanvas = new TCanvas("photss","photss", 800, 800);
  photsCanvas->cd(1);
  phots->Draw();
  photsCanvas->SaveAs("phots.png");

  TCanvas *photswCanvas = new TCanvas("photsw","photsw", 800, 800);
  photswCanvas->cd(1);
  photsw->Draw();
  photswCanvas->SaveAs("photsw.png");

  TCanvas *resultwCanvas = new TCanvas("resultw","resultw", 800, 800);
  resultwCanvas->cd(1);
  resultw->Draw();
  resultwCanvas->SaveAs("resultw.png");


  //  
  */
   
  Printf("ckovTrack = %f", ckovTrack);
  return ckovTrack;

}



/*
double calcRingGeom(double ckovAng, int level)
{
  // Find area covered in the PC acceptance
  // Arguments: ckovAng - cerenkov angle
  //            level   - precision in finding area and portion of ring accepted (multiple of 50)
  //   Returns: area of the ring in cm^2 for given theta ckov
f
  int kN = 50 * level;
  int nPoints = 0;
  double area = 0;

  bool first = kFALSE;

  // this needs to be changed?
  // TVector2
  o2::math_utils::Vector2D<double> pos1;

  for (int i = 0; i < kN; i++) {
    if (!first) {
      pos1 = o2::hmpid::Recon::tracePhot(ckovAng, double(TMath::TwoPi() * (i + 1) / kN)); // find a good trace for the first photon
      if (pos1.X() == -999)
        continue; // no area: open ring

      if (!fParam->isInside(pos1.X(), pos1.Y(), 0)) {
        pos1 = o2::hmpid::Recon::intWithEdge(fMipPos, pos1); // find the very first intersection...
      } else {
        if (!Param::isInDead(1.0f, 1.0f)) // ef : moved method from Param.cxx to h
          nPoints++;                      // photon is accepted if not in dead zone
      }
      first = kTRUE;
      continue;
    }
    o2::math_utils::Vector2D<double> pos2 = o2::hmpid::Recon::tracePhot(ckovAng, double(TMath::TwoPi() * (i + 1) / kN)); // trace the next photon
    if (pos2.X() == -999)
      continue; // no area: open ring
    if (!fParam->isInside(pos2.X(), pos2.Y(), 0)) {
      pos2 = o2::hmpid::Recon::intWithEdge(fMipPos, pos2);
    } else {
      if (!Param::isInDead(pos2.X(), pos2.Y()))
        nPoints++; // photon is accepted if not in dead zone
    }

    area += TMath::Abs((pos1 - fMipPos).X() * (pos2 - fMipPos).Y() - (pos1 - fMipPos).Y() * (pos2 - fMipPos).X()); // add area of the triangle...
    pos1 = pos2;
  }
  //---  find area and length of the ring;
  fRingAcc = (double)nPoints / (double)kN;
  area *= 0.5;
  fRingArea = area;
  return fRingArea;
} // FindRingGeom()


void trs2Lors(o2::math_utils::Vector3D<double> dirCkov, double& thetaCer, double& phiCer) const
{
  // Theta Cerenkov reconstruction
  //  Arguments: dirCkov photon vector in TRS
  //    Returns: thetaCer of photon in LORS
  //               phiCer of photon in LORS

  // TRotation mtheta;
  // mtheta.RotateY(fTrkDir.Theta()); ef : changed to :

  ROOT::Math::Rotation3D mtheta(ROOT::Math::RotationY(fTrkDir.Theta()));

  // TRotation mphi;
  // mphi.RotateZ(fTrkDir.Phi()); ef : changed to :

  ROOT::Math::Rotation3D mphi(ROOT::Math::RotationZ(fTrkDir.Phi()));

  ROOT::Math::Rotation3D mrot = mphi * mtheta;

  // ef : TVector3->Polar3D
  o2::math_utils::Vector3D<double> dirCkovLORS;
  dirCkovLORS = mrot * dirCkov;
  phiCer = dirCkovLORS.Phi();     // actual value of the phi of the photon
  thetaCer = dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon
}



template <typename T> // typename
o2::math_utils::Vector2D<T> tracePhot(double ckovThe, double ckovPhi) const
{
  // Trace a single Ckov photon from emission point somewhere in radiator up to photocathode taking into account ref indexes of materials it travereses
  // Arguments: ckovThe,ckovPhi- photon ckov angles in TRS, [rad]
  //   Returns: distance between photon point on PC and track projection

  double theta, phi;
  o2::math_utils::Vector3D<double> dirTRS; // ef TVector3 -> Polar3D

  Polar3DVector dirLORS;

  // ef SetMagThetaPhi->SetCoordinates

  dirTRS.SetCoordinates(1, ckovThe, ckovPhi); // photon in TRS
  trs2Lors(dirTRS, theta, phi);
  dirLORS.SetCoordinates(1, theta, phi); // photon in LORS
  return traceForward(dirLORS);          // now foward tracing
} // TracePhot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T>
void propagate(const Polar3DVector& dir, o2::math_utils::Vector3D<double>& pos, double z) const
{
  // Finds an intersection point between a line and XY plane shifted along Z.
  // Arguments:  dir,pos   - vector along the line and any point of the line
  //             z         - z coordinate of plain
  //   Returns:  none
  //   On exit:  pos is the position if this intesection if any
  static o2::math_utils::Vector3D<double> nrm(0, 0, 1);
  o2::math_utils::Vector3D<double> pnt(0, 0, z);

  o2::math_utils::Vector3D<double> diff = pnt - pos;
  double sint = 0; //(nrm * diff) / (nrm * dir);
  pos += sint * dir;
} // Propagate()


template <typename T>
const o2::math_utils::Vector2D<T> intWithEdge(o2::math_utils::Vector2D<T> p1, o2::math_utils::Vector2D<T> p2)
{
  // It finds the intersection of the line for 2 points traced as photons
  // and the edge of a given PC
  // Arguments: 2 points obtained tracing the photons
  //   Returns: intersection point with detector (PC) edges

  double xmin = (p1.X() < p2.X()) ? p1.X() : p2.X();
  double xmax = (p1.X() < p2.X()) ? p2.X() : p1.X();
  double ymin = (p1.Y() < p2.Y()) ? p1.Y() : p2.Y();
  double ymax = (p1.Y() < p2.Y()) ? p2.Y() : p1.Y();

  double m = TMath::Tan((p2 - p1).Phi());
  o2::math_utils::Vector2D<double> pint;
  // intersection with low  X
  pint.SetCoordinates((double)(p1.X() + (0 - p1.Y()) / m), 0.);
  if (pint.X() >= 0 && pint.X() <= fParam->sizeAllX() &&
      pint.X() >= xmin && pint.X() <= xmax &&
      pint.Y() >= ymin && pint.Y() <= ymax)
    return pint;
  // intersection with high X
  pint.SetCoordinates((double)(p1.X() + (fParam->sizeAllY() - p1.Y()) / m), (double)(fParam->sizeAllY()));
  if (pint.X() >= 0 && pint.X() <= fParam->sizeAllX() &&
      pint.X() >= xmin && pint.X() <= xmax &&
      pint.Y() >= ymin && pint.Y() <= ymax)
    return pint;
  // intersection with left Y
  pint.SetCoordinates(0., (double)(p1.Y() + m * (0 - p1.X())));
  if (pint.Y() >= 0 && pint.Y() <= fParam->sizeAllY() &&
      pint.Y() >= ymin && pint.Y() <= ymax &&
      pint.X() >= xmin && pint.X() <= xmax)
    return pint;
  // intersection with righ Y
  pint.SetCoordinates((double)(fParam->sizeAllX()), (double)(p1.Y() + m * (fParam->sizeAllX() - p1.X()))); // ef: Set->SetCoordinates
  if (pint.Y() >= 0 && pint.Y() <= fParam->sizeAllY() &&
      pint.Y() >= ymin && pint.Y() <= ymax &&
      pint.X() >= xmin && pint.X() <= xmax)
    return pint;
  return p1;
} // IntWithEdge()v */

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
    if (binEnt > max) max = binEnt;
    Printf("ent i %d || val %f ", i, binEnt);
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
    if (binEnt > max) max = binent;
    Printf("ent i %d || val %f || val2 %f", i, binEnt, binent);

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
  const auto cos_ckov = TMath::Sqrt(p_sq + m*m)/(cos_ckov_denom);

  const auto ckovAngle = TMath::ACos(cos_ckov);

  return ckovAngle;
}


static constexpr double arrWaveLenDefault[30] = {
  162, 164, 166, 168, 170, 172, 174, 176, 178, 180,
  182, 184, 186, 188, 190, 192, 194, 196, 198, 200,
  202, 204, 206, 208, 210, 212, 214, 216, 218, 220};
static constexpr double nm2eV = 1239.842609;

double randomMass() 
{  
  auto index = static_cast<int>(rndInt->Integer(3));
  Printf("randomMass indes = %d", index);
  return masses[index];
}

double randomEnergy()
{
  
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

