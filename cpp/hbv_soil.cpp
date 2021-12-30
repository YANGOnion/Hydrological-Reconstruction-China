#include "Rcpp.h"
using namespace Rcpp;
// [[Rcpp::export]]
List updateState(double Snowpack0,double Meltwater0,double Si0,double Su0p,
                 double P,double Tair,double PET,double LAI,
                 double TP,double TM,double SCF,double CFX,double CFR,double CWH,
                 double mc,double FC,double beta,double LP){
  // Rain/Snow/Melt/Refreezing/interception
  double Su0;
  double Rain;
  double Snow;
  double Melt;
  double Refreezing;
  double Tosoil;
  double Simax;
  double Ptf;
  double Ei;
  double Ea;
  double Recharge;
  if(Tair>TP){
    Rain=P;
    Simax=mc*LAI;
    if(Simax==0){
      Ptf=P;
      Si0=0;
      Ei=0;
    }else{
      Ptf=max(NumericVector::create(0,P-(Simax-Si0)));
      Si0=Si0+P-Ptf;
      Ei=min(NumericVector::create(Si0,PET*pow(Si0/Simax,2/3)));
      Si0=Si0-Ei;
    }
    Snow=0;
  }else{
    Rain=0;
    Ptf=0;
    Snow=P*SCF;
    Ei=0;
  }
  if(Tair>TM){
    Melt=min(NumericVector::create(Snowpack0,CFX*(Tair-TM)));
    Refreezing=0;
  }else{
    Melt=0;
    Refreezing=min(NumericVector::create(Meltwater0,CFR*CFX*(TM-Tair)));
  }
  // update snow storages
  Snowpack0=Snowpack0-Melt+Snow+Refreezing;
  Meltwater0=Meltwater0-Refreezing+Melt;
  Tosoil=max(NumericVector::create(0,Meltwater0-CWH*Snowpack0));
  Meltwater0=Meltwater0-Tosoil;
  // update soil storage
  Su0=Su0p*FC;
  Recharge=(Ptf+Tosoil)*pow(Su0/FC,beta);
  Su0=min(NumericVector::create(FC,Su0+Ptf+Tosoil-Recharge));
  Ea=(PET-Ei)*min(NumericVector::create(1,Su0/LP/FC));
  Ea=min(NumericVector::create(Su0,Ea));
  Su0=Su0-Ea;
  Su0p=Su0/FC;
  return List::create(_["Snowpack0"]=Snowpack0,_["Meltwater0"]=Meltwater0,_["Si0"]=Si0,_["Su0p"]=Su0p,
                      _["Rain"]=Rain,_["Snow"]=Snow,_["Melt"]=Melt,_["Tosoil"]=Tosoil,
                        _["Ei"]=Ei,_["Ea"]=Ea);
}

// [[Rcpp::export]]
DataFrame runModel(DataFrame df,double I_Snowpack0,double I_Meltwater0,double I_Si0,double I_Su0p,
                   double P_TP,double P_TM,double P_SCF,double P_CFX,double P_CFR,double P_CWH,
                   double P_mc,double P_FC,double P_beta,double P_LP){
  // V_: vector, P_: parameter, I_: initial value, S_: state
  // output data
  int len=df.nrows();
  List output;
  NumericVector V_Snowpack0 (len+1);
  NumericVector V_Meltwater0 (len+1);
  NumericVector V_Si0 (len+1);
  NumericVector V_Su0p (len+1);
  NumericVector V_Rain (len);
  NumericVector V_Snow (len);
  NumericVector V_Melt (len);
  NumericVector V_Tosoil (len);
  NumericVector V_Ei (len);
  NumericVector V_Ea (len);
  // input data
  NumericVector V_P=df["prcp"];
  NumericVector V_PET=df["pet"];
  NumericVector V_T=df["tair"];
  NumericVector V_LAI=df["lai"];
  // initial states
  V_Snowpack0(0)=I_Snowpack0;
  V_Meltwater0(0)=I_Meltwater0;
  V_Si0(0)=I_Si0;
  V_Su0p(0)=I_Su0p;
  // update
  for(int i=0;i<len;i++){
    output=updateState(V_Snowpack0(i),V_Meltwater0(i),V_Si0(i),V_Su0p(i),
                       V_P(i),V_T(i),V_PET(i),V_LAI(i),
                       P_TP, P_TM, P_SCF, P_CFX, P_CFR, P_CWH,
                       P_mc, P_FC, P_beta, P_LP);
    V_Snowpack0(i+1)=output["Snowpack0"];
    V_Meltwater0(i+1)=output["Meltwater0"];
    V_Si0(i+1)=output["Si0"];
    V_Su0p(i+1)=output["Su0p"];
    V_Rain(i)=output["Rain"];
    V_Snow(i)=output["Snow"];
    V_Melt(i)=output["Melt"];
    V_Tosoil(i)=output["Tosoil"];
    V_Ei(i)=output["Ei"];
    V_Ea(i)=output["Ea"];
  }
  V_Snowpack0.erase(0);
  V_Meltwater0.erase(0);
  V_Si0.erase(0);
  V_Su0p.erase(0);
  DataFrame df_out=DataFrame::create(_["Snowpack"]=V_Snowpack0,_["Meltwater"]=V_Meltwater0,_["Si"]=V_Si0,_["Sup"]=V_Su0p,
                                     _["Rain"]=V_Rain,_["Snow"]=V_Snow,_["Melt"]=V_Melt,_["Tosoil"]=V_Tosoil,
                                     _["Ei"]=V_Ei,_["Ea"]=V_Ea);
  return df_out;
}

