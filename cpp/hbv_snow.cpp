#include "Rcpp.h"
using namespace Rcpp;
// [[Rcpp::export]]
List updateState(double Snowpack0,double Meltwater0,
                 double P,double Tair,
                 double TP,double TM,double SCF,double CFX,double CFR,double CWH){
  // Rain/Snow/Melt/Refreezing
  double Rain;
  double Snow;
  double Melt;
  double Refreezing;
  double Tosoil;
  if(Tair>TP){
    Rain=P;
    Snow=0;
  }else{
    Rain=0;
    Snow=P*SCF;
  }
  if(Tair>TM){
    Melt=min(NumericVector::create(Snowpack0,CFX*(Tair-TM)));
    Refreezing=0;
  }else{
    Melt=0;
    Refreezing=min(NumericVector::create(Meltwater0,CFR*CFX*(TM-Tair)));
  }
  // update storages
  Snowpack0=Snowpack0-Melt+Snow+Refreezing;
  Meltwater0=Meltwater0-Refreezing+Melt;
  Tosoil=max(NumericVector::create(0,Meltwater0-CWH*Snowpack0));
  Meltwater0=Meltwater0-Tosoil;
  return List::create(_["Snowpack0"]=Snowpack0,_["Meltwater0"]=Meltwater0,
                      _["Rain"]=Rain,_["Melt"]=Melt,_["Refreezing"]=Refreezing,_["Tosoil"]=Tosoil);
}

// [[Rcpp::export]]
DataFrame runModel(DataFrame df,double I_Snowpack0,double I_Meltwater0,
                   double P_TP,double P_TM,double P_SCF,double P_CFX,double P_CFR,double P_CWH){
  // V_: vector, P_: parameter, I_: initial value, S_: state
  // output data
  int len=df.nrows();
  List output;
  NumericVector V_Snowpack0 (len+1);
  NumericVector V_Meltwater0 (len+1);
  NumericVector V_Rain (len);
  NumericVector V_Melt (len);
  NumericVector V_Refreezing (len);
  NumericVector V_Tosoil (len);
  // input data
  NumericVector V_P=df["prcp"];
  NumericVector V_T=df["tair"];
  // initial states
  V_Snowpack0(0)=I_Snowpack0;
  V_Meltwater0(0)=I_Meltwater0;
  // update
  for(int i=0;i<len;i++){
    output=updateState(V_Snowpack0(i),V_Meltwater0(i),
                       V_P(i),V_T(i),
                       P_TP, P_TM, P_SCF, P_CFX, P_CFR, P_CWH);
    V_Snowpack0(i+1)=output["Snowpack0"];
    V_Meltwater0(i+1)=output["Meltwater0"];
    V_Rain(i)=output["Rain"];
    V_Melt(i)=output["Melt"];
    V_Refreezing(i)=output["Refreezing"];
    V_Tosoil(i)=output["Tosoil"];
  }
  V_Snowpack0.erase(0);
  V_Meltwater0.erase(0);
  DataFrame df_out=DataFrame::create(_["Snowpack"]=V_Snowpack0,_["Meltwater"]=V_Meltwater0,
                                     _["Rain"]=V_Rain,_["Melt"]=V_Melt,_["Refreezing"]=V_Refreezing,_["Tosoil"]=V_Tosoil);
  return df_out;
}


