#define BCR_CalibTools_cxx
#include "BCR_CalibTools.h"


// ===========================================================================
// SECTION 1: DEAD Material Corrections
// ===========================================================================
Double_t BCR_CalibTools::DeadCor2(Double_t Ein, Double_t x0) {
  
  // Dano 2007/11/15 Ein = Emeasured
  
  if (x0<0.) x0 = 0.;

  Int_t idx = (Int_t) x0;
  if (idx > 2) idx=2;
  
  Double_t f1 = 1. + a[idx  ]/TMath::Power(Ein,b[idx  ]);
  Double_t f2 = 1. + a[idx+1]/TMath::Power(Ein,b[idx+1]);
  
  Double_t f = f1 + (f2-f1)*(x0-idx);

  return f;

}

Double_t BCR_CalibTools::EcorDeadCor2(Double_t Ein, Double_t x0) {
  return Ein/DeadCor2(Ein,x0);
}
  
Double_t BCR_CalibTools::EexpDeadCor2(Double_t Etrue, Double_t x0, Int_t nIter) {
  
  Double_t Eexp = Etrue;
  Double_t fac = DeadCor2(Etrue,x0);
  for (int ll=0; ll<nIter; ll++) {
    Eexp = fac * Etrue;
    fac  = DeadCor2(Eexp,x0);
  }
  return Eexp;
}



// ===========================================================================
// SECTION 2: Geometry in RCAL
// ===========================================================================

Bool_t BCR_CalibTools::inBeamHole(Double_t xCell, Double_t yCell){
  Bool_t rc = kFALSE;
  if ((xCell >=  xedge[11]   &&  xCell <= xedge[12])
      &&
      (yCell >= yx0edge[20]  && yCell <= yx0edge[21])
      )
    rc = kTRUE;
  return rc;
}  

Double_t BCR_CalibTools::dxCellEdge(Double_t xCell, Double_t yCell){
  
  Double_t rc  = 999.;
  Double_t rcl = 998.;
  Double_t rcr = 997.;
  
  if (inBeamHole(xCell,yCell))  return rc;
  
  if (xCell >=  xedge[11]   &&  xCell <= xedge[12]) {    
    if (yCell>0) {
      rcl = xCell-xedge[11];
      rcr = xCell-(xedge[12]-0.5);      
      if (rcl > -rcr) rc = rcr;
      else            rc = rcl;
    } else {
      rcl = xCell-(xedge[11]+0.5);
      rcr = xCell-(xedge[12]);      
      if (rcl > -rcr) rc = rcr;
      else            rc = rcl;
    }
    return rc;      
  } else if (xCell <  xedge[11]) {
    for (int ll=11; ll>0; ll--) {
      if (xCell >=xedge[ll-1] && xCell<xedge[ll]){
	rcl = xCell-(xedge[ll-1]);
	rcr = xCell-(xedge[ll]);            
	if (rcl > -rcr) rc = rcr;
	else            rc = rcl;
	return rc;
      }
    }
  } else if (xCell >  xedge[12]) {
    for (int ll=12; ll<23; ll++) {
      if (xCell >=xedge[ll] && xCell<xedge[ll+1]){
	rcl = xCell-(xedge[ll]);
	rcr = xCell-(xedge[ll+1]);            
	if (rcl > -rcr) rc = rcr;
	else            rc = rcl;
	return rc;
      }
    }    
  }
  return rc;
}


Int_t BCR_CalibTools::RcalCellXbin(Double_t xCell) {
  Int_t rc=-1;
  if (xCell>=xedge[11]) {
    for (int ll=11; ll<23; ll++) {
      if (xCell>=xedge[ll] && xCell<xedge[ll+1]) return ll;
    }
  } else {
    for (int ll=11; ll>0; ll--) {
      if (xCell>=xedge[ll-1] && xCell<xedge[ll]) return ll-1;
    }
  }
  return rc;
}

Int_t BCR_CalibTools::RcalCellYbinHiGran(Double_t yCell) {
  Int_t rc=-1;
  if (yCell>=yhigran[23]) {
    for (int ll=23; ll<46; ll++) {
      if (yCell>=yhigran[ll] && yCell<yhigran[ll+1]) return ll;
    }
  } else {
    for (int ll=23; ll>0; ll--) {
      if (yCell>=yhigran[ll-1] && yCell<yhigran[ll+1]) return ll-1;
    }
  }
  return rc;
}

Int_t BCR_CalibTools::RcalCellYbinX0Gran(Double_t yCell) {
  Int_t rc=-1;
  if (yCell>= yx0edge[21]) {
    for (int ll=21; ll<30; ll++) {
      if (yCell>= yx0edge[ll]  && yCell< yx0edge[ll+1]) return ll;
    }
  } else {
    for (int ll=21; ll>0; ll--) {
      if (yCell>= yx0edge[ll-1] && yCell<yx0edge[ll]) return ll-1;
    }
  }
  return rc;
}



Double_t BCR_CalibTools::dyCellEdge(Double_t xCell, Double_t yCell){
  
  Double_t rc  = 999.;
  Double_t rcl = 998.;
  Double_t rch = 997.;
  
  if (inBeamHole(xCell,yCell))  return rc;
  


  if (xCell >=  xedge[11]   &&  xCell <= xedge[12]) {    
    if (yCell>yx0edge[21]) {
      for (int ll=21; ll<30; ll++) {
	if (yCell>=yx0edge[ll] && yCell<=yx0edge[ll+1] ) {
	  rch = yCell-yx0edge[ll+1];
	  rcl = yCell-yx0edge[ll];
	  if (rcl > -rch) rc = rch;
	  else            rc = rcl;
	  return rc;
	}
      }
    } else if (yCell<yx0edge[20]) {
      for (int ll=20; ll>0; ll--) {
	if (yCell>=yx0edge[ll-1] && yCell<=yx0edge[ll] ) {
	  rch = yCell-yx0edge[ll];
	  rcl = yCell-yx0edge[ll-1];
	  if (rcl > -rch) rc = rch;
	  else            rc = rcl;
	  return rc;
	}
      }
    }
  } else {
    Int_t iBinX = RcalCellXbin(xCell);
    if (yCell>=ylogran[yhigranbounds[iBinX][1]]) {
      for (int ll = yhigranbounds[iBinX][1]; ll<23; ll++) {
	if (yCell >= ylogran[ll] && yCell < ylogran[ll+1] ) {
	  rch = yCell-ylogran[ll+1];
	  rcl = yCell-ylogran[ll];
	  if (rcl > -rch) rc = rch;
	  else            rc = rcl;
	  return rc;
	} 
      }
    } else if (yCell<=ylogran[yhigranbounds[iBinX][0]]) {
      for (int ll = yhigranbounds[iBinX][0]; ll>0; ll--) {
	if (yCell >= ylogran[ll-1] && yCell < ylogran[ll] ) {
	  rch = yCell-ylogran[ll];
	  rcl = yCell-ylogran[ll-1];
	  if (rcl > -rch) rc = rch;
	  else            rc = rcl;
	  return rc;
	} 
      }
    } else {
      Int_t iBinY = RcalCellYbinHiGran(yCell);
      rch = yCell-yhigran[iBinY+1];
      rcl = yCell-yhigran[iBinY];
      if (rcl > -rch) rc = rch;
      else            rc = rcl;
      return rc;
    }
    
  }
  return rc;
}


// ===========================================================================
// SECTION 3: Fiducial Volumes
// ===========================================================================
Bool_t  BCR_CalibTools::fiducialSRTD(Double_t xpos, Double_t ypos)
{
  Bool_t rc = kFALSE;
 
  if (ypos > 0 ) {
    if (xpos < 10. && xpos > -35.5  
	&& ypos<37) rc = kTRUE;
    if (xpos > 10. && xpos < 34.5  
	&& ypos<33) rc = kTRUE;
  } else {
    if (xpos < -10. && xpos > -35.5  
	&& ypos>-33) rc = kTRUE;
    if (xpos > -10. && xpos < 34.5  
	&& ypos>-37) rc = kTRUE;
  }  
  return rc;
}



// ===========================================================================
// SECTION 4: Non uniformity corrections
// ===========================================================================
Double_t  BCR_CalibTools::gaps_Data_xtop(Double_t pos)
{
  Double_t rc=1.;  
  if (pos <-80. || pos >80.) return rc;
  Int_t ix =((Int_t) pos + 80)/20;
  const Double_t * par = PARxtop[ix];
  Double_t scale = 1.006; // due to extraction of non-unif in iter N-1
  
  rc = 1./scale * (par[0]
		   +par[1]*TMath::Exp(-(pos-par[2])*(pos-par[2])/
				      par[3]/par[3]) 
		   +par[4]*TMath::Exp(-(pos-par[5])*(pos-par[5])/
				      par[6]/par[6])); 
  return rc;
}

Double_t  BCR_CalibTools::gaps_Data_xbot(Double_t pos)
{
  Double_t rc=1.;
  if (pos <-80. || pos >80.) return rc;
  Int_t ix =((Int_t) pos + 80)/20;
  const Double_t * par = PARxbot[ix];
  Double_t scale = 1.006; // due to extraction of non-unif in iter N-1
  
  rc = 1./scale * (par[0]
		   +par[1]*TMath::Exp(-(pos-par[2])*(pos-par[2])/
				      par[3]/par[3]) 
		   +par[4]*TMath::Exp(-(pos-par[5])*(pos-par[5])/
				      par[6]/par[6])); 
  return rc;
}
Double_t  BCR_CalibTools::gaps_Data_yleft(Double_t pos)
{
  Double_t rc=1.;
  if (pos <-85. || pos >85.) return rc;
  Int_t ix =((Int_t) pos + 85)/10;
  const Double_t * par = PARyleft[ix];
  Double_t scale = 1.006; // due to extraction of non-unif in iter N-1
  
  rc = 1./scale * (par[0]
		   +par[1]*TMath::Exp(-(pos-par[2])*(pos-par[2])/
				      par[3]/par[3])); 
  return rc;

}
Double_t  BCR_CalibTools::gaps_Data_yright(Double_t pos)
{
  Double_t rc=1.;
  if (pos <-85. || pos >85.) return rc;
  Int_t ix =((Int_t) pos + 85)/10;
  const Double_t * par = PARyright[ix];
  Double_t scale = 1.006; // due to extraction of non-unif in iter N-1
  
  rc = 1./scale * (par[0]
		   +par[1]*TMath::Exp(-(pos-par[2])*(pos-par[2])/
				      par[3]/par[3])); 
  return rc;
}
Double_t  BCR_CalibTools::gaps_Data_ymid(Double_t pos)
{
  Double_t rc=1.0; // no reasonable parameterization possible
  if (RcalCellXbin(pos) == 11) rc =1.0;
  return rc;
}


Double_t  BCR_CalibTools::gaps_MC_xtop(Double_t pos)
{
  Double_t rc=1.;
  if (pos <-80. || pos >80.) return rc;
  Int_t ix =((Int_t) pos + 80)/20;
  const Double_t * par = PMCxtop[ix];
  
  rc = par[0]
    +par[1]*TMath::Exp(-(pos-par[2])*(pos-par[2])/
		       par[3]/par[3]) 
    +par[4]*TMath::Exp(-(pos-par[5])*(pos-par[5])/
		       par[6]/par[6]); 
  return rc;
}
Double_t  BCR_CalibTools::gaps_MC_xbot(Double_t pos)
{
  Double_t rc=1.;
  if (pos <-80. || pos >80.) return rc;
  Int_t ix =((Int_t) pos + 60)/20;
  const Double_t * par = PMCxbot[ix];
  
  rc = par[0]
    +par[1]*TMath::Exp(-(pos-par[2])*(pos-par[2])/
		       par[3]/par[3]) 
    +par[4]*TMath::Exp(-(pos-par[5])*(pos-par[5])/
		       par[6]/par[6]); 
  return rc;
}
Double_t  BCR_CalibTools::gaps_MC_yleft(Double_t pos)
{
  Double_t rc=1.01; //average bias
  if (pos <-55. || pos >55.) return rc;
  Int_t ix =((Int_t) pos + 55)/10;
  const Double_t * par = PMCyleft[ix];
  
  rc = par[0]
    +par[1]*TMath::Exp(-(pos-par[2])*(pos-par[2])/
		       par[3]/par[3]); 
  return rc;
}
Double_t  BCR_CalibTools::gaps_MC_yright(Double_t pos)
{
  Double_t rc=1.01; //average bias
  if (pos <-65. || pos >65.) return rc;
  Int_t ix =((Int_t) pos + 65)/10;
  const Double_t * par = PMCyright[ix];
  
  rc = par[0]
    +par[1]*TMath::Exp(-(pos-par[2])*(pos-par[2])/
		       par[3]/par[3]); 
  return rc;
}

Double_t  BCR_CalibTools::gaps_MC_ymid(Double_t pos)
{
  Double_t rc=1.0; // no reasonable parameterization possible
  if (RcalCellXbin(pos) == 11) rc =1.0;
  return rc;
}


Double_t  BCR_CalibTools::nonunif(Double_t xpos, Double_t ypos,
				  Int_t MCflag)
{
  Double_t rc   = 1.;
  Double_t xfac = 1.;
  Double_t yfac = 1.;

  Int_t iBinX = RcalCellXbin(xpos); 

  if (MCflag == 0) { // Data

    if (ypos>0) xfac = gaps_Data_xtop(xpos);
    else        xfac = gaps_Data_xbot(xpos);

    if      (iBinX < 11) yfac = gaps_Data_yleft(ypos);
    else if (iBinX > 11) yfac = gaps_Data_yright(ypos);
    else                 yfac = gaps_Data_ymid(ypos);

  } else { // MC

    if (ypos>0) xfac = gaps_MC_xtop(xpos);
    else        xfac = gaps_MC_xbot(xpos);

    if      (iBinX < 11) yfac = gaps_Data_yleft(ypos);
    else if (iBinX > 11) yfac = gaps_Data_yright(ypos);
    else                 yfac = gaps_Data_ymid(ypos);

  }


  rc = xfac *yfac;
  if (rc<0) rc =1.;
  else rc = 1/rc;
  return rc;
}

 
Double_t BCR_CalibTools::gaps_Data_xtop(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_Data_xtop(pos);
}
Double_t BCR_CalibTools::gaps_Data_xbot(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_Data_xbot(pos);
}
Double_t BCR_CalibTools::gaps_Data_yleft(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_Data_yleft(pos);
}
Double_t BCR_CalibTools::gaps_Data_yright(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_Data_yright(pos);
}
Double_t BCR_CalibTools::gaps_Data_ymid(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_Data_ymid(pos);
}

Double_t BCR_CalibTools::gaps_MC_xtop(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_MC_xtop(pos);
}
Double_t BCR_CalibTools::gaps_MC_xbot(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_MC_xbot(pos);
}
Double_t BCR_CalibTools::gaps_MC_yleft(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_MC_yleft(pos);
}
Double_t BCR_CalibTools::gaps_MC_yright(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_MC_yright(pos);
}
Double_t BCR_CalibTools::gaps_MC_ymid(Double_t * xxx, Double_t * par)
{
  Double_t pos = xxx[0]-par[0]; 
  return gaps_MC_ymid(pos);
}

Double_t BCR_CalibTools::ecorr_newscheme(Double_t Einraw, Double_t x0,
					 Double_t xcal, Double_t ycal,
					 Int_t MCflag)
{
  // Einraw:   assumes Escale03 factors have been applied already
  // x0:       number of radiation length passed by the electron candidates
  // xcal,     position of the candidate at the front face of the calorimter
  //  & ycal:    using cal stand alone co-ordinates 
  // MCflag:  =0 Data   !=0 MC          
  Double_t Eunif  = Einraw * nonunif(xcal, ycal, MCflag);
  Double_t Edmcor = EcorDeadCor2(Eunif,x0);
  Double_t Epres  = Edmcor;  // presampler not available (yet)

  return Epres;

}
