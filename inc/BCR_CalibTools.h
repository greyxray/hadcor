#ifndef BCR_CalibTools_h
#define BCR_CalibTools_h


#include <iostream>
#include <fstream>
#include <iomanip>

#include "TMath.h"



class BCR_CalibTools {

 private:
  // parameters for dead material correction
  static const Double_t a[4];
  static const Double_t b[4];

  // geometry of cracks in rcal
  static const Double_t xedge[24];
  static const Double_t ylogran[24];
  static const Double_t yhigran[47];
  static const Int_t    yhigranbounds[23][2];
  static const Double_t yx0edge[31];


  // parametization of non-uniformities
  static const Double_t PARxtop[8][7];    // done
  static const Double_t PARxbot[8][7];    // done
  static const Double_t PARyleft[17][4];  // done
  //  static const Double_t PARymid[8][7];   
  static const Double_t PARyright[17][4]; // done

  static const Double_t PMCxtop[8][7];    // done
  static const Double_t PMCxbot[8][7];    // done
  static const Double_t PMCyleft[11][4];   // done
  //  static const Double_t PMCymid[8][7];
  static const Double_t PMCyright[13][4];  // done


 public:

  // deadmaterial parameterization and correction functions
  static Double_t DeadCor2(Double_t Ein, Double_t x0);
  static Double_t EcorDeadCor2(Double_t Ein, Double_t x0);
  static Double_t EexpDeadCor2(Double_t Etrue, Double_t x0, Int_t nIter=5);

  static Double_t nonunif(Double_t xpos, Double_t ypos,
			  Int_t MCflag=0);

  // parameterizations
  static Double_t gaps_Data_xtop(Double_t xpos);
  static Double_t gaps_Data_xbot(Double_t xpos);
  static Double_t gaps_Data_yleft(Double_t ypos);
  static Double_t gaps_Data_yright(Double_t ypos);
  static Double_t gaps_Data_ymid(Double_t ypos);

  static Double_t gaps_MC_xtop(Double_t xpos);
  static Double_t gaps_MC_xbot(Double_t xpos);
  static Double_t gaps_MC_yleft(Double_t ypos);
  static Double_t gaps_MC_yright(Double_t ypos);
  static Double_t gaps_MC_ymid(Double_t ypos);

  // overloaded functions to interface with TF1
  static Double_t gaps_Data_xtop(Double_t * xxx, Double_t * par);
  static Double_t gaps_Data_xbot(Double_t * xxx, Double_t * par);
  static Double_t gaps_Data_yleft(Double_t * xxx, Double_t * par);
  static Double_t gaps_Data_yright(Double_t * xxx, Double_t * par);
  static Double_t gaps_Data_ymid(Double_t * xxx, Double_t * par);

  static Double_t gaps_MC_xtop(Double_t * xxx, Double_t * par);
  static Double_t gaps_MC_xbot(Double_t * xxx, Double_t * par);
  static Double_t gaps_MC_yleft(Double_t * xxx, Double_t * par);
  static Double_t gaps_MC_yright(Double_t * xxx, Double_t * par);
  static Double_t gaps_MC_ymid(Double_t * xxx, Double_t * par);

  static Double_t ecorr_newscheme(Double_t Einraw, Double_t x0,
				  Double_t xcal, Double_t ycal,
				  Int_t MCflag=0);

  // geometry in RCAL related functions
  static Bool_t inBeamHole(Double_t xCell, Double_t yCell);
  static Double_t dxCellEdge(Double_t xCell, Double_t yCell);
  static Int_t RcalCellXbin(Double_t xCell);
  static Int_t RcalCellYbinHiGran(Double_t yCell);
  static Int_t RcalCellYbinX0Gran(Double_t yCell);
  static Double_t dyCellEdge(Double_t xCell, Double_t yCell);

  // fiducial volumes
  // SRTD
  static Bool_t fiducialSRTD(Double_t xpos, Double_t ypos);
  


};
#endif


#ifdef BCR_CalibTools_cxx

// parameters for dead material correction
const Double_t BCR_CalibTools::a[4] = {0.,-0.12,-0.30,-0.50};
const Double_t BCR_CalibTools::b[4] = {0., 0.5,  0.48, 0.38};


const Double_t BCR_CalibTools::xedge[24] = 
  {-234.945,-214.615,-194.285,-173.955,-153.625,-133.295,
   -112.965, -92.635, -72.305, -51.975, -31.645, -11.315,
   9.515,  29.845,  50.175,  70.505,  90.835, 111.165,
   131.495, 151.825, 172.155, 192.485, 212.815, 233.145};

const Double_t BCR_CalibTools::ylogran[24] = 
  {-230.100,-210.100,-190.100,-170.100,-150.100,-130.100,
   -110.100, -90.100, -70.100, -50.100, -30.100, -10.100,   
   9.900,  29.900,  49.900,  69.900,  89.900, 109.900, 
   129.900, 149.900, 169.900, 189.900, 209.900, 229.900};

const Double_t BCR_CalibTools::yhigran[47] = 
  {-230.100,-220.100,-210.100,-200.100,-190.100,-180.100,
   -170.100,-160.100,-150.100,-140.100,-130.100,-120.100,
   -110.100,-100.100, -90.100, -80.100, -70.100, -60.100, 
   -50.100, -40.100, -30.100, -20.100, -10.100,  -0.100,
   9.900,  19.900,  29.900,  39.900,  49.900,  59.900,
   69.900,  79.900,  89.900,  99.900, 109.900, 119.900, 
   129.900, 139.900, 149.900, 159.900, 169.900, 179.900, 
   189.900, 199.900, 209.900, 219.900, 229.900};

const Int_t    BCR_CalibTools::yhigranbounds[23][2] =
  { {17, 6}, //  0
    {18, 5}, //  1
    {10,13}, //  2
    { 7,16}, //  3 
    { 6,17},  //  4     
    { 5,18},  //  5
    { 4,19},  //  6
    { 3,20},  //  7
    { 3,20},  //  8 
    { 3,20},  //  9
    { 2,21},  // 10
    { 2,17},  // 11
    { 2,21},  // 12
    { 3,20},  // 13
    { 3,20},  // 14
    { 3,20},  // 15
    { 4,19},  // 16
    { 5,18},  // 17
    { 6,17},  // 18
    { 7,16},  // 19
    {10,13},  // 20
    {18,5},  // 21
    {17,5}};  // 22

const Double_t BCR_CalibTools::yx0edge[31] = 
  {-231.80,-211.80,-191.80,-181.80,-171.80,
   -161.80,-151.80,-141.80,-131.80,-121.80,
   -111.80,-101.80, -91.80, -81.80, -71.80,
   -61.80, -51.80, -41.80, -31.80, -21.80, 
   -11.80,
   11.80,  21.80,  31.80,  41.80,  51.80, 
   61.80,  71.80,  81.80,  91.80, 101.80};



  
const Double_t BCR_CalibTools::PARxtop[8][7] = {
  { 1.0042, 0.03041,-72.309,2.56991,0.00000,  0.000,1.00000},
  { 1.0082, 0.02390,-50.970,1.08284,0.00000,  0.000,1.00000},
  { 1.0123, 0.05426,-30.567,1.61250,0.00000,  0.000,1.00000},
  { 1.0024,-0.08997,-11.351,0.42920,0.01542,-12.862,2.05264},
  { 1.0099,-0.29062, 12.077,2.69900,0.27359, 11.993,2.57861},
  { 1.0090, 0.03413, 28.459,2.28075,0.00000,  0.000,1.00000},
  { 1.0063, 0.03273, 48.769,1.75584,0.00000,  0.000,1.00000},
  { 1.0084, 0.02646, 70.363,0.62948,0.00000,  0.000,1.00000}};

const Double_t BCR_CalibTools::PARxbot[8][7] = {
  { 1.0051, 0.01933,-71.317,1.89454,0.00000,  0.000,1.00000},
  { 1.0068, 0.01532,-51.678,1.49334,0.00000,  0.000,1.00000},
  { 1.0038, 0.04931,-30.657,1.51938,0.00000,  0.000,1.00000},
  { 1.0068, 0.02485,-11.390,1.74397,0.00000,  0.000,1.00000},
  { 1.0018,-0.10659,  9.248,0.38522,0.03172,  9.369,1.48736},
  { 1.0065, 0.05503, 28.632,2.13532,0.00000,  0.000,1.00000},
  { 1.0076, 0.01939, 49.342,1.88144,0.00000,  0.000,1.00000},
  { 1.0045, 0.01295, 69.481,1.78747,0.00000,  0.000,1.00000}};


const Double_t BCR_CalibTools::PMCxtop[8][7] = {
  { 1.0000, 0.29974,-67.997,0.00808,0.02029,-67.929,5.61752},
  { 1.0028, 0.02441,-51.808,3.73001,0.00283,-51.805,-2.21385},
  { 1.0049,-0.05603,-30.533,0.48862,0.03101,-30.749,-3.93154},
  { 1.0061,-0.16889,-11.278,0.73852,0.05479,-11.243,2.10461},
  { 1.0042,-0.07920,  6.465,0.33060,0.04024,  8.933,1.10910},
  { 1.0069,-0.05793, 28.724,0.55489,0.02575, 29.676,4.72068},
  { 1.0071,-0.29999, 48.345,2.41681,0.32008, 48.668,-2.56245},
  { 1.0081, 0.02189, 67.000,4.00000,0.00000,  0.000,1.00000}};

  
const Double_t BCR_CalibTools::PMCxbot[8][7] = {
  { 1.0026, 0.02318,-73.000,6.00000,0.00000,  0.000,1.00000},  // OK
  { 1.0022, 0.02262,-50.663,3.87919,0.00000,  0.000,1.00000},  // OK
  { 1.0074, 0.04110,-30.529,2.40906,-0.05743,-30.974,0.66369}, // OK
  { 1.0052, 0.0,    -10.,   0.3,    0.00000,  0.000,1.00000},  // unfitable
  { 1.0050,-0.30000,  9.150,0.55000, 0.0660,   9.330,2.000},    // manually
  { 1.0062, 0.03958, 28.724,3.34354,-0.05890, 29.253,0.66876},  // OK
  { 1.0078, 0.02176, 49.921,5.37736,0.00000,  0.000,1.00000},  // OK
  { 1.0063, 0.01564, 67.405,6.00000,0.00000,  0.000,1.00000}};  // OK





const Double_t BCR_CalibTools::PARyleft[17][4] = {
  { 1.0021, 0.02683,-81.455,5.00000},
  { 1.0140,-0.02987,-70.149,3.04972},
  { 1.0059, 0.02112,-61.271,3.44153},
  { 1.0058,-0.04009,-49.509,2.04531},
  { 1.0126,-0.01956,-39.318,0.64522},
  { 1.0051,-0.04844,-29.355,1.44084},
  { 0.9921, 0.03105,-20.771,5.00000},
  { 1.0055,-0.07584, -9.884,0.83292},
  { 1.0228,-0.05313, -0.186,0.50000},
  { 1.0147,-0.08910,  9.655,1.17097},
  { 1.0218, 0.02006, 18.000,0.90736},
  { 1.0128,-0.05959, 29.156,1.57812},
  { 1.0152, 0.00497, 41.566,0.50000},
  { 1.0074,-0.03064, 48.860,1.81376},
  { 1.0168, 0.01408, 61.289,0.53108},
  { 1.0042,-0.03437, 69.092,2.28947},
  { 0.9962, 0.03263, 81.518,5.00000}};

const Double_t  BCR_CalibTools::PARyright[17][4] = {
  { 1.0109, 0.01761,-81.020,1.62996},
  { 1.0176,-0.03987,-69.204,2.88403},
  { 1.0060, 0.01208,-61.610,1.71291},
  { 1.0101,-0.03409,-49.068,1.60126},
  { 1.0141,-0.00896,-39.798,0.56803},
  { 1.0045,-0.05258,-29.600,1.63741},
  { 1.0135, 0.00178,-18.000,5.00000},
  { 1.0149,-0.08288, -9.845,1.12282},
  { 1.0158,-0.05022, -0.232,0.58159},
  { 1.0063,-0.07067,  9.628,1.13675},
  { 1.0136,-0.02756, 19.086,0.88792},
  { 1.0069,-0.03733, 28.982,1.84618},
  { 1.0133,-0.01916, 38.796,0.99383},
  { 1.0107,-0.03592, 48.598,1.71612},
  { 1.0134,-0.00559, 60.488,0.50000},
  { 1.0204,-0.03875, 70.619,5.00000},
  { 1.0029, 0.01542, 81.925,0.50000}};

const Double_t  BCR_CalibTools::PMCyleft[11][4] = {
  { 1.0067,-0.02588,-48.920,0.50000},
  { 1.0109,-0.03247,-38.840,0.94660},
  { 1.0103,-0.03771,-29.297,1.03686},
  { 1.0121,-0.03313,-19.356,1.12682},
  { 1.0132,-0.05602, -9.918,0.81745},
  { 1.0114,-0.04699, -0.267,0.80050},
  { 1.0124,-0.04514,  9.430,0.97281},
  { 1.0122,-0.03658, 19.105,1.04898},
  { 1.0140,-0.03473, 28.943,1.14280},
  { 1.0121,-0.03083, 38.597,1.03525},
  { 1.0142,-0.03446, 48.572,1.27417}};


const Double_t  BCR_CalibTools::PMCyright[13][4] = {
  { 1.0118,-0.02119,-58.028,1.46502},
  { 1.0117,-0.03061,-49.049,0.91021},
  { 1.0136,-0.03156,-38.987,0.80421},
  { 1.0140,-0.03258,-29.389,1.26435},
  { 1.0175,-0.03465,-19.644,0.90737},
  { 1.0172,-0.05498, -9.913,0.85171},
  { 1.0154,-0.04716, -0.220,1.05345},
  { 1.0155,-0.04826,  9.476,0.86713},
  { 1.0161,-0.03936, 19.261,0.84501},
  { 1.0157,-0.03573, 29.073,0.72712},
  { 1.0118,-0.02865, 38.801,0.99259},
  { 1.0135,-0.02464, 48.446,1.29754},
  { 1.0121,-0.01998, 58.361,0.91396}};



#endif
