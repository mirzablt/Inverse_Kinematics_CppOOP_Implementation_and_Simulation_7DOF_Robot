
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif


#include <math.h>

#include <MatrixLibrary.cpp>
#include <InverseKinematicsLibrary.cpp>

#define u_width 4
#define y_width 1


void InverseKinematicsCppOOP_Outputs_wrapper(const real_T *Pose,
			const real_T *Init,
			real_T *Config,
			real_T *CriterionValue)
{
/* Kod u narednim linijama je implementiran u bloku S-Function Builder u mapi Outputs.
   Sve linije koda koje prethode ovom komentaru su automatski generirane pritiskom 
   tastera Build u bloku S-Function Builder. 

   U mapi Libraries u bloku S-Function Builder su navedene implemetirane biblioteke:
     -InverseKinematicsLibrary;
     -MatrixLibrary.

  Ulazi bloka su:
    -Pose - zeljena poza robota data u formi matrice homogene tranformacije koja
            odredjuje zeljenu poziciju i orijentaciju vrha manipulatora. Ovu matricu
            generira blok planiranja trajektorije (koji nije obuhvacen ovim radom);
    -Init - inicijalna vrijednost vektora zglobnih varijabli;
  i deklarisani u mapi Data Properties/Input ports bloka S-funkcije.

  Izlazi bloka su:
     -Config         - vrijednost vektora zglobnih varijabli robota koje je rezultat 
                       algoritma inverzne kinematike, a za koju vrh manipulatora zauzima
                       zeljenu poziciju i orijentaciju;
     -CriterionValue - vrujednost funkcije kriterija. */

    using namespace std;

    
    // Broj zglobova robota KUKA iiwa 14R820.
    
    int const robotJointsNum = 7;
    
	              
    /* DH parametri robotskog manipulatora KUKA iiwa 14R820  specificirani od 
       strane proizvodjaca. */
    
   	double const a1 = 0;
	double const a2 = 0;
	double const a3 = 0;
	double const a4 = 0;
	double const a5 = 0;
    double const a6 = 0;
	double const a7 = 0;
	double const alpaha1 =  pi/2.0;
	double const alpaha2 = -pi/2.0;
	double const alpaha3 = -pi/2.0;
	double const alpaha4 =  pi/2.0;
	double const alpaha5 =  pi/2.0;
	double const alpaha6 = -pi/2.0;
    double const alpaha7 =  0.0;
	double const d1 = 0.360;
	double const d2 = 0;
	double const d3 = 0.420;
	double const d4 = 0;
	double const d5 = 0.400;
	double const d6 = 0;
    double const d7 = 0.126;
    
    
    // Matrica DH parametara KUKA iiwa 14R820 robota.
    
    HTMatrix DH(robotJointsNum,4);
	DH.data[0][0] = a1;        DH.data[0][1] = alpaha1;           DH.data[0][2] = d1;  
	DH.data[1][0] = a2;        DH.data[1][1] = alpaha2;           DH.data[1][2] = d2; 
	DH.data[2][0] = a3;        DH.data[2][1] = alpaha3;           DH.data[2][2] = d3;  
	DH.data[3][0] = a4;        DH.data[3][1] = alpaha4;           DH.data[3][2] = d4;  
	DH.data[4][0] = a5;        DH.data[4][1] = alpaha5;           DH.data[4][2] = d5;  
	DH.data[5][0] = a6;        DH.data[5][1] = alpaha6;           DH.data[5][2] = d6;
        DH.data[6][0] = a7;        DH.data[6][1] = alpaha7;           DH.data[6][2] = d7;
 
    
    /* desiredPoseHTM - desired pose homogene tranformation matrix. To je matrica 
       homogene transformacije koja odredjuje zeljenu poziciju vrha manipulatora.
       Parametar "Pose" je deklarisan u mapi "Data Properties/Input Ports" i on je
       ulaz bloka S-funkcije.
       "Pose" je vektor-vrsta kojeg je potrebno prevesti u 4x4 matricu "DesiredPoseHTM". */
    
    HTMatrix desiredPoseHTM(4, 4);
    
    int k = 0;
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            desiredPoseHTM.data[i][j] = Pose[i+j+k];
        }
        k += 3;
    };
  
    
   /* startPoint je pocetna iteracija algoritma optimizacije. U svakom koraku 
       simulacije pocetna tacka je jednaka vektoru zglobnih varijabli iz prethodnog
       koraka. Parametar "Init" je deklarisan u mapi "Data Properties/Input Ports" 
       i on je ulaz bloka S-funkcije. */
    
   ColumnVector startPoint(robotJointsNum);
   for (int i = 0; i < robotJointsNum; i++){
        startPoint.data[i][0] = Init[i];
    };
   
    
    /* Racunanje vektora zglobnih varijabli robotCurrentConfiguration, za odabran
       u funkciju kriterija datu sa (2.13), za robotski manipulator cija je kinematika
       odredjena matricom DH parametara A, za zeljenu pozu odredjenu sa matricom 
       homogene transformacije desiredPoseHTM i pocetnu aproksimaciju startPoint
       koja je jednaka vektoru zglobnih varijabli iz prethodne iteracije.
       Parametar "Config" je deklarisan u mapi "Data Properties/Output Ports" i on
       je izlaz bloka S-funkcije i predstavlja rezultat racunanja implementiranog 
       algoritma inverzne kinematike. */
    
    BFGSalgorithm algorithmInstance;
    ColumnVector robotCurrentConfiguration = algorithmInstance.BFGS(&CostFunction::criterion,
                                                                    DH, desiredPoseHTM, startPoint); 

    for (int i = 0; i < robotJointsNum; i++){
        Config[i] = robotCurrentConfiguration.data[i][0];
    };
    
    
    /* Vrijednost funkcije kriterija za robot cija kinematika odredjena matricom
       DH parametara A, za zeljenu poziciju i orijentaciju odredjenu matricom 
       homogene transformacije desiredPoseHTM, i za BFGS algoritmom izracunati
       vektor zglobnih varijabli  joint_Variables vector. Parametar "CriterionValue"
       je deklarisan u mapi "Data Properties/Output Ports" i on je jedan od izlaza
       bloka S-funkcije. */
    
    CostFunction costFunctionInstance;
    CriterionValue[0] =  costFunctionInstance.criterion(DH, desiredPoseHTM, robotCurrentConfiguration);
}


