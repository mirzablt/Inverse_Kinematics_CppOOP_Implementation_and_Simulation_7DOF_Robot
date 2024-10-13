
#include "InverseKinematicsLibrary.h"



//Klasa HTMatrix. Naslijedjuje klasu Matrix.


    // Konstruktori klase HTMatrix.

    HTMatrix::HTMatrix() : Matrix(4, 4) {}
    
    HTMatrix::HTMatrix(int m, int n) : Matrix(m, n) {}
    
    
    // Konstruktor koji omogucava kreiranje HTMatrix objekta iz postojeceg Matrix objekta.
    
    HTMatrix::HTMatrix(const Matrix& matrix) {
            data = matrix.data;
    }

    
    /* Metoda homogeneTransformMatrix racuna  matricu homogene transformacije B, 
       odredjena relacijom (1.14), na osnovu matrice DH parametara A. A.data[i][0] 
       su parametri a; A.data[i][1] su parametri alfa; A.data[i][2] su d parametri;
       A.data[i][3] su theta parammetri (zglobne varijable). */
    
    HTMatrix HTMatrix::homogeneTransformMatrix(const Matrix& A, int i) const {

        HTMatrix B(4, 4);
        
        B.data[0][0] =  cos(A.data[i][3]);                          B.data[0][1] = -sin(A.data[i][3]) * cos(A.data[i][1]);    
        B.data[1][0] =  sin(A.data[i][3]);                          B.data[1][1] =  cos(A.data[i][3]) * cos(A.data[i][1]);    
        B.data[2][0] =  0.0;                                        B.data[2][1] =  sin(A.data[i][1]);                                  
        B.data[3][0] =  0.0;                                        B.data[3][1] =  0.0;                              
    
        B.data[0][2] =  sin(A.data[i][3]) * sin(A.data[i][1]);      B.data[0][3] =  A.data[i][0] * cos(A.data[i][3]);
        B.data[1][2] = -cos(A.data[i][3]) * sin(A.data[i][1]);      B.data[1][3] =  A.data[i][0] * sin(A.data[i][3]);
        B.data[2][2] =  cos(A.data[i][1]);                          B.data[2][3] =  A.data[i][2];
        B.data[3][2] =  0.0;                                        B.data[3][3] =  1.0;
        
        return B; 
    }

    
    /* Jednacina direktne kinematike definise ovisnost pozicije i orijentacije 
       vrha manipulatora od zglobnih varijabli. Metoda directCinematicsMatrixMethod
       daje matricu (4x4) koja odredjuje jednacinu direktne kinematike robota, 
       koja se dobije mnozenjem matrica homogene transf. A prema jednacini (1.10). */
    
    HTMatrix HTMatrix::directCinematicsMatrixMethod(const HTMatrix& A) const {

        HTMatrix M(A.data.size(), A.data[0].size());
        HTMatrix N = homogeneTransformMatrix(A, 0);

        for (int i = 1; ; ){
            M = N * homogeneTransformMatrix(A, i);
            i++;
            if(i >= A.data.size()) break;
                N = M * homogeneTransformMatrix(A, i);
                i++;
            if(i >= A.data.size()) break;	
        }

        if(A.data.size()%2 == 0) {return M;}
        else {return N;}
    } 

    

//Klasa Quaternion.


    //Konstruktori klase Quaternion.
    
    Quaternion::Quaternion() : s(0.0), v(3, 0.0) {}

    Quaternion::Quaternion(double s, const std::vector<double>& v) : s(s), v(v) {}

    

//Klasa HTMatrixQuaternion.Naslijedjuje klase: HTMatrix, Quaternion
  
    
    // Konstruktori.
    
    HTMatrixQuaternion::HTMatrixQuaternion() : HTMatrix(4, 4), Quaternion() {}

    HTMatrixQuaternion::HTMatrixQuaternion(const HTMatrix& htMatrix, const Quaternion& quaternion)
        : HTMatrix(htMatrix), Quaternion(quaternion) {}

        
    /* Orijentacija vrha manipulatora je odredjena matricom rotacije. Ova orijentacija
       se moze izraziti u terminima kvaterniona.
       Metoda rotationMatrixToQuaternion transformise matricu rotacije (koja je jednaka
       gornjoj trougaonoj matrici matrice homogene transf. A) u kvaternion 
       Q = [s, v0, v1, v2] prema relaciji (2.8) uz tretiranje granicnih slucajeva. */
        
    Quaternion HTMatrixQuaternion::rotationMatrixToQuaternion(const HTMatrix& A) const {
        Quaternion Q;

        double trace = A.data[0][0] + A.data[1][1] + A.data[2][2];

        if (trace > h) {

            double s = 0.5 / std::sqrt(trace + 1.0);
            Q.s = 0.25 / s;
            Q.v[0] = (A.data[2][1] - A.data[1][2]) * s;
            Q.v[1] = (A.data[0][2] - A.data[2][0]) * s;
            Q.v[2] = (A.data[1][0] - A.data[0][1]) * s;

        } else {

            if (A.data[0][0] - A.data[1][1] > h && A.data[0][0] - A.data[2][2] > h) {
                double s = 2.0 * std::sqrt(1.0 + A.data[0][0] - A.data[1][1] - A.data[2][2]);
                Q.s = (A.data[2][1] - A.data[1][2]) / s;
                Q.v[0] = 0.25 * s;
                Q.v[1] = (A.data[0][1] + A.data[1][0]) / s;
                Q.v[2] = (A.data[0][2] + A.data[2][0]) / s;

            } else if (A.data[1][1] - A.data[2][2] > h) {

                double s = 2.0 * std::sqrt(1.0 + A.data[1][1] - A.data[0][0] - A.data[2][2]);
                Q.s = (A.data[0][2] - A.data[2][0]) / s;
                Q.v[0] = (A.data[0][1] + A.data[1][0]) / s;
                Q.v[1] = 0.25 * s;
                Q.v[2] = (A.data[1][2] + A.data[2][1]) / s;

            } else {

                double s = 2.0 * std::sqrt(1.0 + A.data[2][2] - A.data[0][0] - A.data[1][1]);
                Q.s = (A.data[1][0] - A.data[0][1]) / s;
                Q.v[0] = (A.data[0][2] + A.data[2][0]) / s;
                Q.v[1] = (A.data[1][2] + A.data[2][1]) / s;
                Q.v[2] = 0.25 * s;
            }
        }
        return Q;
    }

    
    /* Greska orijentacije izrazava odstupanje trenutne orijentacije vrha manipulatora
       od zeljene orijentacije. Ovo odstupanje se moze izraziti u obliku matrice (2.7).
       Ova matrica se dobije na osnovu matrica dp_htm i ap_htm i ona se moze prevesti
       u kvaternion, ciji vektorski dio (relacija 2.12) opisuje  ovu gresku orijentacije
       u obliku (2.12). Naredna metoda daje vektorski dio kvaterniona na osnovu matrica
       dp_htm i ap_htm. */
    
     Quaternion HTMatrixQuaternion::orientationError(const HTMatrix& dp_htm, const HTMatrix& ap_htm) const {

            Quaternion D = rotationMatrixToQuaternion(dp_htm);  // dp_htm - desired pose homogene transform matrix
            Quaternion E = rotationMatrixToQuaternion(ap_htm);  // ap_htm - actual pose homogene transform matrix
            Quaternion delta;

            delta.s = 1.0;
            delta.v[0] = E.s * D.v[0] - D.s * E.v[0] - (-D.v[2] * E.v[1] + D.v[1] * E.v[2]);
            delta.v[1] = E.s * D.v[1] - D.s * E.v[1] - (D.v[2] * E.v[0] - D.v[0] * E.v[2]);
            delta.v[2] = E.s * D.v[2] - D.s * E.v[2] - (-D.v[1] * E.v[0] + D.v[0] * E.v[1]);

            return delta;
        }    



// Klasa CostFunction. Naslijedjuje klase: HTMatrixQuaternion i ColumnVector.

     
     //Konstruktori.
     
    CostFunction::CostFunction() : HTMatrixQuaternion(), ColumnVector(0) {}
 
    CostFunction::CostFunction(const HTMatrixQuaternion& htMatrixQuaternion, const ColumnVector& columnVector)
        : HTMatrixQuaternion(htMatrixQuaternion), ColumnVector(columnVector) {}

        
    /* U metodi criterion je implementirana fukcija cilja (2.13) koja izrazava mjeru
       odstupanja trenutne pozicije i orijentacije vrha manipulatora od njegove zeljene
       pozicije i orijentacije. Funkcija cilja u sebe ukljucuje:
          -positionPart - ovaj dio izrazava rastojanje trenutne pozicije vrha manipulatora
          (koji je funkcija zglobnih varijabli) od zeljene pozicije;
          -orientationPart -  ovaj dio izrazava ostupanje trenutne orijentacije vrha
           manipulatora (koji je funkcija zglobnih varijabli) od trenutne orijentacije.
           Ovo odstupanje se moze izraziti vektorskim dijelom kvaterniona.
       Matrica tezinskih koeficijenata pozicije Mp iz (2.13) je dijagonalna i njeni elementi
       na glavnoj dijagonali su uzeti da su jednaki i iznose: positionWeights. Matrica 
       tezinskih koeficijenata orijentacije Mo iz (2.13) je dijagonalna i njeni elementi
       na glavnoj dijagonalisu uzeti da su jednaki  i iznose: orientationWeights. Zeljena
       pozicija i orijentacija vrha manipulatora se dobiju na osnovu matrice homogene 
       transf. dp_htm.
       Trenutna poz. i orij. se dobiju iz jednacine direktne kinematike, koja se dobije na 
       osnovu dh parametara robota koje cine cine izglobne varijable qi. U funkciji cilja
       'criterion', var.data[i][0] (i=1,2,...,n) su varijable one predstavljaju zglobne
       promjenjive q1, q2, ...,qn.   */
        
    double CostFunction::criterion (HTMatrix& dh, const HTMatrix& dp_htm, const ColumnVector& var){
     
        for (int i = 0; i < dh.data.size(); i++){
            dh.data[i][3] = var.data[i][0];
        }

        HTMatrix HTMatrixInstance;
        HTMatrixQuaternion HTMatrixQuaternionInstance;

        HTMatrix dcm = HTMatrixInstance.directCinematicsMatrixMethod(dh);
        Quaternion delta = HTMatrixQuaternionInstance.orientationError(dp_htm, dcm);

        double positionPart = 0.0;
        double orientationPart = 0.0;

        const double positionWeights = 1;
        const double orientationWeights = 0.01;

       for(int i = 0; i < 3; i++){
          positionPart += positionWeights * pow(dp_htm.data[i][3] - dcm.data[i][3], 2);
          orientationPart += orientationWeights * pow(delta.v[i], 2);
       };

        return positionPart + orientationPart;
    }


    
// Klasa GradientOfFunction. Naslijedjuje klasu CostFunction.u kojoj su definisane metode:

    
    //Konstruktori.
    
    GradientOfFunction::GradientOfFunction() : CostFunction() {}

    GradientOfFunction::GradientOfFunction(const CostFunction& costFunction)
        : CostFunction(costFunction) {}


    /* Metod grad racuna gradijent funkcije f po promjenjivoj 'var'. Kako racunamo gradijent
       funkcije cilja, to ovom metodu moramo proslijediti i parametre analogno kao u slucaju
       metoda u kojem je implementirana funk. cilja.  */
        
    ColumnVector GradientOfFunction::grad (double(*f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                            HTMatrix& dh, const HTMatrix& dpc, const ColumnVector& var) const {
        
        int n = var.data.size();
        ColumnVector var_1(n);
        ColumnVector var_2(n);
        ColumnVector partialDerivatives(n);
        
        for(int i=0; i<var.data.size(); i++){
            var_1 = var;
            var_2 = var;
            var_1.data[i][0] += h;
            var_2.data[i][0] -= h;
            partialDerivatives.data[i][0] = ( f(dh, dpc, var_1) - f(dh, dpc, var_2) ) / (2.0*h);
        }
        return partialDerivatives;	
    }

    
    /* Metoda line implementira lednacinu (2.15), a to je linija po kojoj vrsimo 
       jednodimenzionalno pretrazivanje po parametru alpha. Ona je odredjena tackom 
       point i vektorom pravca direction. */
    
    ColumnVector GradientOfFunction::line(const ColumnVector& point, const ColumnVector& direction, double alpha) const {

        ColumnVector linea(point.data.size());

        for (int i = 0; i < point.data.size(); i++){
            linea.data[i][0] = point.data[i][0] + alpha * direction.data[i][0];
        }
        return linea;
    }

    
    /* derivativeOnLine daje izvod funkcije f na pravoj odredjenom tackom 'point' 
       i pravcem 'direction', za parametar t. */
    
    double GradientOfFunction::derivativeOnLine (double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                         HTMatrix& dh, const HTMatrix& dpc, const ColumnVector& point,
                         const ColumnVector& direction, double t) const {
       
        double df = f(dh, dpc, line(point, direction, t + h)) -
                    f(dh, dpc, line(point, direction, t - h));
    
        return df / (2.0 * h);
    }

    
    /* secondDerivative daje drugi izvod funkcije f na pravoj odredjenom tackom 'point' 
       i pravcem 'direction', za parametar t. */
    
    double GradientOfFunction::secondDerivative (double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                         HTMatrix& dh, const HTMatrix& dpc, const ColumnVector& point,
                         const ColumnVector& direction, double t)const{
    
        double ddf = derivativeOnLine (f, dh, dpc, point, direction, t + h) -
                     derivativeOnLine (f, dh, dpc, point, direction, t - h);

        return ddf / (2.0 * h); 
    } 



//Klasa OnedimensionalLineSearch. Naslijedjuje klasu OnedimensionalLineSearch.
  
    
    //Konstruktori.
    
    OnedimensionalLineSearch::OnedimensionalLineSearch() : GradientOfFunction() {}

    OnedimensionalLineSearch::OnedimensionalLineSearch(const GradientOfFunction& gradientOfFunction)
        : GradientOfFunction(gradientOfFunction) {}


    /* Metod newtonTwoPoints implementira jednodimenzionalno pretrazivanje odredjeno sa (2.24),
       za odredjivanje optimalnog koraka duz linije odredjene pravcem 'direction' i tackom 'point'.
       Kriterij zaustavljanja (eps) je dostizanje minimalne  promjene problemske varijable. 
       Pocetna vrijednost x1k je uzeta proizvoljno. */
        
    ColumnVector OnedimensionalLineSearch::newtonTwoPoints (double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                     HTMatrix& dh, const HTMatrix& dpc, const ColumnVector& point,
                                     const ColumnVector& direction) const {

        const double EPSILON =0.0001;
        int i;
        double xk1;
        double eps;
        
        ColumnVector XK = point;
        ColumnVector optimalStep(1);
        
        double xk = 0.0; 
        double x1k = 0.2; 

        do{
            double dxk = xk - x1k;
            double dDerivativeOnLineK = derivativeOnLine( f, dh, dpc, XK, direction, xk) -
                                        derivativeOnLine( f, dh, dpc, XK, direction, x1k);

            xk1 = x1k + (dxk / (1 - ((derivativeOnLine(f, dh, dpc, XK, direction, xk)/
                         derivativeOnLine(f, dh, dpc, XK, direction, x1k))*
                         (dDerivativeOnLineK /(dxk * secondDerivative( f, dh, dpc, XK, direction, xk))))));

            eps = fabs(xk1 - xk);
            x1k = xk;
            xk = xk1;

        }while (eps > EPSILON);

        optimalStep.data[0][0] = xk1;

        return optimalStep;
    }

    
    /* Metod zoom implementira algoritam 3.6 (str. 61: J. Nocedal, S. J. Wright. Numerical
       Optimization. Springer, 2006.). Ovaj metod se poziva u  metodu strongWolfeLineSearch
       koji implemetira jednodimenzionalno pretrazivanje zasnovano na Wolfe uslovima (2.25). 
       Ova mmetoda i met. strongWolfeLineSearc je implementirana koristeci impl. u Matlab-u:
       https://github.com/hiroyuki-kasai/GDLibrary/blob/master/line_search/strong_wolfe_line_search.m */
    
    void OnedimensionalLineSearch::zoom(double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                        HTMatrix& dh, const HTMatrix& dp_htm, const ColumnVector& x0,
                                        const ColumnVector& direction, double alphaLow, double alphaHigh,
                                        double fx0, double gradx0, double& alphaSol, double& fSol, double& gradSol) {

        const double c1 = 1e-4;
        const double c2 = 0.5;
        int i = 0;
        int maxIter = 5;

        while (true) {

            double alphaCur = 0.5 * (alphaLow + alphaHigh);
            alphaSol = alphaCur;
            ColumnVector xCur = line(x0, direction, alphaCur);

            double fxCur = f(dh, dp_htm, xCur);
            ColumnVector gradxCur = grad(f, dh, dp_htm, xCur);

            fSol = fxCur;
            gradSol = ((~gradxCur) * direction).data[0][0];

            ColumnVector xLow = line(x0, direction, alphaLow);

            double fxLow = f(dh, dp_htm, xLow);

            if ((fxCur > fx0 + c1 * alphaCur * gradx0) || (fxCur >= fxLow)) {
                alphaHigh = alphaCur;
            } else {
                if (fabs(gradSol) <= -c2 * gradx0) {
                    alphaSol = alphaCur;
                    return;
                }
                if (gradSol * (alphaHigh - alphaLow) >= 0) {
                    alphaHigh = alphaLow;
                }
                alphaLow = alphaCur;
            }
            i++;

            if (i > maxIter) {
                alphaSol = alphaCur;
                return;
            }
        }
    }
    
    
    /* Metoda strongWolfeLineSearch implementira jednodimenzionalno pretrazivanje po pravoj 
       odredjenom pravcem "direction" i tackom x0. Ona racuna vrijednost optimalnog koraka
       alphaSol za funkciju cilja f. Ovo jednodimenzionalno pretrazivanje je zasnovano na Wolfe
       uslovima( algoritam 3.5, str. 60: J. Nocedal, S. J. Wright. Numerical Optimization. 
       Springer, 2006.) */

    void OnedimensionalLineSearch::strongWolfeLineSearch(double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                                         HTMatrix& dh, const HTMatrix& dp_htm, const ColumnVector& x0,
                                                         const ColumnVector& direction, double& alphaSol) {

        double fx0 = f(dh, dp_htm, x0);
        ColumnVector gradx0 = grad(f, dh, dp_htm, x0);

        int maxIter = 3;
        double alphaMax = 20;
        double alphaPre = 0;
        double alphaCur = 1;

        double gradInx0OnDirection = ((~gradx0) * direction).data[0][0];
        double fxPrev = fx0;
        double gradxPrev = gradInx0OnDirection;
        int i = 1;

        const double c1 = 1e-4;
        const double c2 = 0.5;

        double fSol;
        double gradSol;

        while (true) {
            ColumnVector xCur = line(x0, direction, alphaCur);

            double fxCur = f(dh, dp_htm, xCur);
            ColumnVector gradxCur = grad(f, dh, dp_htm, xCur);

            fSol = fxCur;
            gradSol = ((~gradxCur) * direction).data[0][0]; 

            if ((fxCur > fx0 + c1 * alphaCur * gradInx0OnDirection) || ((i > 1) && (fxCur >= fxPrev))) {
                zoom(f, dh, dp_htm, x0, direction, alphaPre, alphaCur, fx0, gradInx0OnDirection, alphaSol, fSol, gradSol);
                return;
            }
            if (fabs(gradSol) <= -c2 * gradInx0OnDirection) {
                alphaSol = alphaCur;
                return;
            }
            if (gradSol >= 0) {
                zoom(f, dh, dp_htm, x0, direction, alphaCur, alphaPre, fx0, gradInx0OnDirection, alphaSol, fSol, gradSol);
                return;
            }
            alphaPre = alphaCur;
            fxPrev = fxCur;
            gradxPrev = gradSol;

            if (i > maxIter) {
                alphaSol = alphaCur;
                return;
            }
            const double r = 0.8;
            alphaCur = alphaCur + (alphaMax - alphaCur) * r;
            i++;
        }
    }



// Klasa BFGSalgorithm. Naslijedjuje klasu OnedimensionalLineSearch.

               
    // Konstruktori

    BFGSalgorithm::BFGSalgorithm() : OnedimensionalLineSearch() {}

    BFGSalgorithm::BFGSalgorithm(const OnedimensionalLineSearch& onedimensionalLineSearch)
    : OnedimensionalLineSearch(onedimensionalLineSearch) {}
    

    /* Implementacija algoritma opisanog u 2.5. Algoritam odredjuje vektor zglobnih 
       varijabli za kojeg funkcija cilja f ima minimalnu vrijednost. Ta vrijednost 
       je bliska nuli (EPSILON = 0.0001). U idealnom slucaju ta vrijednost je nula,
       sto odgovara slucaju da se aktuelna pozicija (i orijentacija) i zeljena poz.
       (i orij.) vrha manipulatora poklapaju. 
       Zeljena pozicija i orijentacija vrha manipulatora se dobiju na osnovu matrice 
       homogene transf. dp_htm, koju generira algoritam planiranja trajektorije. Ova
       matrica je ulaz bloka S-funkcije. Trenutna poz. i orij. se dobiju iz jednacine
       direktne kinematike, koja se dobije na osnovu dh parametara robota koje cine 
       zglobne varijable qi.
       Xk.data[i] (i=1,2,...,n) je vektor problemskih varijabli (i predstavlja vektor 
       zglobnih varijabli q1, q2, ...,qn). */
       
    ColumnVector BFGSalgorithm::BFGS ( double(*f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                       HTMatrix& dh, const HTMatrix& dp_htm, const ColumnVector& startPoint){

        const double EPSILON = 0.0001;                    

        int n = startPoint.data.size();

        ColumnVector Xk(n), Xk1(n), DXk(n);
        ColumnVector Gk(n), Gk1(n), DGk(n);
        ColumnVector step(1), Rk(n);
        
        Matrix Hk, Hk1;
        Matrix Mk, Nk;

        double gamma, beta;	
        Xk = startPoint;                                               // Pocetna aproksimacija
        Hk = eye(startPoint.data.size());                              // Pocetna vrijednost matrice H

        double alpha;                                                  // Optimalni korak jednodimenzionalnog pretrazivanja.

        while (norm(grad(f, dh, dp_htm, Xk)) > EPSILON){               // Kriterij zaustavljanja dat sa (2.59).
            
            Gk    =  grad(f, dh, dp_htm, Xk);
            Rk    =  -1 * (Hk * Gk);                                   // Pravac pretrazivanja dat sa (2.60).
            
            //alpha  =  newtonTwoPoints (f, dh, dp_htm, Xk, Rk);       // Jednodimenzionalno pretrazivanje (2.61).
            strongWolfeLineSearch(f, dh, dp_htm, Xk, Rk, alpha);       // Jednodimenzionalno pretrazivanje odredjeno Wolfeovim uslovima

            Xk1   =  Xk + (Rk * alpha);                                // Racunanje nove aproksimacije (2.62).
            Gk1   =  grad(f, dh, dp_htm, Xk1);

            DXk   =  Xk1 - Xk;                                         // Razlika data sa(2.63a).
            DGk   =  Gk1 - Gk;                                         // Razlika data sa(2.63b).

            gamma =  ((~DXk) * DGk).data[0][0];
            beta  =  (((~DGk) * Hk) * DGk).data[0][0];

            Mk    =  (DXk * (~DXk)) * (gamma + beta) * (1/ pow(gamma,2));
            Nk    =  ( ( (Hk * DGk) * (~DXk)) + ( (DXk * (~DGk)) * Hk)) *  (-1.0/gamma);

            Hk1   = (Hk + Mk) + Nk;                                    // Racunanje aproksimacije inverznog Hessiana H(k+1) prema (2.56).

            Hk    =  Hk1;
            Xk    =  Xk1;	
        }
        return Xk1;
    }


