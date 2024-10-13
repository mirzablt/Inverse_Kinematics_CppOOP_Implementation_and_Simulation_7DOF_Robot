
#ifndef INVERSEKINEMATICSLIBRARY_H
#define INVERSEKINEMATICSLIBRARY_H

#include "MatrixLibrary.h"



class HTMatrix : public Matrix {

    public:
        // Konstruktori.
        HTMatrix();
        HTMatrix(int m, int n);
        HTMatrix(const Matrix& matrix);
    
    private:
        /* Matrica homogene transformacije. Cetvrta kolona matrice DH parametara
           je kolona zglobnih varijabli robota i ove zglobne varijable su A.data[i][3]. */
        HTMatrix homogeneTransformMatrix(const Matrix& A, int i) const;
        
    public:
        /* Matrica (4x4) koja odredjuje jednacinu direktne kinematike koja se dobije 
           mnozenjem matrica A prema jednacini (1.10) */    
        HTMatrix directCinematicsMatrixMethod(const HTMatrix& A) const;
};



class Quaternion {
    public:
        double s;
        std::vector<double> v;

    public:
        // Konstruktori.
        Quaternion();
        Quaternion(double s, const std::vector<double>& v);
    
};



class HTMatrixQuaternion : public HTMatrix, public Quaternion {

    public:
        // Konstruktori
        HTMatrixQuaternion();
        HTMatrixQuaternion(const HTMatrix& htMatrix, const Quaternion& quaternion);

    public:
        // Transformacija koja prevodi matricu rotacije u kvaternion.
        Quaternion rotationMatrixToQuaternion(const HTMatrix& A) const;

        // Racunanje greske orijentacije izrazene u formi kvaterniona.
        Quaternion orientationError(const HTMatrix& dp_htm, const HTMatrix& ap_htm) const;    
};



class CostFunction : public HTMatrixQuaternion, public ColumnVector {

    public:
        // Konstruktori.
        CostFunction();
        CostFunction(const HTMatrixQuaternion& htMatrixQuaternion, const ColumnVector& columnVector);

    public:
        // Fukcija cilja (kriterij optimalnosti).
        static double criterion (HTMatrix& dh, const HTMatrix& dp_htm, const ColumnVector& var);
};



class GradientOfFunction : public  CostFunction {

    public:
        // konstruktori.
        GradientOfFunction();
        GradientOfFunction(const CostFunction& costFunction);

    public:
        // Gradijent funkcije f po promjenjivoj 'var', i on predstavlja vektor-kolonu.
        ColumnVector grad (double(*f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                HTMatrix& dh, const HTMatrix& dpc, const ColumnVector& var) const;

    public:
        // Linija odredjena sa (2.15) po kojoj vrsimo jednodimenzionalno pretrazivanje po parametru alpha.
        ColumnVector line (const ColumnVector& point, const ColumnVector& direction, double alpha) const;

        // Izvod funkcije f na pravoj odredjenom tackom 'point' i pravcem 'direction', za parametar t.
        double derivativeOnLine (double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                 HTMatrix& dh, const HTMatrix& dpc, const ColumnVector& point,
                                 const ColumnVector& direction, double t) const;

        // Drugi izvod funkcije f na pravoj odredjenom tackom 'point' i pravcem 'direction', za parametar t.
        double secondDerivative (double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                 HTMatrix& dh, const HTMatrix& dpc, const ColumnVector& point,
                                 const ColumnVector& direction, double t)const; 
};



class OnedimensionalLineSearch : public GradientOfFunction {
    
    public:
        // Konstruktori.
        OnedimensionalLineSearch();
        OnedimensionalLineSearch(const GradientOfFunction& gradientOfFunction);
        
    public:
        // Jednodimenzionalno pretrazivanje odredjeno sa (2.24, za odredjivanje optimalnog koraka prema (2.16).
        ColumnVector newtonTwoPoints (double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                         HTMatrix& dh, const HTMatrix& dpc, const ColumnVector& point,
                                         const ColumnVector& direction) const;
        
    public:
        /* Metod zoom implementira jedan od koraka elgoritma jednodimenzionalnog pretrazivanja
          u metodi strongWolfeLineSearch. */
        void zoom(double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                  HTMatrix& dh, const HTMatrix& dp_htm, const ColumnVector& x0,
                  const ColumnVector& direction, double alphaLow, double alphaHigh,
                  double fx0, double gradx0, double& alphaSol, double& fSol, double& gradSol);

        /* Metod strongWolfeLineSearch implementira jednodimenzionalno pretrazivanje
           zasnovanog na Wolfe uslovima. */ 
        void strongWolfeLineSearch(double( *f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                                   HTMatrix& dh, const HTMatrix& dp_htm, const ColumnVector& x0,
                                   const ColumnVector& direction, double& alphaSol);
    };



class BFGSalgorithm : public OnedimensionalLineSearch {

    public:
        // Konstruktori.
        BFGSalgorithm();
        BFGSalgorithm(const OnedimensionalLineSearch& onedimensionalLineSearch);

    public:
        // BFGS algoritam nelinearne optimizacije.
        ColumnVector BFGS ( double(*f)(HTMatrix&, const HTMatrix&, const ColumnVector&),
                           HTMatrix& dh, const HTMatrix& dp_htm, const ColumnVector& start_point);
};

#endif