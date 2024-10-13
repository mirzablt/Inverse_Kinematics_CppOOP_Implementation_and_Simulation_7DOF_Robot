
#ifndef MATRIXLIBRARY_H
#define MATRIXLIBRARY_H

#include <iostream>
#include <vector>
#include <math.h>

#define h 10e-5
#define pi 3.141592653



class Matrix {
    
    public:
        std::vector<std::vector<double>> data;

    public:
        // Konstruktori.
        Matrix();
        Matrix(int m, int n);

    // Metoda koja vraca dimenzije matrice.
    std::pair<int, int> size() const;

    
    // Operator overloading polimorfizam za operatore: +, -, *, ~
    
    //A + B
    Matrix operator+(const Matrix& B) const;

    // A - B
    Matrix operator-(const Matrix& B) const;

    // A * B
    Matrix operator*(const Matrix& B) const;

    // A * b
    Matrix operator*(double scalar) const;

    // ~A
    Matrix operator~() const;

    // a * B
    friend Matrix operator*(double scalar, const Matrix& A);
    

    // Definicije friend funkcija klase Matrix.
    
    //Jedinicna matrica dimenzija n x n.
    friend Matrix eye(int n);
    
    //Submatrica matrice A.
    friend Matrix getSubmatrix(const Matrix& A, int p, int q);
    
    // Determinanta matrice A.
    friend double det(const Matrix& A, int n);

    // Adjungovana matrica matrice A.
    friend Matrix adjung(Matrix& A);
    
    // Inverzna matrica matrice A.
    friend Matrix inv(Matrix& A);
};



class ColumnVector : public Matrix {
    
    public:
        //Konstruktori.
        ColumnVector();
        ColumnVector(int m);
        ColumnVector(const Matrix& A);

    public:
        // Overloading operator [].
        double& operator[](int index);

        // Euklidska norma vektora.
        double norm(const ColumnVector& v) const;
};

#endif

