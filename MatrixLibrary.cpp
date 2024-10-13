
#include "MatrixLibrary.h"


// Klasa Matrix.


    // Konstruktori.

    Matrix::Matrix() {}

    Matrix::Matrix(int m, int n) : data(m, std::vector<double>(n, 0.0)) {}

    
    // Metoda koja vraca dimenzije  matrice.
    
    std::pair<int, int> Matrix::size() const {
        return {static_cast<int>(data.size()), static_cast<int>(data[0].size())};
    }

    
    // Overloading operatora +, koji omogucava sabiranje matrica A i B u obluku A + B.
    
    Matrix Matrix::operator+(const Matrix& B) const {
        int m = data.size();
        int n = data[0].size();

        Matrix C(m, n);

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                C.data[i][j] = data[i][j] + B.data[i][j];
            }
        }
        return C;
    }

    
    // Overloading operatora -, koji omogucava racunanje razlike matrica A i B u obluku A - B. 
    
    Matrix Matrix::operator-(const Matrix& B) const {
        int m = data.size();
        int n = data[0].size();

        Matrix C(m, n);

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                C.data[i][j] = data[i][j] - B.data[i][j];
            }
        }
        return C;
    }
    
    
   // Overloading operatora *, koji omogucava mnozenje matrica A i B u obluku A * B.
    
   Matrix Matrix::operator*(const Matrix& B) const {
        int m = size().first;
        int n = size().second;
        int p = B.size().second;

        Matrix C(m, p);

        for (int i = 0; i < m; ++i) {
            for (int k = 0; k < p; ++k) {
                double s = 0.0;
                for (int j = 0; j < n; ++j) {
                    s += data[i][j] * B.data[j][k];
                }
                C.data[i][k] = s;
            }
        }
        return C;
    }


    // Overloading operatora *, koji omogucava mnozenje skalara a i matrice B u obluku  B * a.
   
    Matrix Matrix::operator*(double scalar) const {
        int m = size().first;
        int n = size().second;

        Matrix B(m, n);

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                B.data[i][j] = scalar * data[i][j];
            }
        }
        return B;
    }


    // Overloading operatora *, koji omogucava mnozenje skalara a i matrice B u obluku  a * B.
    Matrix operator*(double scalar, const Matrix& A) {
        return A * scalar;
    }

    
    // Overloading operatora  ~, koji omogucava transponovanje matrice A u obliku ~A.
    
    Matrix Matrix::operator~() const {
        int rows = data.size();
        int cols = data[0].size();
        
        Matrix B(cols, rows);
        
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                B.data[j][i] = data[i][j];
            }
        }
        return B;
    }



    /* Definicija friend funkcije klase Matrix koja daje matricu koja se dobije
       izostavljanjem p-te vrste i q-te kolone matrice A. */
    
    Matrix getSubmatrix(const Matrix& A, int p, int q) {
        int n = A.data.size();
        Matrix B(n - 1, n - 1);
        int i = 0, j = 0; 

        for (int row = 0; row < n; ++row) {
            for (int col = 0; col < n; ++col) {
                if (row != p && col != q) {
                    B.data[i][j++] = A.data[row][col];
                    if (j == n - 1) {
                        j = 0;
                        ++i;
                    }
                }
            }
        }
        return B;
    }

    
    // Friend funkcija klase Matrix koja omogucava racunanje jedinicne matrice reda n.
    
    Matrix eye(int n) {
        Matrix A(n, n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A.data[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return A;
    }

    
    // Friend funkcija klase Matrix koja omogucava racunanje determinante matrice A.
    
    double det(const Matrix& A, int n){
        double d = 0;
        if(n == 1)
           return A.data[0][0];  
        int sign = 1;
        for(int i = 0; i < n; i++) {
            d += sign * A.data[0][i] * det(getSubmatrix(A,0,i), n-1);
            sign = -sign;
        }
        return d;
    }

    
    // Friend funkcija klase Matrix koja omogucava racunanje adjungovane matrice.
    
    Matrix adjung (const Matrix& A){

        int n = A.data.size();
        Matrix B(n, n);

        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                B.data[j][i] = pow(-1,i+j) * det(getSubmatrix(A, i, j), n-1);
            }
        }
        return B;		
    }

    
    // Friend funkcija klase Matrix koja omogucava racunanje inverzne matrice .
    
    Matrix inv(const Matrix& A){

        int n = A.data.size();
        Matrix C(n, n);
        Matrix B = adjung(A);

        if(det(A, n) != 0){
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    C.data[i][j] = B.data[i][j] / det(A, n);
                }
            }
        return C;
       }
    }



//Klasa ColumnVector. Naslijedjuje klasu Matrix.
    
    
    // Konstruktori.

    ColumnVector::ColumnVector() : Matrix(7, 1) {}

    ColumnVector::ColumnVector(int m) : Matrix(m, 1) {}
    

    // Konstruktor koji omogucava kreiranje ColumnVector objekta iz postojeceg Matrix objekta.
    
    ColumnVector::ColumnVector(const Matrix& A) {
        
        if (A.data.size() > 0 && A.data[0].size() > 0) {          
            data = A.data;
        }
    }

    
    // Overloading operator [] koji omogucava da vektore obloka v.data[i][0] pisemo u obliku v.data[i].
    
    double& ColumnVector::operator[](int index) {
        return data[index][0];
    }

    
    //Euklidska norma vektora v.
    
    double ColumnVector::norm(const ColumnVector& v) const {
        
        double norm = 0.0;
        for (int i = 0; i < v.data.size(); i++)
            norm += v.data[i][0] * v.data[i][0];

	return sqrt(norm);      
    }



