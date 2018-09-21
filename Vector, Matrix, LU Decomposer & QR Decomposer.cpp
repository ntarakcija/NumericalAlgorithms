#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <limits>
#include <iomanip>

class Vector {
    std::vector<double> v;
    public:
    explicit Vector(int n);
    Vector(std::initializer_list<double> l);
    int NElems() const { return v.size(); };
    double &operator[](int i) { return v[i]; };
    double operator[](int i) const { return v[i]; };
    double &operator()(int i);
    double operator()(int i) const;
    double Norm() const;
    friend double VectorNorm(const Vector &v);
    double GetEpsilon() const { return this->Norm() * std::numeric_limits<double>::epsilon() * 10; };
    void Print(char separator = '\n', double eps = -1) const;
    friend void PrintVector(const Vector &v, char separator, double eps);
    friend Vector operator +(const Vector &v1, const Vector &v2);
    Vector &operator +=(const Vector &v);
    friend Vector operator -(const Vector &v1, const Vector &v2);
    Vector &operator -=(const Vector &v);
    friend Vector operator *(double s, const Vector &v);
    friend Vector operator *(const Vector &v, double s);
    Vector &operator *=(double s);
    friend double operator *(const Vector &v1, const Vector &v2);
    friend Vector operator /(const Vector &v, double s);
    Vector &operator /=(double s);
    void Chop(double eps = -1);
    bool EqualTo(const Vector &v, double eps = -1) const;
};

class Matrix {
    std::vector<std::vector<double>> matrica;
    void RazmjenaKolona(int a, int b);
    void RazmjenaRedova(int a, int b);
    public:
    Matrix() {}
    Matrix(int m, int n);
    Matrix(const Vector &v);
    Matrix(std::initializer_list<std::vector<double>> l);
    int NRows() const { return matrica.size(); }
    int NCols() const { return matrica[0].size(); }
    double *operator[](int i) { return &matrica[i][0]; }
    const double *operator[](int i) const { return &matrica[i][0]; }
    double &operator()(int i, int j);
    double operator()(int i, int j) const;
    double Norm() const;
    friend double MatrixNorm(const Matrix &m);
    double GetEpsilon() const { return this->Norm() * std::numeric_limits<double>::epsilon() * 10; }
    void Print(int width = 10, double eps = -1) const;
    friend void PrintMatrix(const Matrix &m, int width, double eps);
    friend Matrix operator +(const Matrix &m1, const Matrix &m2);
    Matrix &operator +=(const Matrix &m);
    friend Matrix operator -(const Matrix &m1, const Matrix &m2);
    Matrix &operator -=(const Matrix &m);
    friend Matrix operator *(double s, const Matrix &m);
    friend Matrix operator *(const Matrix &m, double s);
    Matrix &operator *=(double s);
    friend Matrix operator *(const Matrix &m1, const Matrix &m2);
    Matrix &operator *=(const Matrix &m);
    friend Vector operator *(const Matrix &m, const Vector &v);
    friend Matrix Transpose(const Matrix &m);
    void Transpose();
    void Chop(double eps = -1);
    bool EqualTo(const Matrix &m, double eps = -1) const;
    friend Matrix LeftDiv(Matrix m1, Matrix m2);
    friend Vector LeftDiv(Matrix m, Vector v);
    friend Matrix operator /(const Matrix &m, double s);
    Matrix &operator /=(double s);
    friend Matrix operator /(Matrix m1, Matrix m2);
    Matrix &operator /=(Matrix m);
    double Det() const;
    friend double Det(Matrix m);
    void Invert();
    friend Matrix Inverse(Matrix m);
    void ReduceToRREF();
    friend Matrix RREF(Matrix m);
    int Rank() const;
    friend int Rank(Matrix m);
};

class LUDecomposer {
    Matrix l, u, lu;
    std::vector<int> w;
    public:
    LUDecomposer(Matrix m);
    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;
    Matrix GetCompactLU() const;
    Matrix GetL() const { return l; }
    Matrix GetU() const { return u; }
    Vector GetPermuation() const;
};

class QRDecomposer {
    public:
    //QRDecomposer(Matrix m);
    //void Solve(const Vector &b, Vector &x) const;
    //Vector Solve(Vector b) const;
    //void Solve(Matrix &b, Matrix &x) const;
    //Matrix Solve(Matrix b) const;
    //Vector MulQWith(Vector v) const;
    //Matrix MulQWith(Matrix m) const;
    //Vector MulQTWith(Vector v) const;
    //Matrix MulQTWith(Matrix m) const;
    //Matrix GetQ() const;
    //Matrix GetR() const;
};

// QRDecomposer ------------------------------------------------------------------------------------


// LUDecomposer ------------------------------------------------------------------------------------
LUDecomposer::LUDecomposer(Matrix m) {
    if(m.NRows() != m.NCols()) throw std::domain_error("Matrix is not square");
    if(std::abs(m.Det()) < m.GetEpsilon()) throw std::domain_error("Matrix is singular");
    
    int n(m.NRows()), p;
    double s, mi;
    w.resize(m.NRows());
    
    for(int j(0); j < n; j++) {
        for(int i(0); i < j; i++) {
            s = m[i][j];
            for(int k(0); k < i - 1; k++)
                s = s - m[i][k] * m[k][j];
            m[i][j] = s;
        }    
        p = j;
        for(int i(j + 1); i < n; i++) {
            s = m[i][j];
            for(int k(0); k < j - 1; k++)
                s = s - m[i][k] * m[k][j];
            m[i][j] = s;
            if(std::abs(s) > std::abs(m[p][j]))
                p = i;
        }
        if(std::abs(m[p][j]) < m.GetEpsilon())
            throw std::domain_error("Matrix is singular");
        if(p != j) {
            for(int a(0); a < m.NCols(); a++) {
                double temp(m[p][a]);
                m[p][a] = m[j][a];
                m[j][a] = temp;
            }
        }
        w[j] = p;
        mi = m[j][j];
        for(int i(j + 1); i < n; i++)
            m[i][j] = m[i][j] / mi;
    }
    
    Matrix l1(m.NRows(), m.NRows());
    Matrix u1(m.NRows(), m.NRows());
    
    for(int i(0); i < m.NRows(); i++)
        for(int j(0); j < m.NRows(); j++) {
            if(j >= i) {
                u1[i][j] = m[i][j];
                if(i == j) l1[i][j] = 1;
            }
            else l1[i][j] = m[i][j];
        }
    
    u = u1;
    l = l1;
}

void LUDecomposer::Solve(const Vector &b, Vector &x) const {
    if(b.NElems() != u.NRows() || x.NElems() != u.NRows()) throw std::domain_error("Incompatible formats");
    
    int n(l.NRows());
    double s;
    Vector y(l.NRows());
    Vector b1(b);
    
    for(int i(0); i < n; i++) {
        s = b[i];
        for(int j(0); j < i - 1; j++)
            s = s - l[i][j] * y[j];
        y[i] = s;
    }

    for(int i(n - 1); i >= 0; i--) {
        s = y[i];
        for(int j(i + 1); j < n; j++)
            s = s - u[i][j] * x[j];
        x[i] = s / u[i][i];
    }
}

Vector LUDecomposer::Solve(Vector b) const {
    Vector x(l.NRows());
    Solve(b, x);
    return x;
}

void LUDecomposer::Solve(Matrix &b, Matrix &x) const {
    if(b.NRows() != u.NRows() || b.NCols() != u.NCols() || x.NRows() != u.NRows() || x.NCols() != u.NCols())
        throw std::domain_error("Incompatible formats");
    
    int n(l.NRows()), m(b.NCols());
    double s;
    Matrix y(n, n);
    
    for(int k(0); k < m; k++) {
        for(int i(0); i < n; i++) {
            s = b[i][k];
            for(int j(0); j < i - 1; j++)
                s = s - l[i][j] * y[j][k];
            y[i][k]= s;
        }
    }
    for(int k(0); k < m; k++) {
        for(int i(n - 1); i >= 0; i--) {
            s = y[i][k];
            for(int j(i + 1); j < n; j++)
                s = s - u[i][j] * x[j][k];
            x[i][k] = s / u[i][i];
        }
    }
}

Matrix LUDecomposer::Solve(Matrix b) const {
    Matrix x(b.NRows(), b.NCols());
    Solve(b, x);
    return x;
}

Matrix LUDecomposer::GetCompactLU() const {
    Matrix m(u.NRows(), u.NCols());
    for(int i(0); i < m.NRows(); i++)
        for(int j(0); j < m.NCols(); j++)
            if(i <= j) m[i][j] = u[i][j];
            else m[i][j]  = l[i][j];
    return m;
}



// Matrix ------------------------------------------------------------------------------------
inline void Matrix::RazmjenaKolona(int a, int b) {
    for(int i(0); i < NRows(); i++) {
        double temp(matrica[i][a]);
        matrica[i][a] = matrica[i][b];
        matrica[i][b] = temp;
    }
}

inline void Matrix::RazmjenaRedova(int a, int b) {
    for(int i(0); i < NCols(); i++) {
        double temp(matrica[a][i]);
        matrica[a][i] = matrica[b][i];
        matrica[b][i] = temp;
    }
}

Matrix::Matrix(int m, int n) {
    if(m <= 0 || n <= 0) throw std::range_error("Bad dimension");
    matrica.resize(m);
    for(int i(0); i < m; i++) {
        matrica[i].resize(n);
        for(int j(0); j < n; j++)
            matrica[i][j] = 0;
    }
}

Matrix::Matrix(const Vector &v) {
    matrica.resize(v.NElems());
    for(unsigned int i(0); i < matrica.size(); i++) {
        matrica[i].resize(1);
        matrica[i][0] = v[i];
    }
}

Matrix::Matrix(std::initializer_list<std::vector<double>> l) {
    if(l.size() <= 0) throw std::range_error("Bad dimension");
    for(auto it(l.begin()); it < l.end(); it++) {
        if(it->size() <= 0) throw std::range_error("Bad dimension");
        if(it->size() != l.begin()->size()) throw std::logic_error("Bad matrix");
    }
    auto it(l.begin());
    matrica.resize(l.size());
    for(unsigned int i(0); i < matrica.size(); i++) {
        matrica[i].resize(it->size());
        for(unsigned int j(0); j < matrica[i].size(); j++)
            matrica[i][j] = (*it)[j];
        it++;
    }
}

inline double &Matrix::operator()(int i, int j) {
    if(i <= 0 || j <= 0 || i > (int)matrica.size() || j > (int)matrica[0].size()) throw std::range_error("Invalid index");
    return matrica[i - 1][j - 1];
}

inline double Matrix::operator()(int i, int j) const {
    if(i <= 0 || j <= 0 || i > (int)matrica.size() || j > (int)matrica[0].size()) throw std::range_error("Invalid index");
    return matrica[i - 1][j - 1];
}

double Matrix::Norm() const {
    double suma(0);
    for(unsigned int i(0); i < matrica.size(); i++)
        for(unsigned int j(0); j < matrica[i].size(); j++)
            suma += matrica[i][j] * matrica[i][j];
    return std::sqrt(suma);
}

double MatrixNorm(const Matrix &m) {
    double suma(0);
    for(int i(0); i < m.NRows(); i++)
        for(int j(0); j < m.NCols(); j++)
            suma += m.matrica[i][j] * m.matrica[i][j];
    return std::sqrt(suma);
}

void Matrix::Print(int width, double eps) const {
    double novi_eps = eps;
    if(eps < 0) novi_eps = GetEpsilon();
    for(unsigned int i(0); i < matrica.size(); i++) {
        for(unsigned int j(0); j < matrica[i].size(); j++) {
            if(std::abs(matrica[i][j]) < novi_eps) std::cout << std::setw(width) << 0;
            else std::cout << std::setw(width) << matrica[i][j];
        }
        std::cout << std::endl;
    }
}

void PrintMatrix(const Matrix &m, int width = 10, double eps = -1) {
    double novi_eps = eps;
    if(eps < 0) novi_eps = m.GetEpsilon();
    for(int i(0); i < m.NRows(); i++) {
        for(int j(0); j < m.NCols(); j++) {
            if(std::abs(m.matrica[i][j]) < novi_eps) std::cout << std::setw(width) << 0;
            else std::cout << std::setw(width) << m.matrica[i][j];
        }
        std::cout << std::endl;
    }
}

Matrix operator +(const Matrix &m1, const Matrix &m2) {
    if((m1.NRows() != m2.NRows()) || (m1.NCols() != m2.NCols())) throw std::domain_error("Incompatible formats");
    Matrix m(m1.NRows(), m1.NCols());
    for(int i(0); i < m1.NRows(); i++)
        for(int j(0); j < m1.NCols(); j++)
            m.matrica[i][j] = m1.matrica[i][j] + m2.matrica[i][j];
    return m;
}

Matrix &Matrix::operator +=(const Matrix &m) {
     if((m.NRows() != NRows()) || (m.NCols() != NCols())) throw std::domain_error("Incompatible formats");
     for(int i(0); i < NRows(); i++)
        for(int j(0); j < NCols(); j++)
            matrica[i][j] += m.matrica[i][j];
    return *this;
}

Matrix operator -(const Matrix &m1, const Matrix &m2) {
    if((m1.NRows() != m2.NRows()) || (m1.NCols() != m2.NCols())) throw std::domain_error("Incompatible formats");
    Matrix m(m1.NRows(), m1.NCols());
    for(int i(0); i < m1.NRows(); i++)
        for(int j(0); j < m1.NCols(); j++)
            m.matrica[i][j] = m1.matrica[i][j] - m2.matrica[i][j];
    return m;
}

Matrix &Matrix::operator -=(const Matrix &m) {
     if((m.NRows() != NRows()) || (m.NCols() != NCols())) throw std::domain_error("Incompatible formats");
     for(int i(0); i < NRows(); i++)
        for(int j(0); j < NCols(); j++)
            matrica[i][j] -= m.matrica[i][j];
    return *this;
}

Matrix operator *(double s, const Matrix &m) {
    Matrix mat(m.NRows(), m.NCols());
     for(int i(0); i < m.NRows(); i++)
        for(int j(0); j < m.NCols(); j++)
            mat[i][j] = m.matrica[i][j] * s;
    return mat;
}

Matrix operator *(const Matrix &m, double s) {
    Matrix mat(m.NRows(), m.NCols());
     for(int i(0); i < m.NRows(); i++)
        for(int j(0); j < m.NCols(); j++)
            mat[i][j] = m.matrica[i][j] * s;
    return mat;
}

Matrix &Matrix::operator *=(double s) {
    for(int i(0); i < NRows(); i++)
        for(int j(0); j < NCols(); j++)
            matrica[i][j] *= s;
    return *this;
}

Matrix operator *(const Matrix &m1, const Matrix &m2) {
    if(m1.NCols() != m2.NRows()) throw std::domain_error("Incompatible formats");
    Matrix m(m1.NRows(), m2.NCols());
    for(int i(0); i < m1.NRows(); i++)
        for(int j(0); j < m2.NCols(); j++)
            for(int k(0); k < m1.NCols(); k++)
                m[i][j] += m1[i][k] * m2[k][j];
    return m;
}

Matrix &Matrix::operator *=(const Matrix &m) {
    if(NCols() != m.NRows()) throw std::domain_error("Incompatible formats");
    *this = *this * m;
    return *this;
}

Vector operator *(const Matrix &m, const Vector &v) {
    if(m.NCols() != v.NElems()) throw std::domain_error("Incompatible formats");
    Matrix m1(v.NElems(), 1);
    
    Vector vek(m.NRows());
    for(int i(0); i < m.NRows(); i++)
            for(int k(0); k < m.NCols(); k++)
                vek[i] += m.matrica[i][k] * v[k];
    return vek;
}

Matrix Transpose(const Matrix &m) {
    Matrix mat(m.NCols(), m.NRows());
    for(int i(0); i < m.NRows(); i++)
        for(int j(0); j < m.NCols(); j++)
            mat.matrica[j][i] = m.matrica[i][j];
    return mat;
}

void Matrix::Transpose() {
    if(NRows() == NCols()) {
        for(int i(0); i < NRows(); i++)
            for(int j(i + 1); j < NCols(); j++) {
                double temp(matrica[i][j]);
                matrica[i][j] = matrica[j][i];
                matrica[j][i] = temp;
            }
    }
    else {
        Matrix mat(NCols(), NRows());
        for(int i(0); i < NRows(); i++)
            for(int j(0); j < NCols(); j++)
                mat.matrica[j][i] = matrica[i][j];
        *this = mat;
    }
}

void Matrix::Chop(double eps) {
    if(eps < 0) eps = GetEpsilon();
    for(int i(0); i < NRows(); i++)
        for(int j(0); j < NCols(); j++)
            if(std::abs(matrica[i][j]) < eps) matrica[i][j] = 0;
}

bool Matrix::EqualTo(const Matrix &m, double eps) const {
    if(m.NCols() != NCols() || m.NRows() != NRows()) return false;
    double novi_eps(eps);
    if(eps < 0) novi_eps = GetEpsilon();
    for(int i(0); i < NRows(); i++)
        for(int j(0); j < NCols(); j++)
            if(std::abs(matrica[i][j] - m.matrica[i][j]) > novi_eps) return false;
    return true;
}

Matrix LeftDiv(Matrix m1, Matrix m2) {
    if(m1.NCols() != m1.NRows()) throw std::domain_error("Divisor matrix is not square");
    if(m2.NRows() != m1.NRows()) throw std::domain_error("Incompatible formats");
    
    int n(m1.NRows()), m(m2.NCols()), p;
    double mi;
    Matrix x(n, m);
    
    for(int k(0); k < n; k++) {
        p = k;
        for(int i(k + 1); i < n; i++)
            if(std::abs(m1.matrica[i][k]) > std::abs(m1.matrica[p][k]))
                p = i;
        if(std::abs(m1.matrica[p][k]) < m1.GetEpsilon()) throw std::domain_error("Divisor matrix is singular");
        if(p != k) {
            m1.RazmjenaRedova(p, k);
            m2.RazmjenaRedova(p, k);
        }
        for(int i(k + 1); i < n; i++) {
            mi = m1.matrica[i][k] / m1.matrica[k][k];
            for(int j(k + 1); j < n; j++)
                m1.matrica[i][j] = m1.matrica[i][j] - mi * m1.matrica[k][j];
            for(int j(0); j < m; j++)
                m2.matrica[i][j] = m2.matrica[i][j] - mi * m2.matrica[k][j];
        }
    }
    
    double s;
    for(int k(0); k < m; k++) {
        for(int i(n - 1); i >= 0; i--) {
            s = m2.matrica[i][k];
            for(int j(i + 1); j < n; j++)
                s = s - m1.matrica[i][j] * x.matrica[j][k];
            x.matrica[i][k] = s / m1.matrica[i][i];
        } 
    }
    
    return x;
}

Vector LeftDiv(Matrix m, Vector v) {
    Matrix pomocna(v.NElems(), 1);
    for(int i(0); i < v.NElems(); i++)
        pomocna[i][0] = v[i];
    pomocna = LeftDiv(m, pomocna);
    for(int i(0); i < v.NElems(); i++)
        v[i] = pomocna[i][0];
    return v;
}

Matrix operator /(const Matrix &m, double s) {
    if(std::abs(s) < m.GetEpsilon()) throw std::domain_error("Division by zero");
    Matrix mat(m);
    for(int i(0); i < mat.NRows(); i++)
        for(int j(0); j < mat.NCols(); j++)
            mat.matrica[i][j] /= s;
    return mat;
}

Matrix &Matrix::operator /=(double s) {
    if(std::abs(s) < GetEpsilon()) throw std::domain_error("Division by zero");
    for(int i(0); i < NRows(); i++)
        for(int j(0); j < NCols(); j++)
            matrica[i][j] /= s;
    return *this;
}

Matrix operator /(Matrix m1, Matrix m2) {
    m1 /= m2;
    return m1;
}

Matrix &Matrix::operator /=(Matrix m) {
    if(m.NCols() != m.NRows()) throw std::domain_error("Divisor matrix is not square");
    if(std::abs(m.Det()) < m.GetEpsilon()) throw std::domain_error("Divisor matrix is singular");
    if(NCols() != m.NRows()) throw std::domain_error("Incompatible formats");
    
    int n(NRows()), mm(m.NCols()), p;
    double mi;
    Matrix x(mm, n);
    
    for(int k(0); k < n; k++) {
        p = k;
        for(int i(k + 1); i < n; i++)
            if(std::abs(matrica[k][i]) > std::abs(matrica[k][p]))
                p = i;
        if(std::abs(matrica[k][p]) < GetEpsilon()) throw std::domain_error("Divisor matrix is singular");
        if(p != k) {
            RazmjenaKolona(p, k);
            m.RazmjenaKolona(p, k);
        }
        for(int i(k + 1); i < n; i++) {
            mi = matrica[k][i] / matrica[k][k];
            for(int j(k + 1); j < n; j++)
                matrica[j][i] = matrica[j][i] - mi * matrica[j][k];
            for(int j(0); j < mm; j++)
                m.matrica[j][i] = m.matrica[j][i] - mi * m.matrica[j][k];
        }
    }
    
    double s;
    for(int k(0); k < mm; k++) {
        for(int i(n - 1); i >= 0; i--) {
            s = m.matrica[k][i];
            for(int j(i + 1); j < n; j++)
                s = s - matrica[j][i] * x.matrica[k][j];
            x.matrica[k][i] = s / matrica[i][i];
        } 
    }
    
    x.Transpose();
    *this = x;
    return *this;
}

double Matrix::Det() const {
    if(NCols() != NRows()) throw std::domain_error("Matrix is not square");
    Matrix m = *this;
    double d(1), mi;
    int n(NRows()), p;
    for(int k(0); k < n; k++) {
        p = k;
        for(int i(k + 1); i < n; i++)
            if(std::abs(m.matrica[i][k]) > std::abs(m.matrica[p][k])) p = i;
        if(std::abs(m.matrica[p][k]) < GetEpsilon()) return 0;
        if(p != k) {
            for(int a(0); a < n; a++) {
                double temp(m.matrica[k][a]);
                m.matrica[k][a] = m.matrica[p][a];
                m.matrica[p][a] = temp;
            }
            d = -1 * d;
        }
        d *= m.matrica[k][k];
        for(int i(k + 1); i < n; i++) {
            mi = m.matrica[i][k] / m.matrica[k][k];
            for(int j(k + 1); j < n; j++)
                m.matrica[i][j] = m.matrica[i][j] - mi * m.matrica[k][j];
        }
    }
    return d;
}

double Det(Matrix m) {
    if(m.NCols() != m.NRows()) throw std::domain_error("Matrix is not square");
   double d(1), mi;
    int n(m.NRows()), p;
    for(int k(0); k < n; k++) {
        p = k;
        for(int i(k + 1); i < n; i++)
            if(std::abs(m.matrica[i][k]) > std::abs(m.matrica[p][k])) p = i;
        if(std::abs(m.matrica[p][k]) < m.GetEpsilon()) return 0;
        if(p != k) {
            for(int a(0); a < n; a++) {
                double temp(m.matrica[k][a]);
                m.matrica[k][a] = m.matrica[p][a];
                m.matrica[p][a] = temp;
            }
            d = -1 * d;
        }
        d *= m.matrica[k][k];
        for(int i(k + 1); i < n; i++) {
            mi = m.matrica[i][k] / m.matrica[k][k];
            for(int j(k + 1); j < n; j++)
                m.matrica[i][j] = m.matrica[i][j] - mi * m.matrica[k][j];
        }
    }
    return d;
}

void Matrix::Invert() {
    if(NCols() != NRows()) throw std::domain_error("Matrix is not square");
    int n(NCols());
    std::vector<double> w(n);
    int p;
    
    for(int k(0); k < n; k++) {
        p = k;
        for(int i(k + 1); i < n; i++)
            if(std::abs(matrica[i][k]) > std::abs(matrica[p][k])) p = i;
        if(std::abs(matrica[p][k]) < GetEpsilon()) throw std::domain_error("Matrix is singular");
        if(p != k) {
            // Razmijeni k-ti i p-ti red u matrici A
            for(int a(0); a < n; a++) {
                double temp(matrica[k][a]);
                matrica[k][a] = matrica[p][a];
                matrica[p][a] = temp;
            }
        }
        w[k] = p;
        double mi(matrica[k][k]);
        matrica[k][k] = 1;
        for(int j(0); j < n; j++)
            matrica[k][j] /=  mi;
        for(int i(0); i < n; i++) {
            if(i != k) {
                mi = matrica[i][k];
                matrica[i][k] = 0;
                for(int j(0); j < n; j++)
                    matrica[i][j] = matrica[i][j] - mi * matrica[k][j];
            }
        }
    }
    
    for(int j(n - 1); j >= 0; j--) {
        p = w[j];
        if(p != j)
            for(int i(0); i < n; i++) {
                double temp(matrica[i][p]);
                matrica[i][p] = matrica[i][j];
                matrica[i][j] = temp;
            }
    }
}

Matrix Inverse(Matrix m) {
    m.Invert();
    return m;
}

void Matrix::ReduceToRREF() {
    int m(NRows()), n(NCols()), k(-1), l(-1), p(0);
    double mi(0), v;
    std::vector<bool> w(n);
    for(int j(0); j < n; j++)
        w[j] = false;
    while(k < m && l < n) {
        l += 1;
        k += 1;
        v = 0;
        while(v < GetEpsilon() && l < n) {
            p = k;
            for(int i(k); i < m; i++)
                if(std::abs(matrica[i][l]) > v) {
                    v = std::abs(matrica[i][l]);
                    p = i;
                }
            if(v < GetEpsilon())
                l += 1;
        }
        if(l < n) {
            w[l] = true;
            if(p != k) {
                // Razmijeni k-ti i p-ti red u matrici
                RazmjenaRedova(p, k);
            }
            mi = matrica[k][l];
            for(int j(l); j < n; j++)
                matrica[k][j] = matrica[k][j] / mi;
            for(int i(0); i < m; i++)
                if(i != k) {
                    mi = matrica[i][l];
                    for(int j(l); j < n; j++)
                        matrica[i][j] = matrica[i][j] - mi * matrica[k][j];
                }
        }
    }
}

Matrix RREF(Matrix m) {
    m.ReduceToRREF();
    return m;
}

int Matrix::Rank() const {
    Matrix m(*this);
    int mm(m.NRows()), n(m.NCols()), k(-1), l(-1), p(0);
    double mi(0), v;
    std::vector<bool> w(n);
    for(int j(0); j < n; j++)
        w[j] = false;
    while(k < mm && l < n) {
        l += 1;
        k += 1;
        v = 0;
        while(v < m.GetEpsilon() && l < n) {
            p = k;
            for(int i(k); i < mm; i++)
                if(std::abs(m.matrica[i][l]) > v) {
                    v = std::abs(m.matrica[i][l]);
                    p = i;
                }
            if(v < m.GetEpsilon())
                l += 1;
        }
        if(l < n) {
            w[l] = true;
            if(p != k) {
                // Razmijeni k-ti i p-ti red u matrici
                for(int a(0); a < n; a++) {
                    double temp(m.matrica[k][a]);
                    m.matrica[k][a] = m.matrica[p][a];
                    m.matrica[p][a] = temp;
                }
            }
            mi = m.matrica[k][l];
            for(int j(l); j < n; j++)
                m.matrica[k][j] /= mi;
            for(int i(0); i < mm; i++)
                if(i != k) {
                    mi = m.matrica[i][l];
                    for(int j(l); j < n; j++)
                        m.matrica[i][j] = m.matrica[i][j] - mi * m.matrica[k][j];
                }
        }
    }
    return k;
}

int Rank(Matrix m) {
    return m.Rank();    
}

// Vector ------------------------------------------------------------------------------------

inline Vector::Vector(int n) {
    if(n <= 0) throw std::range_error("Bad dimension");
    v.resize(n);
    for(unsigned int i(0); i < v.size(); i++) v[i] = 0;
}

Vector::Vector(std::initializer_list<double> l) {
    if(l.size() == 0) throw std::range_error("Bad dimension");
    v.resize(l.size());
    auto it(l.begin());
    for(unsigned int i(0); i < v.size(); i++) {
        v[i] = *it;
        it++;
    }
}

inline double &Vector::operator()(int i) {
    if(i <= 0 || (unsigned)i > v.size()) throw std::range_error("Invalid index");
    else return v[i - 1];
}

inline double Vector::operator()(int i) const {
    if(i <= 0 || (unsigned)i > v.size()) throw std::range_error("Invalid index");
    else return v[i - 1];
}

double Vector::Norm() const {
    double suma(0);
    for(unsigned int i(0); i < v.size(); i++) suma += v[i] * v[i];
    return std::sqrt(suma);
}

double VectorNorm(const Vector &v) {
    double suma(0);
    for(int i(0); i < v.NElems(); i++) suma += v[i] * v[i];
    return std::sqrt(suma);
}

void Vector::Print(char separator, double eps) const {
    double novi_eps(eps);
    if(eps < 0) novi_eps = GetEpsilon();
    unsigned int i(0);
    for(i = 0; i < v.size() - 1; i++)
        if(std::abs(v[i]) < novi_eps) std::cout << 0 << separator;
        else std::cout << v[i] << separator;
    if(std::abs(v[i]) < novi_eps) std::cout << 0;
    else std::cout << v[i];
    if(separator == '\n') std::cout << separator;
}

void PrintVector(const Vector &v, char separator = '\n', double eps = -1) {
    double novi_eps(eps);
    if(eps < 0) novi_eps = v.GetEpsilon();
    int i(0);
    for(i = 0; i < v.NElems() - 1; i++)
        if(std::abs(v[i]) < novi_eps) std::cout << 0 << separator;
        else std::cout << v[i] << separator;
    if(std::abs(v[i]) < novi_eps) std::cout << 0;
    else std::cout << v[i];
    if(separator == '\n') std::cout << separator;
}

Vector operator +(const Vector &v1, const Vector &v2) {
    if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3(v1.NElems());
    for(int i(0); i < v3.NElems(); i++) v3[i] = v1[i] + v2[i];
    return v3;
}

Vector &Vector::operator +=(const Vector &v) {
    if(NElems() != v.NElems()) throw std::domain_error("Incompatible formats");
    for(int i(0); i < v.NElems(); i++) this->v[i] += v[i];
    return *this;
}

Vector operator -(const Vector &v1, const Vector &v2) {
    if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3(v1.NElems());
    for(int i(0); i < v3.NElems(); i++) v3[i] = v1[i] - v2[i];
    return v3;
}

Vector &Vector::operator -=(const Vector &v) {
    if(NElems() != v.NElems()) throw std::domain_error("Incompatible formats");
    for(int i(0); i < v.NElems(); i++) this->v[i] -= v[i];
    return *this;
}

Vector operator *(double s, const Vector &v) {
    Vector v1(v.NElems());
    for(int i(0); i < v.NElems(); i++) v1[i] = v[i];
    for(int i(0); i < v.NElems(); i++) v1[i] *= s;
    return v1;
}

Vector operator *(const Vector &v, double s) {
    Vector v1(v.NElems());
    for(int i(0); i < v.NElems(); i++) v1[i] = v[i];
    for(int i(0); i < v.NElems(); i++) v1[i] *= s;
    return v1;
}

Vector &Vector::operator *=(double s) {
    for(int i(0); i < NElems(); i++) v[i] *= s;
    return *this;
}

double operator *(const Vector &v1, const Vector &v2) {
    if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
    double suma(0);
    for(int i(0); i < v1.NElems(); i++) suma += v1[i] * v2[i];
    return suma;
}

Vector operator /(const Vector &v, double s) {
    if(std::abs(s) < v.GetEpsilon()) throw std::domain_error("Division by zero");
    Vector v1(v.NElems());
    for(int i(0); i < v.NElems(); i++) v1[i] = v[i];
    for(int i(0); i < v.NElems(); i++) v1[i] /= s;
    return v1;
}

Vector &Vector::operator /=(double s) {
    if(std::abs(s) < GetEpsilon()) throw std::domain_error("Division by zero");
    for(int i(0); i < NElems(); i++) v[i] /= s;
    return *this;
}

void Vector::Chop(double eps) {
    double novi_eps = eps;
    if(eps < 0) novi_eps = GetEpsilon();
    for(unsigned int i(0); i < v.size(); i++)
        if(std::abs(v[i]) < novi_eps) v[i] = 0;
}

bool Vector::EqualTo(const Vector &v, double eps) const {
    if((unsigned)v.NElems() != this->v.size()) return false;
    double novi_eps = eps;
    if(eps < 0) novi_eps = GetEpsilon();
    for(unsigned int i(0); i < this->v.size(); i++)
        if(std::abs(v[i] - this->v[i]) > novi_eps) return false;
    return true;
}

// -------------------------------------------------------------------------------------------

int main () {
    
    // TESTIRANJE KLASA VECTOR I MATRIX
    {
    //TESTIRANJE KLASE VECTOR
        // Testiranje obicnog konstrutktora i metode Print
        Vector v1(5);
        v1.Print(',');
        std::cout << std::endl;
        try {
            Vector v3(0);
        }
        catch(std::range_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje sekvencijskog konstruktora i funkcije PrintVector
        Vector v2{1, 2, 3, 4, 5};
        v2.Print(',');
        std::cout << std::endl;
        try {
            Vector v3{};
        }
        catch(std::range_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije NElems
        std::cout << "Broj elemenata v1: " << v1.NElems() << std::endl;
        std::cout << "Broj elemenata v2: " << v2.NElems() << std::endl;
        
        // Testiranje operatora []
        std::cout << "Vrijednost v2[3] je: " << v2[3] << std::endl;
        v1[3] = 10;
        std::cout << "Nova vrijednost od v2[3] je: " << v2[3] << std::endl;
        
        // Testiranje operatora ()
        std::cout << "Vrijednost v2(3) je: " << v2(3) << std::endl;
        v2(3) = 10;
        std::cout << "Nova vrijednost v2(3) je: " << v2(3) << std::endl;
        try {
            std::cout << v2(0);
        }
        catch(std::range_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje metode Norm i prijateljske funkcije
        std::cout << "Euklidska norma vektora v2 je: " << v2.Norm() << std::endl;
        std::cout << "Euklidska norma vektora v2 je (friend funkcija): " << VectorNorm(v2) << std::endl;
        
        // Testiranje metode GetEpsion
        std::cout << "Euklidska norma pomnozena sa masinskim epsilonom: " << v2.GetEpsilon() << std::endl;
        
        // Testiranje sabiranja vektora (+)
        Vector v3{10, 15, 20, 25, 30};
        Vector v4(v2 + v3);
        std::cout << "Zbir vektora (+) v2 i v3 je vektor: ";
        v4.Print(',');
        std::cout << std::endl;
        try {
            Vector v5(5), v6(10);
            Vector v7(v5 + v6); // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje sabiranja vektora (+=)
        v2 += v3;
        std::cout << "Zbir vektora (+=) v2 i v3 je vektor v2: ";
        v4.Print(',');
        std::cout << std::endl;
        try {
            Vector v5(5), v6(10);
            Vector v7(v5 + v6); // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje oduzimanja vektora (-)
        std::cout << "Razlika vektora (-) v2 i v3 je vektor: ";
        v4 = v2 - v3;
        v4.Print(',');
        std::cout << std::endl;
        try {
            Vector v5(5), v6(10);
            Vector v7(v5 - v6); // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje oduzimanja vektora (-=)
        std::cout << "Razlika vektora (-=) v2 i v3 je vektor v2: ";
        v2 -= v3;
        v2.Print(',');
        std::cout << std::endl;
        try {
            Vector v5(5), v6(10);
            Vector v7(v5 - v6); // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje mnozenja skalara i vektora (*)
        Vector v5{10, 20, 30, 40, 50};
        v5 = 5 * v5;
        std::cout << "Vektor v5 ponozen sa 5 (*) je: ";
        v5.Print(',');
        std::cout << std::endl;
        
        // Testiranje mnozenja vektora i skalara
        v5 = v5 * 5;
        std::cout << "Vektor v5 ponozen sa 5 (*) je: ";
        v5.Print(',');
        std::cout << std::endl;
        
        // Testiranje mnozenja vektora i skalara
        v5 *= 5;
        std::cout << "Vektor v5 ponozen sa 5 (*=) je: ";
        v5.Print(',');
        std::cout << std::endl;
        
        // Testiranje skalarnog proizvoda vektora
        Vector v6{10, 10, 10, 10, 10}, v7{5, 5, 5, 5, 5};
        std::cout << "Skalarni proizvod vektora v6 i v7 je: ";
        std::cout << v6 * v7 << std::endl;
        try {
            Vector v8(5), v9(10);
            std::cout << v8 * v9; // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje dijeljenja vektora skalarom (/)
        v6 = v6 / 2;
        std::cout << "Vektor v6 podijeljen (/) sa 2 je: ";
        v6.Print(',');
        std::cout << std::endl;
        
        // Testiranje dijeljenja vektora skalarom (/=)
        v6 /= 2;
        std::cout << "Vektor v6 podijeljen (/=) sa 2 je: ";
        v6.Print(',');
        std::cout << std::endl;
        try {
            Vector v10(5);
            v10 /= 0; // Treba baciti izuzetak zbog dijeljenja s nulom
        }
        catch(std::range_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        std::cout << std::endl;

    // TESTIRANJE KLASE MATRIX
        // Testiranje konstruktora sa dva parametra i metode Print
        Matrix m1(3,3);
        m1.Print();
        try {
            Matrix m2(0, 0); // Treba baciti izuzetak zbog neregularnih dimenzija
        }
        catch(std::range_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje konstruktora sa vektorom i funkcije PrintMatrix
        Matrix m2(v7);
        PrintMatrix(m2);
        
        // Testiranje sekvencijskog konstruktora
        Matrix m3{{5, 5, 5}, {10, 10, 10}, {15, 15, 15}};
        m3.Print();
        
        // Testiranje metoda za velicinu matrice
        std::cout << "Matrica m1 je dimenzija: " << m1.NRows() << "x" << m1.NCols() << std::endl;
        std::cout << "Matrica m2 je dimenzija: " << m2.NRows() << "x" << m2.NCols() << std::endl;
        std::cout << "Matrica m3 je dimenzija: " << m3.NRows() << "x" << m3.NCols() << std::endl;
        
        // Testiranje operatora []
        std::cout << "Vrijednost m3[1][1] je: " << m3[1][1] << std::endl;
        m3[1][1] = 100;
        std::cout << "Nova vrijednost m3[1][1] je: " << m3[1][1] << std::endl;
        
        // Testiranje operatora ()
        std::cout << "Vrijednost m3(1,1) je: " << m3(1,1) << std::endl;
        m3(1,1) = 10;
        std::cout << "Nova vrijednost m3(1,1) je: " << m3(1,1) << std::endl;
        
        // Testiranje metode Norm i funkcije MatrixNorm
        std::cout << "Frobeniusova norma matrice m3 je: " << m3.Norm() << std::endl;
        std::cout << "Frobeniusova norma matrice m3 je: " << MatrixNorm(m3) << std::endl;
        
        // Testiranje metode GetEpsilon
        std::cout << "Frobeniusova norma pomnozena sa masinskim epsilonom: " << m3.GetEpsilon() << std::endl;
        
        // Testiranje sabiranja matrica (+)
        Matrix m4{{5, 5, 5}, {10, 10, 10}, {15, 15, 15}};
        Matrix m5{{5, 5, 5}, {10, 10, 10}, {15, 15, 15}};
        std::cout << "Zbir matrica m4 i m5 je:" << std::endl;
        PrintMatrix(m4 + m5);
        try {
            Matrix m6 (m2 + m4); // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje sabiranja matrica (+=)
        std::cout << "Zbir matrica m4 i m5 je matrica m4:" << std::endl;
        m4 += m5;
        PrintMatrix(m4);
        try {
            m2 += m4; // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje oduzimanja matrica (-)
        std::cout << "Razlika matrica m4 i m5 je:" << std::endl;
        PrintMatrix(m4 - m5);
        try {
            Matrix m6 (m2 - m4); // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje oduzimanja matrica (-=)
        std::cout << "Razlika matrica m4 i m5 je matrica m4:" << std::endl;
        m4 -= m5;
        PrintMatrix(m4);
        try {
            m2 -= m4; // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje mnozenja skalara matricom
        std::cout << "Broj 5 pomnozen matricom m5 je matrica m4:" << std::endl;
        m4 = 5 * m5;
        PrintMatrix(m4);
        
        // Testiranje mnozenja matrice skalarom (*)
        std::cout << "Matrica m4 pomnozena sa 5 je matrica m4:" << std::endl;
        m4 = m4 * 5;
        PrintMatrix(m4);
        
        // Testiranje mnozenja matrice skalarom (*=)
        std::cout << "Matrica m4 pomnozena sa 5 je matrica m4:" << std::endl;
        m4 *= 5;
        PrintMatrix(m4);
        
        // Testiranje mnozenja matrica (*=)
        Matrix m6{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}, {1, 1, 1}}, m7{{2, 2, 2}, {2, 2, 2}, {2, 2, 2}};
        std::cout << "Proizvod matrica m6 i m7 je:" << std::endl;
        PrintMatrix(m6 * m7);
        try {
            m2 = m2 * m4; // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje mnozenja matrica (*=)
        std::cout << "Proizvod matrica m6 i m7 je matrica m6:" << std::endl;
        m6 *= m7;
        PrintMatrix(m6);
        try {
            m2 *= m4; // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje mnozenja matrice i vektora
        Matrix m8{{5, 5, 5}, {5, 5, 5}, {5, 5, 5}, {5, 5, 5}};
        Vector v10{5, 5, 5};
        std::cout << "Matrica m8 pomnozena vektorom 10 je vektor: ";
        PrintVector(m8 * v10, ',');
        std::cout << std::endl;
        try {
            Matrix m9{{5, 5, 5}, {5, 5, 5}, {5, 5, 5}, {5, 5, 5}};
            Vector v11{5, 5, 5, 5};
            PrintVector(m9 * v11); // Treba baciti izuzetak zbog razlicitih dimenzija
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje transponovanja matrice (kvadratna)
        std::cout << "Matrica m9 nakon transponovanja je:" << std::endl;
        Matrix m9{{1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}};
        PrintMatrix(Transpose(m9));
        
        // Testiranje transponovanja matrice (nije kvadratna)
        std::cout << "Matrica m10 nakon transponovanja je:" << std::endl;
        Matrix m10{{1, 2, 3}, {1, 2, 3}, {1, 2, 3}, {1, 2, 3}};
        PrintMatrix(Transpose(m10));
        
        // Testiranje transponovanja matrice funkcijom clanicom (kvadratna)
        m9.Transpose();
        std::cout << "Matrica m9 nakon transponovanja je:" << std::endl;
        PrintMatrix(m9);
        
        // Testiranje transponovanja matrice funkcijom clanicom (nije kvadratna)
        m10.Transpose();
        std::cout << "Matrica m10 nakon transponovanja je:" << std::endl;
        PrintMatrix(m10);
    }
    
    
    // TESTIRANJE KLASE VEKTOR
    {
        // Testiranje funkcije Chop
        Vector v1{1, 2, 3, 0.000001, 0.00000002, 0.1, 5};
        v1.Chop(0.001);
        std::cout << "Vektor v1 nakon funkcije Chop: ";
        v1.Print(' ');
        std::cout << std::endl;
        
        // Testiranje funckije EqualTo
        Vector v2{1, 2, 3, 0, 0, 0.1, 5};
        if(v1.EqualTo(v2)) std::cout << "Vektori v1 i v2 su isti." << std::endl;
        else std::cout << "Vektori v1 i v2 nisu isti." << std::endl;
        v2[0] = 10;
        if(v1.EqualTo(v2)) std::cout << "Vektori v1 i v2 su isti." << std::endl;
        else std::cout << "Vektori v1 i v2 nisu isti." << std::endl;
        Vector v3(5);
        if(v2.EqualTo(v3)) std::cout << "Vektori v2 i v3 su isti." << std::endl;
        else std::cout << "Vektori v2 i v3 nisu isti." << std::endl;
    }
    
    // TESTIRANJE KLASE MATRICA
    {
        // Testiranje funkcije Chop
        Matrix m1{{5, 5, 0.00005}, {10, 10, 0.000000025}, {15, 15, 15}};
        m1.Chop(0.0001);
        std::cout << "Matrica m1 nakon funkcije Chop:" << std::endl;
        m1.Print(10);
        
        // Testiranje funkcije EqualTo
        Matrix m2{{5, 5, 0}, {10, 10, 0}, {15, 15, 15}};
        if(m1.EqualTo(m2)) std::cout << "Matrice m1 i m2 su iste." << std::endl;
        else std::cout << "Matrice m1 i m2 nisu iste." << std::endl;
        m2[0][0] = 10;
        if(m1.EqualTo(m2)) std::cout << "Matrice m1 i m2 su iste." << std::endl;
        else std::cout << "Matrice m1 i m2 nisu iste." << std::endl;
        Matrix m3(2, 2);
        if(m2.EqualTo(m3)) std::cout << "Matrice m2 i m3 su iste." << std::endl;
        else std::cout << "Matrice m2 i m3 nisu iste." << std::endl;
        
        // Testiranje funkcije Det (friend)
        Matrix m4{{1, 1}, {2, 1}};
        std::cout << "Determinanta matrice m4 je: " << Det(m4) << std::endl;
        try {
            // Treba baciti izuzetak jer nije kvadratna
            Matrix m5{{1, 1}, {2, 1}, {2, 3}}; 
            std::cout << Det(m5);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije Det (članica)
        Matrix m5{{1, 1}, {2, 1}};
        std::cout << "Determinanta matrice m5 je: " << m5.Det() << std::endl;
        try {
            // Treba baciti izuzetak jer nije kvadratna
            Matrix m6{{1, 1}, {2, 1}, {2, 3}}; 
            std::cout << m6.Det();
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije Invert (članica)
        Matrix m6{{1, 3, 1}, {4, 5, 6}, {7, 8, 9}};
        std::cout << "Inverzna matrica matrice m6 je matrica:\n";
        m6.Invert();
        m6.Print(10);
        try {
            // Treba baciti izuzetak jer je determinanta 0
            Matrix m7{{1, 1}, {1, 1}};
            m7.Invert();
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak jer nije kvadratna
            Matrix m7{{1, 3, 1}, {4, 5, 6}};
            m7.Invert();
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije Invert (friend)
        Matrix m7{{1, 3, 1}, {4, 5, 6}, {7, 8, 9}};
        std::cout << "Inverzna matrica matrice m6 je matrica:\n";
        Inverse(m7).Print();
        try {
            // Treba baciti izuzetak jer je determinanta 0
            Matrix m8{{1, 1}, {1, 1}};
            Inverse(m8);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak jer nije kvadratna
            Matrix m8{{1, 3, 1}, {4, 5, 6}};
            Inverse(m8);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije LeftDiv
        Matrix m8{{1, 3, 1}, {4, 5, 6}, {7, 8, 9}};
        Matrix m9{{1, 3, 1}, {4, 5, 6}, {7, 8, 9}};
        std::cout << "Rezultat lijevog dijeljenja matrica m8 i m9 je matrica:\n";
        LeftDiv(m8, m9).Print();
        try {
            // Treba baciti izuzetak jer je determinanta 0
            Matrix m10{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
            LeftDiv(m10, m9);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak jer nije kvadratna
            Matrix m10{{1, 3}, {4, 5}, {2, 3}};
            LeftDiv(m10, m9);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak jer su dimenzije različite
            Matrix m10{{1, 3}, {4, 5}};
            LeftDiv(m10, m9);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje dijeljenja matrice brojem (/)
        Matrix m10{{10, 10, 10}, {10, 10, 10}, {10, 10, 10}};
        std::cout << "Matrica m10 podijeljena sa 5 je:\n";
        m10 = m10 / 5;
        m10.Print();
        
        // Testiranje dijeljenja matrice brojem (/=)
        std::cout << "Matrica m10 podijeljena sa 5 je:\n";
        m10 /= 5;
        m10.Print();
        
        // Testiranje desnog matričnog dijeljenja
        Matrix m11{{1, 3, 1}, {4, 5, 6}, {7, 8, 9}};
        Matrix m12{{1, 3, 1}, {4, 5, 6}, {7, 8, 9}};
        std::cout << "Rezultat lijevog dijeljenja matrica m8 i m9 je matrica:\n";
        LeftDiv(m12, m11).Print();
        try {
            // Treba baciti izuzetak jer je determinanta 0
            Matrix m13{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
            LeftDiv(m13, m12);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak jer nije kvadratna
            Matrix m13{{1, 3}, {4, 5}, {2, 3}};
            LeftDiv(m13, m12);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak jer su dimenzije različite
            Matrix m13{{1, 3}, {4, 5}};
            LeftDiv(m13, m12);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije ReduceToRREF (članica)
        Matrix m13{{5, 5, 5}, {10.5, 10.5, 10.5}, {20, 20, 20}};
        m13.ReduceToRREF();
        m13.Print();
        
        // Testiranje funkcije ReduceToRREF (friend)
        Matrix m14{{5, 5, 5}, {10.5, 10.5, 10.5}, {20, 20, 20}};
        m14 = RREF(m14);
        m14.Print();
        
        // Testiranje funkcije Rank (članica)
        std::cout << "Rang matrice m14 je: " << m14.Rank() << std::endl;
      
        // Testiranje funkcije Rank (friend)
        std::cout << "Rang matrice m14 je: " << Rank(m14) << std::endl;
    }
        
    // TESTIRANJE KLASE LU DECOMPOSER
    {
        // Testiranje konstruktora
        Matrix m1{{11, 9, 24, 2}, {1, 5, 2, 6}, {3, 17, 18, 1}, {2, 5, 7, 1}};
        LUDecomposer lud1(m1);
        try {
            // Treba baciti izuzetak jer matrica nije kvadratna
            Matrix m2{{1, 1}, {1, 1}, {10, 2}};
            LUDecomposer lud2(m2);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak jer je matrica singularna
            Matrix m2{{1, 1}, {1, 1}};
            LUDecomposer lud2(m2);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije GetL
        std::cout << "Matrica l matrice m15 je:\n";
        lud1.GetL().Print();
        
        // Testiranje funkcije GetU
        std::cout << "Matrica u matrice m15 je:\n";
        lud1.GetU().Print();
        
        // Testiranje funkcije Solve (sa dva vektora)
        Vector b{1, 2, 3, 4}, x(4);
        lud1.Solve(b, x);
        std::cout << "Rjesenja sistema su: ";
        PrintVector(x, ' ');
        std::cout << std::endl;
        Vector c{10, -11, -14, 17};
        lud1.Solve(c, x);
        std::cout << "Rjesenja sistema su: ";
        PrintVector(x, ' ');
        std::cout << std::endl;
        try {
            // Treba baciti izuzetak zbog pogrešne dimenzije vektora d
            Vector d{1, 2, 3, 4, 5}, y(4);
            lud1.Solve(d, y);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak zbog pogrešne dimenzije vektora x
            Vector d{1, 2, 3, 4}, y(5);
            lud1.Solve(d, y);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije Solve (sa jednim vektorom)
        std::cout << "Rjesenja sistema su: ";
        PrintVector(lud1.Solve(b), ' ');
        std::cout << std::endl;
        std::cout << "Rjesenja sistema su: ";
        PrintVector(lud1.Solve(c), ' ');
        std::cout << std::endl;
        try {
            // Treba baciti izuzetak zbog pogrešne dimenzije vektora d
            Vector d{1, 2, 3, 4, 5};
            lud1.Solve(d);
        }
        catch(std::domain_error izuzetak) {
            std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije Solve (sa dvije matrice)
        
        Matrix m2{{1, 5, 1, 1}, {2, 2, 2, 2}, {2, 2, 2, 2}, {2, 2, 2, 2}};
        Matrix matx2(4, 4);
        std::cout << "Rjesenja sistema su:\n";
        lud1.Solve(m2, matx2);
        matx2.Print();
        Matrix m3{{1, 5, 5, 13}, {21, 5, 23, 2}, {2, 5, 2, 2}, {2, 5, 2, 2}};
        Matrix matx3(4, 4);
        std::cout << "Rjesenja sistema su:\n";
        lud1.Solve(m3, matx3);
        matx3.Print();
        try {
            // Treba baciti izuzetak zbog pogrešne dimenzije matrice m4
            Matrix m4{{1, 5, 1, 1}, {2, 2, 2, 2}, {2, 2, 2, 2}};
            lud1.Solve(m4, matx2);
        }
        catch(std::domain_error izuzetak) {
          std::cout << izuzetak.what() << std::endl;
        }
        try {
            // Treba baciti izuzetak zbog pogrešne dimenzije matrice matx4
            Matrix m4{{1, 5, 1, 1}, {2, 2, 2, 2}, {2, 2, 2, 2}, {2, 3, 4, 5}}, matx4(2, 2);
            lud1.Solve(m4, matx4);
        }
        catch(std::domain_error izuzetak) {
          std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funkcije Solve (sa jednom matricom)
        std::cout << "Rjesenja sistema su: ";
        PrintMatrix(lud1.Solve(m2));
        std::cout << std::endl;
        std::cout << "Rjesenja sistema su: ";
        PrintMatrix(lud1.Solve(m3));
        std::cout << std::endl;
        try {
            // Treba baciti izuzetak zbog pogrešne dimenzije matrice m4
            Matrix m4{{1, 5, 1, 1}, {2, 2, 2, 2}, {2, 2, 2, 2}};
            lud1.Solve(m4, matx2);
        }
        catch(std::domain_error izuzetak) {
          std::cout << izuzetak.what() << std::endl;
        }
        
        // Testiranje funckije GetCompactLU
        std::cout << "Matrica U:\n";
        lud1.GetU().Print();
        std::cout << "Matrica L:\n";
        lud1.GetL().Print();
        std::cout << "Rezultat funkcije GetCompactLU:\n";
        PrintMatrix(lud1.GetCompactLU());
        
    }
    
    Matrix B{{2,3,5},{2,3,7},{4,1,8}};
    LUDecomposer lu2(B);
    Vector x{1,2,4};
    Vector rez=B*x;
    lu2.Solve(rez,x);
    x.Print();
    std::cout<<x.EqualTo(LeftDiv(B, rez));
    std::cout << std::endl;
    Matrix C{{1, 2},{2, 1}};
    Vector c{1, 2};
	return 0;
}