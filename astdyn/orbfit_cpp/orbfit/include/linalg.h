#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace orbfit {

// ─────────────────────────────────────────────────────────────
//  3-vector
// ─────────────────────────────────────────────────────────────
struct Vec3 {
    double x, y, z;

    Vec3(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s)      const { return {x*s, y*s, z*s}; }
    Vec3 operator/(double s)      const { return {x/s, y/s, z/s}; }
    Vec3& operator+=(const Vec3& o) { x+=o.x; y+=o.y; z+=o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x-=o.x; y-=o.y; z-=o.z; return *this; }

    double dot(const Vec3& o)  const { return x*o.x + y*o.y + z*o.z; }
    Vec3   cross(const Vec3& o) const {
        return { y*o.z - z*o.y,
                 z*o.x - x*o.z,
                 x*o.y - y*o.x };
    }
    double norm()  const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
    Vec3   unit()  const { double n=norm(); return (n>0)?(*this/n):Vec3{}; }

    double& operator[](int i) {
        if(i==0) return x; if(i==1) return y; return z;
    }
    double operator[](int i) const {
        if(i==0) return x; if(i==1) return y; return z;
    }
};

inline Vec3 operator*(double s, const Vec3& v) { return v*s; }

// ─────────────────────────────────────────────────────────────
//  Dense matrix (row-major)
// ─────────────────────────────────────────────────────────────
class Matrix {
public:
    int rows, cols;
    std::vector<double> data;

    Matrix() : rows(0), cols(0) {}
    Matrix(int r, int c, double val=0.0) : rows(r), cols(c), data(r*c, val) {}

    double& operator()(int i, int j)       { return data[i*cols+j]; }
    double  operator()(int i, int j) const { return data[i*cols+j]; }

    Matrix operator+(const Matrix& o) const {
        Matrix r(rows,cols);
        for(int i=0;i<rows*cols;i++) r.data[i]=data[i]+o.data[i];
        return r;
    }
    Matrix operator-(const Matrix& o) const {
        Matrix r(rows,cols);
        for(int i=0;i<rows*cols;i++) r.data[i]=data[i]-o.data[i];
        return r;
    }
    Matrix operator*(double s) const {
        Matrix r(rows,cols);
        for(int i=0;i<rows*cols;i++) r.data[i]=data[i]*s;
        return r;
    }
    Matrix operator*(const Matrix& B) const {
        if(cols!=B.rows) throw std::runtime_error("Matrix size mismatch in multiply");
        Matrix C(rows,B.cols,0.0);
        for(int i=0;i<rows;i++)
            for(int k=0;k<cols;k++)
                for(int j=0;j<B.cols;j++)
                    C(i,j)+=(*this)(i,k)*B(k,j);
        return C;
    }
    Matrix transpose() const {
        Matrix T(cols,rows);
        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                T(j,i)=(*this)(i,j);
        return T;
    }

    static Matrix identity(int n) {
        Matrix I(n,n,0.0);
        for(int i=0;i<n;i++) I(i,i)=1.0;
        return I;
    }

    void print(const char* name="") const {
        std::cout << name << " (" << rows << "x" << cols << "):\n";
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++)
                std::cout << std::setw(14) << std::scientific << std::setprecision(6) << (*this)(i,j);
            std::cout << "\n";
        }
    }
};

// ─────────────────────────────────────────────────────────────
//  Cholesky decomposition  (A = L L^T,  A symmetric positive definite)
//  Returns lower-triangular L. Throws if not SPD.
// ─────────────────────────────────────────────────────────────
inline Matrix cholesky(const Matrix& A) {
    int n = A.rows;
    Matrix L(n,n,0.0);
    for(int i=0;i<n;i++){
        for(int j=0;j<=i;j++){
            double s=A(i,j);
            for(int k=0;k<j;k++) s-=L(i,k)*L(j,k);
            if(i==j){
                if(s<=0.0) throw std::runtime_error("Cholesky: matrix not positive definite");
                L(i,j)=std::sqrt(s);
            } else {
                L(i,j)=s/L(j,j);
            }
        }
    }
    return L;
}

// ─────────────────────────────────────────────────────────────
//  LU decomposition with partial pivoting
//  Solves A x = b.
// ─────────────────────────────────────────────────────────────
inline std::vector<double> solve_linear(Matrix A, std::vector<double> b) {
    int n = A.rows;
    std::vector<int> piv(n);
    for(int i=0;i<n;i++) piv[i]=i;

    for(int k=0;k<n;k++){
        // find pivot
        int p=k;
        double mx=std::abs(A(k,k));
        for(int i=k+1;i<n;i++)
            if(std::abs(A(i,k))>mx){ mx=std::abs(A(i,k)); p=i; }
        if(mx<1e-15) throw std::runtime_error("Singular matrix in solve_linear");
        std::swap(piv[k],piv[p]);
        for(int j=0;j<n;j++) std::swap(A(k,j),A(p,j));
        std::swap(b[k],b[p]);
        // eliminate
        for(int i=k+1;i<n;i++){
            double f=A(i,k)/A(k,k);
            for(int j=k;j<n;j++) A(i,j)-=f*A(k,j);
            b[i]-=f*b[k];
        }
    }
    // back-substitution
    std::vector<double> x(n);
    for(int i=n-1;i>=0;i--){
        x[i]=b[i];
        for(int j=i+1;j<n;j++) x[i]-=A(i,j)*x[j];
        x[i]/=A(i,i);
    }
    return x;
}

// ─────────────────────────────────────────────────────────────
//  Matrix inverse via LU (small matrices ≤ 6x6)
// ─────────────────────────────────────────────────────────────
inline Matrix invert(const Matrix& A) {
    int n=A.rows;
    Matrix inv(n,n);
    for(int j=0;j<n;j++){
        std::vector<double> e(n,0.0); e[j]=1.0;
        auto col=solve_linear(A,e);
        for(int i=0;i<n;i++) inv(i,j)=col[i];
    }
    return inv;
}

} // namespace orbfit
