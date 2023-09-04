#pragma once
#include<cmath>
#include<iostream>
#include<cassert>

// #define double floaat

template <int n> struct vec{ // 向量
    double a[n + 1] = {0};
    double& operator[](int i)      {assert(i >= 0 && i <= n); return a[i];}
    double  operator[](int i) const{assert(i >= 0 && i <= n); return a[i];}
    double norm2() const{return *this * *this;}
    double norm() const{return sqrt(norm2());}
    vec<n> normalize() const{return (*this) / norm();}
};

template <int n> vec<n> operator+ (const vec<n> &x, const vec<n> &y){
    vec<n> res = x;
    for(int i = 1; i <= n; i++)
        res[i] += y[i];
    return res;
}

template <int n> vec<n> operator- (const vec<n> & x, const vec<n> &y){
    vec<n> res = x;
    for(int i = 1; i <= n; i++)
        res[i] -= y[i];
    return res;
}

template <int n> vec<n> operator- (const vec<n> & x){
    vec<n> res = x;
    for(int i = 1; i <= n; i++)
        res[i] = -res[i];
    return res;
}

template <int n> double operator* (const vec<n> &x, const vec<n> &y){
    double ans = 0;
    for(int i = 1; i <= n; i++)
        ans += x[i] * y[i];
    return ans;
}

template <int n> vec<n> operator* (const double &x, const vec<n> &y){
    vec<n> res = y;
    for(int i = 1; i <= n; i++)
        res[i] *= x;
    return res;
}

template <int n> vec<n> operator* (const vec<n> &x, const double &y){
    vec<n> res = x;
    for(int i = 1; i <= n; i++)
        res[i] *= y;
    return res;
}

template <int n> vec<n> operator/ (const vec<n> &x, const double y){
    vec<n> res = x;
    for(int i = 1; i <= n; i++)
        res[i] /= y;
    return res;
}

template <int n> std::ostream& operator<< (std::ostream& out, const vec<n>& x){
    for(int i = 1; i <= n; i++)
        out << x[i] <<" ";
    return out;
}

template <int n> bool operator!= (const vec<n> &x, const vec<n>& y){
    for(int i = 1; i <= n; i++)
        if(x[i] != y[i])
            return true;
    return false;
}

template <int n> struct dt;

template <int nro, int nco> struct mat{ // 矩阵
    vec<nco> row[nro] = {{}};

    vec<nro>& operator[](int i)      {assert(i >= 0 && i <= nro); return row[i];}
    vec<nro>  operator[](int i) const{assert(i >= 0 && i <= nro); return row[i];}

    vec<nro> col(const int id)const{
        assert(id >= 0 && id <= nco);
        vec<nro> res;
        for(int i = 1; i <= nro; i++)
            res[i] = row[i][id];
        return res;
    }

    void set_col(const int id, const vec<nro> &x){
        assert(id >= 0 && id  <= nco);
        for(int i = 1; i <= nro; i++)
            row[i][id] = x[i];
    }
    
    static mat<nro, nco> identity(){
        mat<nro, nco> res;
        for(int i = 0; i <= nro; i++)
            for(int j = 0; j <= nco; j++)
                res[i][j] = (i == j);
        return res;
    }

    double det() const{
        return dt<nco>::det(*this);
    }

    mat<nro - 1, nco - 1> del(const int x, const int y){
        mat<nro - 1, nco - 1> res;
        for(int i = 1; i <= nro; i++)
            for(int j = 1; i <= nco; i++)
                res[i][j] = row[i < x ? i : i + 1][j < y ? j : j + 1];
        return res;
    }

    double cofactor(const int x, const int y) const{
        return del(x, y).det() * ((x + y) & 1 ? 1 : -1);
    }

    mat<nro, nco> adjugate() const{
        mat<nro, nco> res;
        for(int i = 1; i <= nro; i++)
            for(int j = 1; j <= nco; j++)
                res[i][j] = cofactor(i, j);
        return res;
    }

    mat<nro, nco> transpose() const{
        mat<nco, nro> res;
        for(int i = 1; i <= nco; i++)
            res[i] = this->col(i);
        return res;
    }

};

template<int nro, int nco> vec<nro> operator*(const mat<nro, nco>& x, const vec<nco>& y){
    vec<nro> res;
    for(int i = 1; i <= nro; i++)
        res[i] += x[i] * y;
    return res;
}

template<int nro, int nco1, int nco2> mat<nro, nco2> operator*(const mat<nro, nco1>& x, const mat<nco1, nco2>& y){
    mat<nro, nco2> res;
    for(int i = 1; i <= nro; i++)
        for(int j = 1; j <= nco2; j++)
            res[i][j] = x[i] * y.col(j);
    return res;
}

template<int nro, int nco> mat<nro, nco> operator*(const mat<nro, nco>& x, const double& y){
    mat<nro, nco> res;
    for(int i = 1; i <= nro; i++)
        for(int j = 1; j <= nco; j++)
            res[i][j] = x[i][j] * y;
    return res;
}

template<int nro, int nco> mat<nro, nco> operator/(const mat<nro, nco>& x, const double& y){
    mat<nro, nco> res;
    for(int i = 1; i <= nro; i++)
        for(int j = 1; j <= nco; j++)
            res[i][j] = x[i][j] / y;
    return res;
}

template<int nro, int nco> mat<nro, nco> operator+(const mat<nro, nco>& x, const mat<nro, nco>& y){
    mat<nro, nco> res;
    for(int i = 1; i <= nro; i++)
        for(int j = 1; j <= nco; j++)
            res[i][j] = x[i][j] + y[i][j];
    return res;
}

template<int nro, int nco> mat<nro, nco> operator-(const mat<nro, nco>& x, const mat<nro, nco>& y){
    mat<nro, nco> res;
    for(int i = 1; i <= nro; i++)
        for(int j = 1; j <= nco; j++)
            res[i][j] = x[i][j] - y[i][j];
    return res;
}

template<int nro, int nco> std::ostream& operator<<(std::ostream out, const mat<nro, nco>& x){
    for(int i = 1; i <= nro; i++)
        out << x[i] << "\n";
    return out;
}

template<int n> struct dt{
    static double det(const mat<n, n>& x){
        double res = 0;
        for(int i = 1; i <= n; i++)
            res += x[1][i] * x.cofactor(1, i);
        return res;
    }
};

template<>struct dt<1>{
    static double det(const mat<1, 1>& x){
        return x[1][1];
    }
};