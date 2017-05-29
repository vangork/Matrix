#ifndef MATRIX_H_LY
#define MATRIX_H_LY

#include <iostream>
#include <cstring>
#include <cassert>
#include <cmath>
#include <vector>
#include "clapack.h"

namespace Matrix_ly
{

    enum __stype {FLOAT, DOUBLE, COMPLEX, DOUBLECOMPLEX};

    template<class T> 
        class storage_type;

    template<>
        class storage_type<double> 
        {
            public:
                static int get_type() {return DOUBLE;}
        };

    template<>
        class storage_type<float>
        {
            public:
                static int get_type() {return FLOAT;}
        };

    template<>
        class storage_type<complex>
        {
            public:
                static int get_type() {return COMPLEX;}
        };

    template<>
        class storage_type<doublecomplex>
        {
            public:
                static int get_type() {return DOUBLECOMPLEX;}
        };


    //template<class T>
    //struct mat_ref
    //{
    //    T * _data;
    //    int _ref_counter;
    //};

    template<typename _Tp>
        bool operator<(const std::complex<_Tp> &__x, const std::complex<_Tp> &__y)
        {
            _Tp _a = __x.real() * __x.real() + __x.imag() * __x.imag();
            _Tp _b = __y.real() * __y.real() + __y.imag() * __y.imag();
            return _a < _b;
        }

    template<typename _Tp>
        bool operator>(const std::complex<_Tp> &__x, const std::complex<_Tp> &__y)
        {
            _Tp _a = __x.real() * __x.real() + __x.imag() * __x.imag();
            _Tp _b = __y.real() * __y.real() + __y.imag() * __y.imag();
            return _a > _b;
        }

    template<typename _Tp>
        bool operator<=(const std::complex<_Tp> &__x, const std::complex<_Tp> &__y)
        {
            _Tp _a = __x.real() * __x.real() + __x.imag() * __x.imag();
            _Tp _b = __y.real() * __y.real() + __y.imag() * __y.imag();
            return _a <= _b;
        }

    template<typename _Tp>
        bool operator>=(const std::complex<_Tp> &__x, const std::complex<_Tp> &__y)
        {
            _Tp _a = __x.real() * __x.real() + __x.imag() * __x.imag();
            _Tp _b = __y.real() * __y.real() + __y.imag() * __y.imag();
            return _a >= _b;
        }

    template<typename _Tp>
        std::ostream & operator<<(std::ostream &__os, const std::complex<_Tp> &__z)
        {
            //if(abs(m.r) > ZERO)
            //    os << m.r;
            //if(abs(m.r) > ZERO && abs(m.i) > ZERO)
            //    os << " + ";
            //if(abs(m.i) > ZERO)
            //    os << m.i << "i";
            __os << __z.real() << " + " << __z.imag() << "i";
            return __os;
        }

    template<class _Tp>
        class vec
        {
            public:

                vec(const int __l = 0)
                    : _M_len(__l)
                {
                    _M_data = new _Tp[__l];
                    memset(_M_data, 0, sizeof(_Tp) * __l);
                }

                vec(const _Tp* __d, const int __l)
                    :_M_len(__l)
                {
                    _M_data = new _Tp[__l];
                    memcpy(_M_data, __d, sizeof(_Tp) * __l);
                }

                vec(const vec<_Tp>& __z)
                {
                    _M_len = __z.len();
                    _M_data = new _Tp[_M_len];
                    memcpy(_M_data, __z.data(), sizeof(_Tp) * _M_len);
                }

                ~vec()
                {
                    delete[] _M_data;
                }

            public:
                void set(const _Tp* __d, const int __l)
                {
                    delete[] _M_data; 
                    _M_len = __l;
                    _M_data = new _Tp[__l];
                    memcpy(_M_data, __d, sizeof(_Tp) * __l);
                }

                void set(const int __l = 0)
                {
                    delete[] _M_data; 
                    _M_len = __l;
                    _M_data = new _Tp[__l];
                    memset(_M_data, 0, sizeof(_Tp) * __l);
                }

                int len() const 
                { return _M_len; }

                _Tp* data() const 
                { return _M_data; }

                bool is_inbound(const int __t) const
                { return __t >= 0 && __t < _M_len; }

                void zeros()
                { memset(_M_data, 0, sizeof(_Tp) * _M_len); }

                _Tp& operator[](const int __t)
                {
                    assert(is_inbound(__t));
                    return _M_data[__t];
                }

                const _Tp& operator[](const int __t) const
                {
                    assert(is_inbound(__t));
                    return _M_data[__t];
                }

            public:
                vec<_Tp> operator+(const vec<_Tp>&);

                vec<_Tp> operator-(const vec<_Tp>&);

                template<class _Up>
                    vec<_Tp>& operator=(const vec<_Up>&);

                vec<_Tp>& operator=(const vec<_Tp>&);

                template<class _Up>
                    vec<_Tp> operator^(const _Up);

                vec<doublecomplex> operator^(const doublecomplex&);

                vec<_Tp> exp();

                _Tp max();

            private:
                int _M_len;
                _Tp* _M_data;
        };

    template<class _Tp>
        vec<_Tp> vec<_Tp>::operator+(const vec<_Tp>& __z)
        { 
            //daxpy to be modified
            if(_M_len != __z.len())
            {
                std::cout << "Vector dimensions must agree" << std::endl;
                assert(0);
            }
            vec<_Tp> __r(_M_len);
            for(int __i = 0; __i < _M_len ; ++__i)
            {
                __r[__i] = _M_data[__i] + __z[__i];
            }
            return __r;
        }

    template<class _Tp>
        vec<_Tp> vec<_Tp>::operator-(const vec<_Tp>& __z)
        { 
            //daxpy to be modified
            if(_M_len != __z.len())
            {
                std::cout << "Vector dimensions must agree" << std::endl;
                assert(0);
            }
            vec<_Tp> __r(_M_len);
            for(int __i = 0; __i < _M_len ; ++__i)
            {
                __r[__i] = _M_data[__i] - __z[__i];
            }
            return __r;
        }

    template<class _Tp>
        vec<_Tp>& vec<_Tp>::operator=(const vec<_Tp>& __z)
        {
            if(this == &__z)
                return *this;
            delete[] _M_data;
            _M_len = __z.len();
            _M_data = new _Tp[_M_len];
            memcpy(_M_data, __z.data(), sizeof(_Tp) * _M_len);
            return *this;
        }

    template<class _Tp>
        template<class _Up>
        vec<_Tp>& vec<_Tp>::operator=(const vec<_Up>& __z)
        {
            if(this == &__z)
                return *this;
            delete[] _M_data;
            _M_len = __z.len();
            _M_data = new _Tp[_M_len];
            for(int __i = 0; __i < _M_len ; ++__i)
                _M_data[__i] = __z[__i];
            return *this;
        }

    template<class _Tp>
        template<class _Up>
        vec<_Tp> vec<_Tp>::operator^(const _Up __t)
        {
            vec<_Tp> __r(_M_len);
            for(int __i = 0; __i < _M_len; ++__i)
                __r[__i] = //::exp(std::log(_M_data[__i]) * __t);
                    //std::exp(std::log(_M_data[__i]) * __t);
                    std::pow(_M_data[__i], __t);
            return __r;
        }

    template<class _Tp>
        vec<doublecomplex> vec<_Tp>::operator^(const doublecomplex& __t)
        {
            vec<doublecomplex> __r(_M_len);
            for(int __i = 0; __i < _M_len; ++__i)
                __r[__i] = //::exp(std::log(_M_data[__i]) * __t);
                    //std::exp(std::log(_M_data[__i]) * __t);
                    std::pow(_M_data[__i], __t);
            return __r;
        }

    template<class _Tp>
        vec<_Tp> vec<_Tp>::exp()
        {
            vec<_Tp> __r(_M_len);
            for(int __i = 0; __i < _M_len; ++__i)
                __r[__i] = std::exp(_M_data[__i]);
            return __r;
        }

    template<class _Tp>
        _Tp vec<_Tp>::max()
        {
            if(_M_len == 0)
            {
                std::cout << "Empty Vector have no max element" << std::endl;
                assert(0);
            }
            _Tp __r = *_M_data;
            for(int __i = 1; __i < _M_len; ++__i)
                if(__r < _M_data[__i])
                    __r = _M_data[__i];
            return __r;
        }

    template<class _Tp>
        std::ostream & operator<<(std::ostream &__os, const vec<_Tp> &__z)
        {
            int __l = __z.len();
            __os << "[" ;
            int __i;
            for(__i = 0; __i < __l - 1; ++__i)
                __os << __z[__i] << "    ";
            if(__l)
                __os << __z[__i];
            __os << "]";
            return __os;
        }

    template<class T>
        class mat
        {
            public:
                mat():_row(0),_col(0){_data = new T[0];}

                mat(int row, int col):_row(row),_col(col)
                {
                    assert(row > 0 && col > 0 && "Invalid row or col numbers.\n");
                    _data = new T[row * col];
                    memset(_data, 0, sizeof(T) * _row * _col);
                }

                mat(const T *data, int row, int col, bool is_col_major = true);

                mat(const std::vector<T> &v, int n, bool is_lower_matrix = true);

                mat(const mat<T> &m);

                void set(const T *data, int row, int col, bool is_col_major = true);

                void set(int row, int col);
                
                void set(const std::vector<T> &v, int n, bool is_lower_matrix = true);

                ~mat()
                {
                    delete[] _data;
                }

                int row() const {return _row;}

                int col() const {return _col;}

                T *data() const {return _data;}

                void zeros()
                {
                    memset(_data, 0, sizeof(T) * _row * _col);
                }

                mat<T> operator+(const mat<T> &m);

                mat<T> operator-(const mat<T> &m);

                mat<T> operator-();

                mat<T> operator*(const mat<T> &m);

                mat<T> operator*(const double &t);

                friend mat<T> operator*(const double &t, mat<T> &m)
                {
                    mat<T> ret = m * t;
                    return ret;
                }

                mat<doublecomplex> operator*(const doublecomplex &t) ;

                friend mat<doublecomplex> operator*(const doublecomplex &t, mat<T> &m)
                {
                    mat<T> ret = m * t;
                    return ret;
                }

                //mat<T> operator/(const mat<T> &m);
                mat<T> backslash(const mat<T> &m);

                mat<T>& operator=(const mat<T> &m);

                mat<doublecomplex> operator^(const double &d);

                T& operator()(int m,int n);

                const T& operator()(int m,int n) const;

                void eig(vec<doublecomplex> &eval, mat<doublecomplex> &evec);

                mat<T> inv();

                mat<doublecomplex> expm();

                //mat<doublecomplex> to_dcmat(); 

                void diag(const vec<T> &v);

                T sum();
                mat<T> sum(bool isRow);

                vec<T> to_vec();

                //template <class typedef> friend std::ostream &operator<<(
                //        std::ostream &os), const mat<typedef> &m);

            private:
                bool is_inbound(int m, int n) const
                {
                    return m >= 0 && m < _row && n >= 0 && n < _col;
                }
                static bool SmallEnough(double a, double tol = 1e-64) { return fabs(a) < tol; }

            private:
                T * _data;
                int _row;
                int _col;
        };

    template <class T>
        std::ostream &operator<<(std::ostream &os, const mat<T> &m);

    template<class T>
        mat<T>::mat(const T *data, int row, int col, bool is_col_major)
        {
            assert(row > 0 && col > 0 && "Invalid row or col numbers.\n");
            _row = row;
            _col = col;
            _data = new T[row * col];
            if(is_col_major)
                memcpy(_data ,data, sizeof(T) * _row * _col);
            else
                for(int i = 0; i < _col; i++)
                    for (int j = 0; j < _row; j++)
                        _data[i * _row + j] = data[j * _col + i];
        }

    template<class T>     
        mat<T>::mat(const std::vector<T> &v, int n, bool is_lower_matrix)
        {
            assert(v.size() == n * (n + 1) / 2 && "Illegal matrix deimention.\n");
            _row = n;
            _col = n;
            _data = new T[n * n];
            for(int i = 0; i < n; i++)
                for(int j = 0; j < n; j++)
                {
                    if(i > j)
                        _data[i * n + j] = _data[j * n + i];
                    else
                        _data[i * n + j] = v[i + j*(j+1)/2];
                }

        }

    template<class T>
        mat<T>::mat(const mat<T> &m)
        {
            _row = m.row();
            _col = m.col();
            _data = new T[_row * _col];
            memcpy(_data, m.data(), sizeof(T) * _row * _col);
        }

    template<class T>
        void mat<T>::set(const T *data, int row, int col, bool is_col_major)
        {
            assert(row > 0 && col > 0 && "Invalid row or col numbers.\n");
            delete[] _data;
            _row = row;
            _col = col;
            _data = new T[_row * _col];
            if(is_col_major)
                memcpy(_data ,data, sizeof(T) * _row * _col);
            else
                for(int i = 0; i < _col; i++)
                    for (int j = 0; j < _row; j++)
                        _data[i * _row + j] = data[j * _col + i];
        }

    template<class T>     
        void mat<T>::set(const std::vector<T> &v, int n, bool is_lower_matrix)
        {
            assert(v.size() == (size_t)(n * (n + 1) / 2) && "Illegal matrix deimention");
            delete[] _data;
            _row = n;
            _col = n;
            _data = new T[n * n];
            for(int i = 0; i < n; i++)
                for(int j = 0; j < n; j++)
                {
                    if(i > j)
                        _data[i * n + j] = _data[j * n + i];
                    else
                        _data[i * n + j] = v[i + j*(j+1)/2];
                }
        }

    template<class T>
        void mat<T>::set(int row, int col)
        {
            assert(row > 0 && col > 0 && "Invalid row or col numbers.\n");
            delete[] _data;
            _row = row;
            _col = col;
            _data = new T[_row * _col];
            memset(_data, 0, sizeof(T) * _row * _col);
        }

    template<class T>
        mat<T> mat<T>::operator+(const mat<T> &m)
        { 
            //daxpy to be modified
            assert(_row == m.row() && _col == m.col() && "Matrix dimensions must agree.\n");
            mat<T> r(_row, _col);
            for(int i = 0; i < _row * _col ; i++)
            {
                (r.data())[i] = _data[i] + (m.data())[i];
            }
            return r;
        }

    template<class T>
        mat<T> mat<T>::operator-(const mat<T> &m)
        {
            assert(_row == m.row() && _col == m.col() && "Matrix dimensions must agree.\n");
            mat<T> r(_row, _col);
            for(int i = 0; i < _row * _col ; i++)
            {
                (r.data())[i] = _data[i] - (m.data())[i];
            }
            return r;
        }

    //saxpy to be modified
    template<class T>
        mat<T> mat<T>::operator-()
        {
            mat<T> r(_row, _col);
            for(int i = 0; i < _row * _col ; ++i)
            {
                (r.data())[i] = -(data())[i];
            }
            return r;
        }

    template<class T>
        mat<T> mat<T>::operator*(const mat<T> &m)
        {
            assert(_col == m.row() && "Matrix dimensions must agree.\n"); 
            int M = _row;
            int K = _col;
            int N = m.col();
            char no_trans = 'N';
            mat<T> r(M, N);

            int t = storage_type<T>::get_type();
            switch(t)
            {
                case DOUBLE:
                    {
                        double alpha = 1.0;
                        double beta = 0.0;
                        dgemm_(&no_trans, &no_trans, &M, &N, &K,
                                &alpha, (double *)_data, &M,
                                (double *)(m.data()), &K,
                                &beta, (double *)(r.data()), &M);
                        break;
                    }
                    //case COMPLEX:
                    //    {
                    //        complex<float> alpha(1.0, 0.0), beta(0.0, 0.0);
                    //        cgemm_(&no_trans, &no_trans, &M, &N, &K,
                    //               &alpha, (complex *)_data, &M,
                    //               (complex *)m.data(), &K,
                    //               &beta, (complex *)c, &M);
                    //        break;
                    //    }
                case DOUBLECOMPLEX:
                    {
                        doublecomplex alpha(1.0, 0.0), beta(0.0, 0.0);
                        zgemm_(&no_trans, &no_trans, &M, &N, &K,
                                &alpha, (doublecomplex *)_data, &M,
                                (doublecomplex *)m.data(), &K,
                                &beta, (doublecomplex *)(r.data()), &M);
                        break;
                    }

            }
            return r;
        }

    template<class T>
        mat<T> mat<T>::operator*(const double &t)
        {
            mat<T> r(_row, _col);
            for(int i = 0; i < _col * _row; i++)
                (r.data())[i] = _data[i] * t;
            return r;
        }

    template<class T>
        mat<doublecomplex> mat<T>::operator*(const doublecomplex &t)
        {
            mat<doublecomplex> r(_row, _col);
            for(int i = 0; i < _col * _row; i++)
                (r.data())[i] = _data[i] * t;
            return r;
        }

    template<class T>
        mat<T> mat<T>::backslash(const mat<T> &m)
        {
            assert(_row == m.row() && "Matrix dimensions must agree.\n"); 
            int M = _row;
            int N = _col;
            int nrhs = m.col();
            T *a = new T[M * N];
            memcpy(a, _data, sizeof(T) * M * N);
            int ldb = M > N ? M : N;
            T *b = new T[ldb * nrhs];
            memcpy(b, m.data(), sizeof(T) * M * nrhs);

            
            int *jpvt = new int[N];//ipiv(xGESV) the same dimension
            memset(jpvt, 0, sizeof(int) * N);

            T *work ;
            int info = 0;

            int t = storage_type<T>::get_type();
            switch(t)
            {
                case DOUBLE:
                    {
                        if(M != N)
                        {
                            //int lwork = ;
                            //work = new T[lwork];
                            //char trans = 'N';
                            //dgels_(&trans, &M, &N, &nrhs, (double *)a, &M,
                            //        (double *)b, &ldb,
                            //        (double *)work, &lwork, &info);

                            //double rcond=DLAMCH('P')
                            //dgelss_(&M, &N, &nrhs, (double *)a, &M,
                            //        (double *)b, &ldb, s, &rcond, &rank,
                            //        (double *)work, &lwork, &info);

                            char cmach = 'P'; 
                            double rcond = dlamch_(&cmach);
                            int rank;
                            //std::cout << "rcond = " << rcond << std::endl;

                            int lwork = 2 * M + 2 * N + nrhs;//MAX( MN+3*N+1, 2*MN+NRHS )
                            work = new T[lwork];

                            dgelsy_(&M, &N, &nrhs, (double *)a, &M, 
                                    (double *)b, &ldb, jpvt, &rcond, &rank, 
                                    (double *)work, &lwork, &info);
                            //std::cout << "rank = " << rank << std::endl;
                            
                            
                            
                            delete[] work;
                        }
                        else
                        {
                            //func dgesv_ solve the equation Ax = B (x must be square)
                            //DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
                            dgesv_(&M, &nrhs, (double *)a, &M, jpvt, (double *)b, &ldb, &info);
                        }
                        break;
                    }
                case DOUBLECOMPLEX:
                    {
                        if(M != N)
                        {
                            char cmach = 'P'; 
                            double rcond = dlamch_(&cmach);
                            int rank;
                            //std::cout << "rcond = " << rcond << std::endl;

                            int lwork = M + 2 * M + nrhs; //MN + MAX( 2*MN, N+1, MN+NRHS )
                            work = new T[lwork];
                            double *rwork = new double[2 * N];

                            //ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, RWORK, INFO )
                            zgelsy_(&M, &N, &nrhs, (doublecomplex *)a, &M,
                                    (doublecomplex *)b, &ldb, jpvt, &rcond, &rank,
                                    (doublecomplex *)work, &lwork, rwork, &info);
                            //std::cout << "rank = " << rank << std::endl;
                            delete[] work;
                            delete[] rwork;
                        }
                        else
                        {
                            zgesv_(&M, &nrhs, (doublecomplex *)a, &M, jpvt, (doublecomplex *)b, &ldb, &info);
                        }
                        break;
                    }
            }
            assert((info >= 0) && "BackSlash error!");
            if(info > 0)
                std::cout << "Warning: Matrix is singular to working precision. \n";
            mat<T> ret(b, N, nrhs);
            delete[] a;
            delete[] b;
            delete[] jpvt;
            
            return ret;
        }

    template<class T>
        mat<T>& mat<T>::operator=(const mat<T> & m)
        {
            if(this == &m)
                return *this;
            delete[] _data;
            _row = m.row();
            _col = m.col();
            _data = new T[_row * _col];
            memcpy(_data, m.data(), sizeof(T) * _row * _col);
            return *this;
        }

    template<class T>
        mat<doublecomplex> mat<T>::operator^(const double &d)
        {
            assert(_row == _col && "Matrix must be square.\n");
            mat<doublecomplex> q;
            vec<doublecomplex> lamda;
            eig(lamda, q);
            //std::cout << "lamda = " << lamda << std::endl;
            //std::cout << "q = " << q << std::endl;

            lamda = lamda ^ d;
            mat<doublecomplex> mlamda;
            mlamda.diag(lamda);
            mat<doublecomplex> r;
            r = q * mlamda;
            r = r * q.inv();
            return r;
        }

    template<class T>
        T& mat<T>::operator()(int m,int n)
        {
            assert(is_inbound(m, n) && "Index exceeds matrix dimensions.\n");
            return _data[n * _row + m];
        }

    template<class T>
        const T& mat<T>::operator()(int m,int n) const
        {
            assert(is_inbound(m, n) && "Index exceeds matrix dimensions.\n");
            return _data[n * _row + m];
        }
 
    template<class T>
        void mat<T>::eig(vec<doublecomplex> &eval, mat<doublecomplex> &evec)
        {
            assert(_row == _col && "Matrix must be square.\n");

            //parameters for *geev
            char jobvl = 'N';
            char jobvr = 'V';
            int n = _row;
            T *vl = new T[1];
            int ldvl = 1;
            int ldvr = n;
            int lwork = n * 4;
            T *work = new T[lwork];
            int info = 0;
            T *a = new T[n * n];
            memcpy(a, _data, sizeof(T) * n * n);

            int t = storage_type<T>::get_type();
            switch(t)
            {
                case DOUBLE:
                    {
                        double *vr = new double[ldvr * n];
                        double *wr = new double[n];
                        double *wi = new double[n];
                        dgeev_(&jobvl, &jobvr, &n, (double *)a, &n, 
                                (double *)wr, (double *)wi,
                                (double *)vl, &ldvl,
                                (double *)vr, &ldvr, 
                                (double *)work, &lwork, &info);

                        assert(info == 0 && "Matrix is singular to working precision.");

                        eval.set(n);
                        evec.set(ldvr, n);
                        for(int i = 0; i < n;)
                        {
                            eval[i].real() = wr[i];
                            eval[i].imag() = wi[i];

                            if(SmallEnough(wi[i]))
                            {
                                for (int j = 0; j < ldvr; j++)
                                {
                                    evec(j,i).real() = vr[ldvr * i + j];
                                    evec(j,i).imag() = 0.0;
                                }
                                i++;
                            }
                            else
                            {
                                for (int j = 0; j < ldvr; j++)
                                {
                                    evec(j,i).real() = vr[ldvr * i + j];
                                    evec(j,i).imag() = vr[ldvr * (i + 1) + j];
                                    evec(j,i + 1).real() = vr[ldvr * i + j];
                                    evec(j,i + 1).imag() = - vr[ldvr * (i + 1) + j];
                                }
                                eval[i + 1].real() = wr[i + 1];
                                eval[i + 1].imag() = wi[i + 1];
                                i = i + 2;
                            }
                        }
                        
                        
                        delete[] wr;
                        delete[] wi;
                        delete[] vr;
                        break;
                    }
                case DOUBLECOMPLEX:
                    {
                        doublecomplex *vr = new doublecomplex[ldvr * n];
                        doublecomplex *w = new doublecomplex[n];
                        doublecomplex *rwork = new doublecomplex[2 * n];
                        zgeev_(&jobvl, &jobvr, &n,(doublecomplex *)a, &n, 
                                (doublecomplex *)w, (doublecomplex *)vl, &ldvl, 
                                (doublecomplex *)vr, &ldvr, (doublecomplex *)work,
                                &lwork, (double *)rwork, &info);
                        eval.set(w, n);
                        evec.set((doublecomplex *)vr, n, n);
                        delete[] vr;
                        delete[] w;
                        delete[] rwork;

                        break;
                    }
                default:
                    std::cout << "The type of input data have no defination" << std::endl;
            }
            delete[] a;
            delete[] vl;
            delete[] work;

        }


    template<class T>
        mat<T> mat<T>::inv()
        {
            assert(_row == _col && "Matrix must be square.\n");
            int n = _row;
            int lda = n;
            int *ipiv = new int[n];
            int lwork = n * n;
            T *work = new T[lwork];
            int info = 0;
            T *a = new T[n * n];
            memcpy(a, _data, sizeof(T) * n * n);

            int t = storage_type<T>::get_type();
            switch(t)
            {
                case DOUBLE:
                    {
                        dgetrf_(&n, &n, (double *)a, &lda, ipiv, &info);
                        assert(!info);
                        dgetri_(&n, (double *)a, &lda, ipiv, (double *)work, &lwork, &info); 
                        break;
                    }
                case DOUBLECOMPLEX:
                    {
                        zgetrf_(&n, &n, (doublecomplex *)a, &lda,
                                ipiv, &info);
                        assert(!info);
                        zgetri_(&n, (doublecomplex *)a, &lda, ipiv, 
                                (doublecomplex *)work, &lwork, &info);
                        break;
                    }
                default:
                    std::cout << "The type of input data have no defination" << std::endl;
            }
            assert(info == 0 && "Matrix is singular to working precision.");
            mat<T> m(a, n, n);
            delete[] ipiv;
            delete[] work;
            delete[] a;
            return m;
        }

    template<class T>
        mat<doublecomplex> mat<T>::expm()
        {
            assert(_row == _col && "Matrix must be square.\n");
            mat<doublecomplex> q;
            vec<doublecomplex> lamda;
            eig(lamda, q);
            lamda = lamda.exp();
            mat<doublecomplex> mlamda;
            mlamda.diag(lamda);
            mat<doublecomplex> r;
            r = q * mlamda;
            r = r * q.inv();
            return r;
        }

    //template<class T>
    //mat<doublecomplex> mat<T>::to_dcmat()
    //{
    //    int t = storage_type<T>::get_type();
    //    doublecomplex *data = new doublecomplex[_row * _col];
    //    switch(t)
    //    {
    //    case DOUBLE:
    //        {
    //            for(int n = 0; n < _row * _col; n++)
    //                {
    //                    data[n].r = _data[n];
    //                    data[n].i = 0.0;
    //                }
    //        }
    //    }
    //    mat<doublecomplex> m(data, _row, _col);
    //    delete[] data;
    //    return m;
    //}

    template<class T>
        void mat<T>::diag(const vec<T> &v)
        {
            delete[] _data;
            _row = v.len();
            _col = _row;
            _data = new T[_row * _col];
            memset(_data, 0, sizeof(T) * _row * _col);
            for(int i = 0; i < _row; i++)
            {
                _data[i * _col + i] = v[i];
            }
        }

    template<class T>
        T mat<T>::sum()
        {
            T s = 0;
            for(int i = 0; i < _row * _col; ++i)
            {
                s += _data[i];
            }
            return s;
        }

    template<class T>
        mat<T> mat<T>::sum(bool isRow)
        {
            if(isRow)
            {
                mat<T> s(1, _col);
                for(int i = 0; i < _col; i++)
                    for(int j = 0; j < _row; j++)
                        s(0, i) += (*this)(i, j);
                return s;
            }
            else
            {
                mat<T> s(_row, 1);
                for(int i = 0; i < _row; i++)
                    for(int j = 0; j < _col; j++)
                        s(i, 0) += (*this)(i, j);
                return s;
            }
        }

    template <class T>
        vec<T> mat<T>::to_vec()
        {
            if(_row == 1 || _col == 1)
            {
                int s = _row * _col;
                vec<T> r(s);
                for(int i = 0; i < s; ++i)
                {
                    r[i] = _data[i];
                }
                return r;
            }
            else
            {
                assert(0 && "This Matrix does not have this method\n");
            }
        }


    template <class T>
        std::ostream &operator<<(std::ostream &os, const mat<T> &m)
        {
            int row = m.row();
            int col = m.col();

            if(!row)
                os << "[]";
            else
            {
                os << "[" ;
                for(int i = 0; i < row; i++)
                {
                    if(i > 0)
                        os << " ";
                    for(int j = 0; j < col; j++)
                    { 
                        T d = m(i,j);
                        //if((abs(d)) < ZERO)
                        //    os << 0;
                        //else
                        os << d;
                        if(j != col - 1)
                            os << "    ";
                    }
                    if(i == row - 1)
                        os << "]";
                    else
                        os << ";" << std::endl;
                }
            }
            return os;
        }

    typedef vec<double> drvec;
    typedef vec<doublecomplex> dcvec;
    typedef mat<double> drmat;
    typedef mat<doublecomplex> dcmat;


    void dcmat2dr(const dcmat &__m, drmat &__n);
    void dcvec2dr(const dcvec &__m, drvec &__n);
}


#endif

