#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include <cassert>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
// you need to add the following libraries to your project : gsl, gslcblas

class Vec3 {
private:
    float mVals[3];
public:
    Vec3() {}
    Vec3( float x , float y , float z ) {
       mVals[0] = x; mVals[1] = y; mVals[2] = z;
    }
    float & operator [] (unsigned int c) { return mVals[c]; }
    float operator [] (unsigned int c) const { return mVals[c]; }
    void operator = (Vec3 const & other) {
       mVals[0] = other[0] ; mVals[1] = other[1]; mVals[2] = other[2];
    }
    float squareLength() const {
       return mVals[0]*mVals[0] + mVals[1]*mVals[1] + mVals[2]*mVals[2];
    }
    float length() const { return sqrt( squareLength() ); }
    void normalize() { float L = length(); mVals[0] /= L; mVals[1] /= L; mVals[2] /= L; }
    static float dot( Vec3 const & a , Vec3 const & b ) {
       //Fonction à compléter
       float res = 0.0f; //Faire le calcul ici et assigner le resultat à res

        for(int i =0; i<3; i++)
        {
            res += a[i]*b[i];
        }
        return res; 
    }
    static Vec3 cross( Vec3 const & a , Vec3 const & b ) {
       //Fonction à compléter
        Vec3 res; //Faire le calcul ici et assigner le resultat à res
        res[0] = a[1] * b[2] - a[2] * b[1];
        res[1] = (a[0] * b[2] - a[2] * b[0]) * (-1);
        res[2] = a[0] * b[1] - a[1] * b[0]; 
        return res; //A remplacer par le résultat du produit vectoriel
    }
    void operator += (Vec3 const & other) {
        mVals[0] += other[0];
        mVals[1] += other[1];
        mVals[2] += other[2];
    }
    void operator -= (Vec3 const & other) {
        mVals[0] -= other[0];
        mVals[1] -= other[1];
        mVals[2] -= other[2];
    }
    void operator *= (float s) {
        mVals[0] *= s;
        mVals[1] *= s;
        mVals[2] *= s;
    }
    void operator /= (float s) {
        mVals[0] /= s;
        mVals[1] /= s;
        mVals[2] /= s;
    }

    bool operator ==(Vec3 const &other)
    {
        return (mVals[0] == other[0]) && (mVals[1] == other[1]) && (mVals[2] == other[2]);
    }
    static
    Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return  Vec3( ThetaPhiR[2] *(cos(ThetaPhiR[0]) * cos(ThetaPhiR[1])) , ThetaPhiR[2] *(sin(ThetaPhiR[0]) * cos(ThetaPhiR[1])) , ThetaPhiR[2] * sin(ThetaPhiR[1]) );
    }

    static
    Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
        return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
    }

    static
    Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
        float R = xyz.length();
        float phi = asin( xyz[2] / R );
        float theta = atan2( xyz[1] , xyz[0] );
        return Vec3( theta , phi , R );
}
};

static inline Vec3 operator + (Vec3 const & a , Vec3 const & b) {
   return Vec3(a[0]+b[0] , a[1]+b[1] , a[2]+b[2]);
}
static inline Vec3 operator - (Vec3 const & a , Vec3 const & b) {
   return Vec3(a[0]-b[0] , a[1]-b[1] , a[2]-b[2]);
}
static inline Vec3 operator * (float a , Vec3 const & b) {
   return Vec3(a*b[0] , a*b[1] , a*b[2]);
}
static inline Vec3 operator / (Vec3 const &  a , float b) {
   return Vec3(a[0]/b , a[1]/b , a[2]/b);
}
static inline std::ostream & operator << (std::ostream & s , Vec3 const & p) {
    s << p[0] << " " << p[1] << " " << p[2];
    return s;
}
static inline std::istream & operator >> (std::istream & s , Vec3 & p) {
    s >> p[0] >> p[1] >> p[2];
    return s;
}




class Mat3
{
public:
    ////////////         CONSTRUCTORS          //////////////
    Mat3()
    {
        vals[0] = 0;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = 0;
        vals[4] = 0;
        vals[5] = 0;
        vals[6] = 0;
        vals[7] = 0;
        vals[8] = 0;
    }
    Mat3( float v1 , float v2 , float v3 , float v4 , float v5 , float v6 , float v7 , float v8 , float v9)
    {
        vals[0] = v1;
        vals[1] = v2;
        vals[2] = v3;
        vals[3] = v4;
        vals[4] = v5;
        vals[5] = v6;
        vals[6] = v7;
        vals[7] = v8;
        vals[8] = v9;
    }
    Mat3( const Mat3 & m )
    {
        for(int i = 0 ; i < 3 ; ++i )
            for(int j = 0 ; j < 3 ; ++j )
                (*this)(i,j) = m(i,j);
    }


    Vec3 operator * (const Vec3 & p) // computes m.p
    {
        //Fonction à completer
        //Pour acceder a un element de la matrice (*this)(i,j) et du point p[i]
        return Vec3(); //A remplacer par le résultat de la multiplication
    }

    Mat3 operator * (const Mat3 & m2)
    {
        //Fonction à completer
        //Pour acceder a un element de la premiere matrice (*this)(i,j) et de la deuxième m2(k,l)

        return Mat3(); //A remplacer par le résultat de la multiplication
    }

    bool isnan() const {
        return std::isnan(vals[0]) || std::isnan(vals[1]) || std::isnan(vals[2])
                 || std::isnan(vals[3]) || std::isnan(vals[4]) || std::isnan(vals[5])
                 || std::isnan(vals[6]) || std::isnan(vals[7]) || std::isnan(vals[8]);
    }

    void operator = (const Mat3 & m)
    {
        for(int i = 0 ; i < 3 ; ++i )
            for(int j = 0 ; j < 3 ; ++j )
                (*this)(i,j) = m(i,j);
    }

    void operator += (const Mat3 & m)
    {
        for(int i = 0 ; i < 3 ; ++i )
            for(int j = 0 ; j < 3 ; ++j )
                (*this)(i,j) += m(i,j);
    }
    void operator -= (const Mat3 & m)
    {
        for(int i = 0 ; i < 3 ; ++i )
            for(int j = 0 ; j < 3 ; ++j )
                (*this)(i,j) -= m(i,j);
    }
    void operator /= (double s)
    {
        for( unsigned int c = 0 ; c < 9 ; ++c )
            vals[c] /= s;
    }

    Mat3 operator - (const Mat3 & m2)
    {
        return Mat3( (*this)(0,0)-m2(0,0) , (*this)(0,1)-m2(0,1) , (*this)(0,2)-m2(0,2) , (*this)(1,0)-m2(1,0) , (*this)(1,1)-m2(1,1) , (*this)(1,2)-m2(1,2) , (*this)(2,0)-m2(2,0) , (*this)(2,1)-m2(2,1) , (*this)(2,2)-m2(2,2) );
    }
    Mat3 operator + (const Mat3 & m2)
    {
        return Mat3( (*this)(0,0)+m2(0,0) , (*this)(0,1)+m2(0,1) , (*this)(0,2)+m2(0,2) , (*this)(1,0)+m2(1,0) , (*this)(1,1)+m2(1,1) , (*this)(1,2)+m2(1,2) , (*this)(2,0)+m2(2,0) , (*this)(2,1)+m2(2,1) , (*this)(2,2)+m2(2,2) );
    }

    Mat3 operator / (float s)
    {
        return Mat3( (*this)(0,0)/s , (*this)(0,1)/s , (*this)(0,2)/s , (*this)(1,0)/s , (*this)(1,1)/s , (*this)(1,2)/s , (*this)(2,0)/s , (*this)(2,1)/s , (*this)(2,2)/s );
    }
    Mat3 operator * (float s)
    {
        return Mat3( (*this)(0,0)*s , (*this)(0,1)*s , (*this)(0,2)*s , (*this)(1,0)*s , (*this)(1,1)*s , (*this)(1,2)*s , (*this)(2,0)*s , (*this)(2,1)*s , (*this)(2,2)*s );
    }

    ////////        ACCESS TO COORDINATES      /////////
    float operator () (unsigned int i , unsigned int j) const
    { return vals[3*i + j]; }
    float & operator () (unsigned int i , unsigned int j)
    { return vals[3*i + j]; }

    ////////        BASICS       /////////
    inline float sqrnorm()
    {
        return vals[0]*vals[0] + vals[1]*vals[1] + vals[2]*vals[2]
                + vals[3]*vals[3] + vals[4]*vals[4] + vals[5]*vals[5]
                + vals[6]*vals[6] +  vals[7]*vals[7] + vals[8]*vals[8];
    }

    inline float norm()
    { return sqrt( sqrnorm() ); }

    inline float determinant() const
    {
        return vals[0] * ( vals[4] * vals[8] - vals[7] * vals[5] )
                - vals[1] * ( vals[3] * vals[8] - vals[6] * vals[5] )
                + vals[2] * ( vals[3] * vals[7] - vals[6] * vals[4] );
    }

    static
    Mat3 inverse( Mat3 const & m , double defaultValueForInverseSingularValue = 0.0 )
    {
        float det = m.determinant();
        if( fabs(det) != 0.0 )
        {
            return Mat3( m(1,1)*m(2,2) - m(2,1)*m(1,2) , m(0,2)*m(2,1) - m(0,1)*m(2,2) , m(0,1)*m(1,2) - m(0,2)*m(1,1) ,
                             m(1,2)*m(2,0) - m(1,0)*m(2,2) , m(0,0)*m(2,2) - m(0,2)*m(2,0) , m(0,2)*m(1,0) - m(0,0)*m(1,2) ,
                             m(1,0)*m(2,1) - m(1,1)*m(2,0) , m(0,1)*m(2,0) - m(0,0)*m(2,1) , m(0,0)*m(1,1) - m(0,1)*m(1,0) ) / det ;
        }

        // otherwise:
        Mat3 U ; float sx ; float sy ; float sz ; Mat3 Vt;
        m.SVD(U,sx,sy,sz,Vt);
        float sxInv = sx == 0.0 ? 1.0 / sx : defaultValueForInverseSingularValue;
        float syInv = sy == 0.0 ? 1.0 / sy : defaultValueForInverseSingularValue;
        float szInv = sz == 0.0 ? 1.0 / sz : defaultValueForInverseSingularValue;
        return Vt.getTranspose() * Mat3::diag(sxInv , syInv , szInv) * U.getTranspose();
    }

    void SVD( Mat3 & U , float & sx , float & sy , float & sz , Mat3 & Vt ) const
    {
        gsl_matrix * u = gsl_matrix_alloc(3,3);
        for(unsigned int i = 0 ; i < 3; ++i)
            for(unsigned int j = 0 ; j < 3; ++j)
                gsl_matrix_set( u , i , j , (*this)(i,j) );

        gsl_matrix * v = gsl_matrix_alloc(3,3);
        gsl_vector * s = gsl_vector_alloc(3);
        gsl_vector * work = gsl_vector_alloc(3);

        gsl_linalg_SV_decomp (u,
                              v,
                              s,
                              work);

        sx = s->data[0];
        sy = s->data[1];
        sz = s->data[2];
        for(unsigned int i = 0 ; i < 3; ++i)
        {
            for(unsigned int j = 0 ; j < 3; ++j)
            {
                U(i,j) = gsl_matrix_get( u , i , j );
                Vt(i,j) = gsl_matrix_get( v , j , i );
            }
        }
        assert( sx >= sy );
        assert( sy >= sz );

        gsl_matrix_free(u);
        gsl_matrix_free(v);
        gsl_vector_free(s);
        gsl_vector_free(work);

        // a transformation float is given as R.B.S.Bt, R = rotation , B = local basis (rotation matrix), S = scales in the basis B
        // it can be obtained from the svd decomposition of float = U Sigma Vt :
        // B = V
        // S = Sigma
        // R = U.Vt
    }


    inline float trace() const
    { return vals[0] + vals[4] + vals[8]; }

    ////////        TRANSPOSE       /////////
    inline
    void transpose()
    {
        float xy = vals[1] , xz = vals[2] , yz = vals[5];
        vals[1] = vals[3];
        vals[3] = xy;
        vals[2] = vals[6];
        vals[6] = xz;
        vals[5] = vals[7];
        vals[7] = yz;
    }
    Mat3 getTranspose() const
    {
        return Mat3(vals[0],vals[3],vals[6],vals[1],vals[4],vals[7],vals[2],vals[5],vals[8]);
    }

    // ---------- ROTATION <-> AXIS/ANGLE ---------- //
    template< class point_t >
    void getAxisAndAngleFromRotationMatrix( point_t & axis , float & angle )
    {
        angle = acos( (trace() - 1.f) / 2.f );
        axis[0] = vals[7] - vals[5];
        axis[1] = vals[2] - vals[6];
        axis[2] = vals[3] - vals[1];
        axis.normalize();
    }

    template< class point_t >
    inline static
    Mat3 getRotationMatrixFromAxisAndAngle( const point_t & axis , float angle )
    {
        Mat3 w = vectorial(axis);
        return Identity() +  w * std::sin(angle) +  w * w * ((1.0) - std::cos(angle));
    }

    // ---------- STATIC STANDARD MATRICES ---------- //
    inline static Mat3 Identity()
    {  return Mat3(1,0,0  ,  0,1,0  ,  0,0,1);  }

    inline static Mat3 Zero()
    {  return Mat3(0,0,0  ,  0,0,0  ,  0,0,0);  }

    template< typename T2 >
    inline static Mat3 diag( T2 x , T2 y ,T2 z )
    {  return Mat3(x,0,0  ,  0,y,0  ,  0,0,z);  }


    template< class point_t >
    inline static Mat3 getFromCols(const point_t & c1 , const point_t & c2 , const point_t & c3)
    {
        // 0 1 2
        // 3 4 5
        // 6 7 8
        return Mat3( c1[0] , c2[0] , c3[0] ,
                         c1[1] , c2[1] , c3[1] ,
                         c1[2] , c2[2] , c3[2] );
    }
    template< class point_t >
    inline static Mat3 getFromRows(const point_t & r1 , const point_t & r2 , const point_t & r3)
    {
        // 0 1 2
        // 3 4 5
        // 6 7 8
        return Mat3( r1[0] , r1[1] , r1[2] ,
                         r2[0] , r2[1] , r2[2] ,
                         r3[0] , r3[1] , r3[2] );
    }

    Mat3 operator - () const
    {
        return Mat3( - vals[0],- vals[1],- vals[2],- vals[3],- vals[4],- vals[5],- vals[6],- vals[7],- vals[8] );
    }


private:
    float vals[9];
    // will be noted as :
    // 0 1 2
    // 3 4 5
    // 6 7 8
};



inline static
Mat3 operator * (float s , const Mat3 & m)
{
    return Mat3( m(0,0)*s , m(0,1)*s , m(0,2)*s , m(1,0)*s , m(1,1)*s , m(1,2)*s , m(2,0)*s , m(2,1)*s , m(2,2)*s );
}


inline static std::ostream & operator << (std::ostream & s , Mat3 const & m)
{
    s << m(0,0) << " \t" << m(0,1) << " \t" << m(0,2) << std::endl << m(1,0) << " \t" << m(1,1) << " \t" << m(1,2) << std::endl << m(2,0) << " \t" << m(2,1) << " \t" << m(2,2) << std::endl;
    return s;
}


#endif
