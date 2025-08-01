// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matrix34.h"
#include "vector2.h"
#include "random.h"


Vector3 Matrix34::rotationAxis() const
{
    real x = value(2,1) - value(1,2);
    real y = value(0,2) - value(2,0);
    real z = value(1,0) - value(0,1);
    real n = std::sqrt(x*x + y*y + z*z);
    if ( n > REAL_EPSILON )
        return Vector3(x/n, y/n, z/n);
    // if 'n==0', the matrix is diagonal, and the rotation angle is 0 or M_PI
    if ( value(2, 2) > 0 )
        return Vector3(0, 0, 1);
    else if ( value(1, 1) > 0 )
        return Vector3(0, 1, 0);
    else
        return Vector3(1, 0, 0);
}

real Matrix34::rotationAngle() const
{
    // trace() = 1 + 2 * cos(angle)
    return std::acos(0.5*trace()-0.5);
}


void Matrix34::getEulerAngles(real& a, real& b, real& c) const
{
    real cb = std::sqrt(square(value(0,0)) + square(value(1,0)));
    
    b = std::atan2(-value(2,0), cb);
    
    if ( cb != 0 ) {
        a = std::atan2(value(1,0), value(0,0));
        c = std::atan2(value(2,1), value(2,2));
    }
    else {
        a = 0;
        c = std::atan2(-value(1,1), value(2,1));
    }
}




/// returns a rotation of angle PI around axis Z
Matrix34 Matrix34::rotation180()
{
    return Matrix34(-1, 0, 0, 0, -1, 0, 0, 0, 1);
}


Matrix34 Matrix34::flipX()
{
    return Matrix34(-1, 0, 0, 0, 1, 0, 0, 0, 1);
}


Matrix34 Matrix34::align111()
{
    real X = 1.0/M_SQRT3;
    real Y = -sqrt((2+M_SQRT3)/6.0);
    real Z = (3-M_SQRT3)/6.0;
    return Matrix34(X,Y,Z,X,Z,Y,X,X,X);
}


Matrix34 Matrix34::rotationAroundX(const real angle)
{
    real c = std::cos(angle);
    real s = std::sin(angle);
    return Matrix34(1, 0, 0, 0, c, s, 0, -s, c);
}

Matrix34 Matrix34::rotationAroundY(const real angle)
{
    real c = std::cos(angle);
    real s = std::sin(angle);
    return Matrix34(c, 0, -s, 0, 1, 0, s, 0, c);
}

Matrix34 Matrix34::rotationAroundZ(const real angle)
{
    real c = std::cos(angle);
    real s = std::sin(angle);
    return Matrix34(c, s, 0, -s, c, 0, 0, 0, 1);
}


Matrix34 Matrix34::rotationAroundPrincipalAxis(index_t i, const real angle)
{
    real c = std::cos(angle);
    real s = std::sin(angle);
    
    i %= 3;
    index_t j = (i+1)%3;
    index_t k = (i+2)%3;
    
    Matrix34 res(0, 1);
    res(j,j) = c;
    res(k,j) = s;
    res(j,k) = -s;
    res(k,k) = c;
    return res;
}


Matrix34 Matrix34::rotationFromAngles(const real a[3])
{
    real ca = std::cos(a[0]), sa = std::sin(a[0]);
    real cb = std::cos(a[1]), sb = std::sin(a[1]);
    real cc = std::cos(a[2]), sc = std::sin(a[2]);
    
    Matrix34 res;

    res(0,0) =  ca*cb;
    res(1,0) =  sa*cb;
    res(2,0) = -sb;
    
    res(0,1) =  ca*sb*sc - sa*cc;
    res(1,1) =  sa*sb*sc + ca*cc;
    res(2,1) =  cb*sc;
    
    res(0,2) =  ca*sb*cc + sa*sc;
    res(1,2) =  sa*sb*cc - ca*sc;
    res(2,2) =  cb*cc;

    return res;
}


Matrix34 Matrix34::rotationAroundAxisEuler(const real a[3])
{
    real ca = std::cos(a[0]), sa = std::sin(a[0]), ca1 = 1 - ca;
    real cb = std::cos(a[1]), sb = std::sin(a[1]);
    real cc = std::cos(a[2]), sc = std::sin(a[2]);
    
    real sacc      = sa * cc,         sasc    = sa * sc;
    real saccsb    = sacc * sb,       sacccb  = sacc * cb;
    real ccccca1   = cc * cc * ca1,   ccscca1 = cc * sc * ca1;
    real cbccccca1 = cb * ccccca1;
    
    Matrix34 res;
    res(0,0) = cb * cbccccca1 + ca;
    res(0,1) = sb * cbccccca1 - sasc;
    res(0,2) = cb * ccscca1   + saccsb;
    
    res(1,0) = sb * cbccccca1 + sasc;
    res(1,1) = ca - cb * cbccccca1 + ccccca1;
    res(1,2) = sb * ccscca1   - sacccb;
    
    res(2,0) = cb * ccscca1 - saccsb;
    res(2,1) = sb * ccscca1 + sacccb;
    res(2,2) = 1 - ccccca1;

    return res;
}


Matrix34 Matrix34::randomRotation()
{
    //James Arvo, Fast random rotation matrices. in Graphics Gems 3.
    real u1 = M_PI * RNG.sreal();
    real u2 = M_PI * RNG.sreal();
    real u3 = RNG.preal();
    real uu = std::sqrt(u3);
    Vector3 V( uu*std::cos(u2), uu*std::sin(u2), std::sqrt(1-u3) );
    return householder(V) * rotationAroundZ(std::cos(u1), std::sin(u1));
}


Matrix34 Matrix34::randomRotation(real angle)
{
    return rotationAroundAxis(Vector3::randU(), std::cos(angle), std::sin(angle));
}


Matrix34 Matrix34::rotationToVector(const Vector3& vec)
{
    Matrix34 res;
    Vector3 X, Y, Z = normalize(vec);
    Z.orthonormal(X, Y);
    res.setColumns(Z, X, Y);
    return res;
}


/**
 Set a matrix that transform (1,0,0) into X, and (0,0,1) into Z
 */
Matrix34 Matrix34::rotationToVectors(Vector3 X, Vector3 Z)
{
    X.normalize();
    Z = normalize(Z - dot(Z, X) * X);
    Vector3 Y = cross(Z, X);

    Matrix34 res;
    res.setColumns(X, Y, Z);
    return res;
}


Matrix34 Matrix34::randomRotationToVector(const Vector3& vec)
{
    Matrix34 res;
    Vector3 X, Y, Z = normalize(vec);
    Z.orthonormal(X, Y);
#if ( 0 )
    real a = M_PI * RNG.sreal();
    real c = std::cos(a), s = std::sin(a);
    res.setColumns(Z, X*c+Y*s, Y*c-X*s);
#else
    real C, S;
    RNG.urand2(C, S);
    res.setColumns(Z, X*C+Y*S, Y*C-X*S);
#endif
    return res;
}

