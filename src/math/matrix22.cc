// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matrix22.h"
#include "random.h"


real Matrix22::rotationAngle() const
{
    return std::atan2(val[1], val[0]);
}


/// returns a rotation of angle PI around axis Z
Matrix22 Matrix22::rotation180()
{
    return Matrix22(-1.0, 0.0, -1.0, 0.0);
}

Matrix22 Matrix22::flipX()
{
    return Matrix22(-1, 0, 0, 1);
}


Matrix22 Matrix22::align111()
{
    return Matrix22(M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, -M_SQRT1_2);
}

Matrix22 Matrix22::randomRotation()
{
    real a = M_PI*RNG.sreal();
    real c = std::cos(a);
    real s = std::sin(a);
    return Matrix22(c, s, -s, c);
}


Matrix22 Matrix22::randomRotation(real angle)
{
    real c = std::cos(angle);
    real s = std::sin(angle) * RNG.sflip();
    return Matrix22(c, s, -s, c);
}


Matrix22 Matrix22::rotationToVector(const Vector2& vec)
{
    real n = vec.norm();
    real c = vec.XX / n;
    real s = vec.YY / n;
    return Matrix22(c, s, -s, c);
}

