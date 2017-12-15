/**
 * @file src/cdsem/point3.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__POINT3__HEADER_INCLUDED
#define eSCATTER__CDSEM__POINT3__HEADER_INCLUDED

namespace legacy_thomas {

struct point3 {
    point3();
    point3(double _x, double _y, double _z);
    point3 operator+(const point3&) const;
    point3 operator-(const point3&) const;
    point3& operator+=(const point3&);
    point3& operator-=(const point3&);
    point3 operator*(double) const;
    point3 operator/(double) const;
    point3& operator*=(double);
    point3& operator/=(double);
    double norm() const;
    double x, y, z;
};

double dot_product(const point3&, const point3&);
point3 cross_product(const point3&, const point3&);

point3 rotate(const point3&, const point3& theta);

}

#endif
