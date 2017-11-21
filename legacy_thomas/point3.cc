/**
 * @file src/cdsem/point3.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include "point3.hh"
#include <cmath>

namespace legacy_thomas {

point3::point3() : point3(0, 0, 0) {
};

point3::point3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {
}

point3 point3::operator+(const point3& p) const {
    return point3(x+p.x, y+p.y, z+p.z);
}

point3 point3::operator-(const point3& p) const {
    return point3(x-p.x, y-p.y, z-p.z);
}

point3& point3::operator+=(const point3& p) {
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
}

point3& point3::operator-=(const point3& p) {
    x -= p.x;
    y -= p.y;
    z -= p.z;
    return *this;
}

point3 point3::operator*(double s) const {
    return point3(x*s, y*s, z*s);
}

point3 point3::operator/(double s) const {
    return point3(x/s, y/s, z/s);
}

point3& point3::operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

point3& point3::operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

double point3::norm() const {
    return std::sqrt(x*x+y*y+z*z);
}

double dot_product(const point3& p1, const point3& p2) {
    return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
}

point3 cross_product(const point3& p1, const point3& p2) {
    return point3(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x);
}

point3 rotate(const point3& p, const point3& t) {
    const double theta = t.norm();
    if(theta == 0)
        return p;
    const point3 rot_axis = t/theta;
    const point3 unit_u = rot_axis*dot_product(p, rot_axis);
    const point3 unit_v = p-unit_u;
    const point3 unit_w = cross_product(rot_axis, unit_v);
    return unit_u+unit_v*std::cos(theta)+unit_w*std::sin(theta);
}

}