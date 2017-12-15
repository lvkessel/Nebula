/**
 * @file src/cdsem/triangle.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__TRIANGLE__HEADER_INCLUDED
#define eSCATTER__CDSEM__TRIANGLE__HEADER_INCLUDED

#include "point3.hh"

namespace legacy_thomas {

struct triangle {
    enum id_enum : int {
        NOP = -128,
        TERMINATOR = -127,
        DETECTOR = -126,
        DETECTOR_LT50 = -125,
        DETECTOR_GE50 = -124,
        VACUUM = -123,
        MIRROR = -122
    };

    triangle(const point3& _A, const point3& _B, const point3& _C, int _in, int _out);

    point3 edge1() const;
    point3 edge2() const;
    point3 normal() const;

    point3 A, B, C;
    int in, out;
};

triangle rotate(const triangle&, const point3& theta);

}

#endif
