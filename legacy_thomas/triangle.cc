/**
 * @file src/cdsem/triangle.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include "triangle.hh"

namespace legacy_thomas {

triangle::triangle(const point3& _A, const point3& _B, const point3& _C, int _in, int _out)
: A(_A), B(_B), C(_C), in(_in), out(_out) {
}

point3 triangle::edge1() const {
    return B-A;
}

point3 triangle::edge2() const {
    return C-A;
}

point3 triangle::normal() const {
    return cross_product(edge1(), edge2());
}

triangle rotate(const triangle& t, const point3& theta) {
    return triangle(rotate(t.A, theta), rotate(t.B, theta), rotate(t.C, theta), t.in, t.out);
}

}