/*!
 * @file src/common/polynomial.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__POLYNOMIAL__HEADER_INCLUDED
#define eSCATTER__COMMON__POLYNOMIAL__HEADER_INCLUDED

#include <initializer_list>
#include <vector>

#include "archive.hh"

/*!
 *
 */
class polynomial {
public:
    polynomial(double x0 = 0);
    polynomial(const std::initializer_list<double>& ai, double x0 = 0);
    polynomial operator-() const;
    polynomial operator+(double y) const;
    polynomial operator-(double y) const;
    polynomial operator*(double y) const;
    polynomial operator/(double y) const;
    polynomial& operator+=(double y);
    polynomial& operator-=(double y);
    polynomial& operator*=(double y);
    polynomial& operator/=(double y);
    double operator()(double x) const;
    double operator[](int i) const;
    polynomial& differentiate();
    polynomial& integrate(double a0 = 0);
    int order() const;
    double origin() const;
    polynomial& set_coeff(int i, double ai);
    polynomial& transform(double x0);
private:
    std::vector<double> _coeff_vec;
    double _origin;
};

polynomial operator+(double y, const polynomial&);
polynomial operator*(double y, const polynomial&);

archive::ostream& operator<<(archive::ostream&, const polynomial&);
archive::istream& operator>>(archive::istream&, polynomial&);

#endif
