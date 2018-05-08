/*!
 * @file src/common/spline.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__SPLINE__HEADER_INCLUDED
#define eSCATTER__COMMON__SPLINE__HEADER_INCLUDED

#include <map>

#include "archive.hh"
#include "polynomial.hh"

/*!
 *
 */
class spline {
friend archive::ostream& operator<<(archive::ostream&, const spline&);
friend archive::istream& operator>>(archive::istream&, spline&);
public:
    using const_iterator = std::map<double,polynomial>::const_iterator;
    static spline linear(const std::map<double,double>& xy_map);
    static spline akima(const std::map<double,double>& xy_map);
    double operator()(double x) const;
    spline operator-() const;
    spline operator+(double y) const;
    spline operator-(double y) const;
    spline operator*(double y) const;
    spline operator/(double y) const;
    spline& operator+=(double y);
    spline& operator-=(double y);
    spline& operator*=(double y);
    spline& operator/=(double y);
    const_iterator cbegin() const;
    const_iterator cend() const;
    spline& differentiate();
    const std::pair<double,double>& domain() const;
    spline& integrate(double x);
private:
    spline() = default;
    std::map<double,polynomial> _poly_map;
    std::pair<double,double> _domain_pair;
};

spline operator+(double y, const spline&);
spline operator*(double y, const spline&);

archive::ostream& operator<<(archive::ostream&, const spline&);
archive::istream& operator>>(archive::istream&, spline&);

#endif
