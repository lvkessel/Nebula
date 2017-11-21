/**
 * @file src/common/spline.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "spline.hh"

#include <algorithm>
#include <cmath>
#include <deque>

spline spline::linear(const std::map<double,double>& xy_map) {
    spline s;
    if(xy_map.empty()) {
        return spline();
    } else if(xy_map.size() == 1) {
        const double x = xy_map.cbegin()->first;
        const double y = xy_map.cbegin()->second;
        s._poly_map[0] = polynomial({y}, x);
        s._domain_pair = std::make_pair(x, x);
        return s;
    }
    std::deque<double> x, y;
    for(auto cit = xy_map.cbegin(); cit != xy_map.cend(); cit++) {
        x.push_back(cit->first);
        y.push_back(cit->second);
    }
    for(size_t i = 0; i < xy_map.size()-1; i++)
        s._poly_map[x[i]] = polynomial(x[i])
            .set_coeff(0, y[i])
            .set_coeff(1, (y[i+1]-y[i])/(x[i+1]-x[i]));
    s._domain_pair.first = xy_map.cbegin()->first;
    s._domain_pair.second = xy_map.crbegin()->first;
    return s;
}

spline spline::akima(const std::map<double,double>& xy_map) {
    spline s;
    if(xy_map.empty()) {
        return spline();
    } else if(xy_map.size() == 1) {
        const double x = xy_map.cbegin()->first;
        const double y = xy_map.cbegin()->second;
        s._poly_map[0] = polynomial({y}, x);
        s._domain_pair = std::make_pair(x, x);
        return s;
    } else if(xy_map.size() == 2) {
        return linear(xy_map);
    }
     std::deque<double> x, y, m;
     for(auto cit = xy_map.cbegin(); cit != xy_map.cend(); cit++) {
         x.push_back(cit->first);
         y.push_back(cit->second);
     }
     for(size_t i = 0; i < xy_map.size()-1; i++)
         m.push_back((y[i+1]-y[i])/(x[i+1]-x[i]));
     for(size_t i = 0; i < 2; i++) {
         x.push_front(2*x[0]-x[1]);
         y.push_front(y[0]-(2.0*m[0]-m[1])*(x[1]-x[0]));
         m.push_front((y[1]-y[0])/(x[1]-x[0]));
     }
     for(size_t i = xy_map.size()-1; i < xy_map.size()+2; i++) {
         x.push_back(-x[i]+x[i+1]+x[i+2]);
         y.push_back(y[i+2]+(2.0*m[i+1]-m[i])*(x[i+3]-x[i+2]));
         m.push_back((y[i+3]-y[i+2])/(x[i+3]-x[i+2]));
     }
     for(size_t i = 2; i < xy_map.size()+1; i++) {
         double t1 = m[i];
         if((m[i+1] != m[i]) || (m[i-1] != m[i-2]))
             t1 = (std::fabs(m[i+1]-m[i])*m[i-1]+std::fabs(m[i-1]-m[i-2])*m[i])/(std::fabs(m[i+1]-m[i])+std::fabs(m[i-1]-m[i-2]));
         double t2 = m[i+1];
         if((m[i+2] != m[i+1]) || (m[i] != m[i-1]))
             t2 = (std::fabs(m[i+2]-m[i+1])*m[i]+std::fabs(m[i]-m[i-1])*m[i+1])/(std::fabs(m[i+2]-m[i+1])+std::fabs(m[i]-m[i-1]));
         s._poly_map[x[i]] = polynomial(x[i])
             .set_coeff(0, y[i])
             .set_coeff(1, t1)
             .set_coeff(2, (3.0*m[i]-2.0*t1-t2)/(x[i+1]-x[i]))
             .set_coeff(3, (t1+t2-2.0*m[i])/((x[i+1]-x[i])*(x[i+1]-x[i])));
     }
    s._domain_pair.first = xy_map.cbegin()->first;
    s._domain_pair.second = xy_map.crbegin()->first;
    return s;
}

double spline::operator()(double x) const {
    if(_poly_map.empty())
        return 0;
    auto cit = _poly_map.cbegin();
    if(x >= _poly_map.crbegin()->first)
        cit = std::prev(_poly_map.cend());
    else if(x > _poly_map.cbegin()->first)
        cit = std::prev(_poly_map.lower_bound(x));
    return cit->second(x);
}

spline spline::operator-() const {
    return spline(*this)*(-1);
}

spline spline::operator+(double y) const {
    return spline(*this) += y;
}

spline spline::operator-(double y) const {
    return spline(*this) -= y;
}

spline spline::operator*(double y) const {
    return spline(*this) *= y;
}

spline spline::operator/(double y) const {
    return spline(*this) /= y;
}

spline& spline::operator+=(double y) {
    for(auto it = _poly_map.begin(); it != _poly_map.end(); it++)
        it->second += y;
    return *this;
}

spline& spline::operator-=(double y) {
    return (*this) += (-y);
}

spline& spline::operator*=(double y) {
    for(auto it = _poly_map.begin(); it != _poly_map.end(); it++)
        it->second *= y;
    return *this;
}

spline& spline::operator/=(double y) {
    return (*this) *= (1.0/y);
}

spline::const_iterator spline::cbegin() const {
    return _poly_map.cbegin();
}

spline::const_iterator spline::cend() const {
    return _poly_map.cend();
}

spline& spline::differentiate() {
    for(auto it = _poly_map.begin(); it != _poly_map.end(); it++)
        it->second.differentiate();
    return *this;
}

const std::pair<double,double>& spline::domain() const {
    return _domain_pair;
}

spline& spline::integrate(double x) {
    double sum = 0;
    for(auto it = _poly_map.begin(); it != _poly_map.end(); it++) {
        it->second.integrate(sum);
        if(std::next(it) != _poly_map.end())
            sum = it->second(std::next(it)->first);
    }
    return (*this) -= (*this)(x);
}

spline operator+(double y, const spline& s) {
    return s+y;
}

spline operator*(double y, const spline& s) {
    return s*y;
}

archive::ostream& operator<<(archive::ostream& oa, const spline& obj) {
    oa.put_float64(obj._domain_pair.first).put_float64(obj._domain_pair.second);
    oa.put_uint32(obj._poly_map.size());
    for(auto cit = obj._poly_map.cbegin(); cit != obj._poly_map.cend(); cit++) {
        oa.put_float64(cit->first);
        oa << cit->second;
    }
    return oa;
}

archive::istream& operator>>(archive::istream& ia, spline& obj) {
    ia.get_float64(obj._domain_pair.first).get_float64(obj._domain_pair.second);
    obj._poly_map.clear();
    uint32_t n;
    ia.get_uint32(n);
    for(uint32_t i = 0; i < n; i++) {
        double x;
        ia.get_float64(x);
        ia >> obj._poly_map[x];
    }
    return ia;
}
