/**
 * @file src/common/polynomial.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "polynomial.hh"
#include <cinttypes>
#include <cmath>

polynomial::polynomial(double x0) : polynomial({}, x0) {
}

polynomial::polynomial(const std::initializer_list<double>& ai, double x0) {
    for(auto it = ai.begin(); it != ai.end(); it++)
        _coeff_vec.push_back(*it);
    _origin = x0;
}

polynomial polynomial::operator-() const {
    return polynomial(*this)*(-1);
}

polynomial polynomial::operator+(double y) const {
    return polynomial(*this) += y;
}

polynomial polynomial::operator-(double y) const {
    return polynomial(*this) -= y;
}

polynomial polynomial::operator*(double y) const {
    return polynomial(*this) *= y;
}

polynomial polynomial::operator/(double y) const {
    return polynomial(*this) /= y;
}

polynomial& polynomial::operator+=(double y) {
    _coeff_vec[0] += y;
    return *this;
}

polynomial& polynomial::operator-=(double y) {
    return (*this) += (-y);
}

polynomial& polynomial::operator*=(double y) {
    for(int i = order(); i >= 0; i--)
        _coeff_vec[i] *= y;
    return *this;
}

polynomial& polynomial::operator/=(double y) {
    return (*this) *= (1.0/y);
}

double polynomial::operator()(double x) const {
    if(_coeff_vec.empty())
        return 0;
    double y = 0;
    for(int i = order(); i > 0; i--) {
        y += _coeff_vec[i];
        y *= x-_origin;
    }
    y += _coeff_vec[0];
    return y;
}

double polynomial::operator[](int i) const {
    if(_coeff_vec.empty() || (i < 0) || (i > order()))
        return 0;
    return _coeff_vec[i];
}

polynomial& polynomial::differentiate() {
    if(order() == 0) {
        _coeff_vec[0] = 0;
    } else {
        for(int i = 1; i <= order(); i++)
            _coeff_vec[i] *= i;
        _coeff_vec.erase(_coeff_vec.begin());
    }
    return *this;
}

polynomial& polynomial::integrate(double a0) {
    for(int i = 0; i <= order(); i++)
        _coeff_vec[i] /= i+1;
    _coeff_vec.insert(_coeff_vec.begin(), a0);
    return *this;
}

int polynomial::order() const {
    if(_coeff_vec.empty())
        return 0;
    return _coeff_vec.size()-1;
}

double polynomial::origin() const {
    return _origin;
}

polynomial& polynomial::set_coeff(int i, double ai) {
    if(i < 0)
        return *this;
    if((i > order()) && (ai == 0))
        return *this;
    if(_coeff_vec.empty() || (i > order()))
        _coeff_vec.resize(i+1, 0);
    _coeff_vec[i] = ai;
    return *this;
}

polynomial& polynomial::transform(double x0) {
    auto bico = [](int n, int k)->uint64_t {
        uint64_t cnk = 1;
        for(int i = 1; i <= k; i++) {
            cnk *= n+1-i;
            cnk /= i;
        }
        return cnk;
    };
    for(int i = 0; i <= order(); i++) {
        polynomial ai(_origin);
        for(int j = 0; j <= order()-i; j++)
            ai.set_coeff(j, _coeff_vec[i+j]*bico(i+j, i));
        _coeff_vec[i] = ai(x0);
    }
    _origin = x0;
    return *this;
}

polynomial operator+(double y, const polynomial& p) {
    return p+y;
}

polynomial operator*(double y, const polynomial& p) {
    return p*y;
}

archive::ostream& operator<<(archive::ostream& oa, const polynomial& obj) {
    oa.put_float64(obj.origin());
    const uint32_t n = obj.order();
    oa.put_uint32(n);
    for(uint32_t i = 0; i <= n; i++)
        oa.put_float64(obj[i]);
    return oa;
}

archive::istream& operator>>(archive::istream& ia, polynomial& obj) {
    double x0;
    ia.get_float64(x0);
    obj = polynomial(x0);
    uint32_t n;
    ia.get_uint32(n);
    for(uint32_t i = 0; i <= n; i++) {
        double ai;
        ia.get_float64(ai);
        obj.set_coeff(i, ai);
    }
    return ia;
}
