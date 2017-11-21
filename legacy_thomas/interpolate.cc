/**
 * @file src/common/interpolate.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "interpolate.hh"

double interpolate(const std::map<double,double>& xy_map, double x) {
    if(xy_map.empty())
        return 0;
    if(xy_map.size() == 1)
        return xy_map.cbegin()->second;

    auto cit = xy_map.upper_bound(x);
    if(cit == xy_map.cbegin())
        cit++;
    else if(cit == xy_map.cend())
        cit--;

    const double x2 = cit->first;
    const double y2 = cit->second;
    cit--;
    const double x1 = cit->first;
    const double y1 = cit->second;

    const double t = (x-x1)/(x2-x1);
    return (1.0-t)*y1 + t*y2;
}

double interpolate(const std::map<double,std::map<double,double>>& xyz_map, double x, double y) {
    if(xyz_map.empty())
        return 0;
    if(xyz_map.size() == 1)
        return interpolate(xyz_map.cbegin()->second, y);

    auto cit = xyz_map.upper_bound(x);
    if(cit == xyz_map.cbegin())
        cit++;
    else if(cit == xyz_map.cend())
        cit--;

    const double x2 = cit->first;
    const double z2 = interpolate(cit->second, y);
    cit--;
    const double x1 = cit->first;
    const double z1 = interpolate(cit->second, y);

    const double t = (x-x1)/(x2-x1);
    return (1.0-t)*z1 + t*z2;
}
