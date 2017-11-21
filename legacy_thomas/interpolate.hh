/**
 * @file src/common/interpolate.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__INTERPOLATE__HEADER_INCLUDED
#define eSCATTER__COMMON__INTERPOLATE__HEADER_INCLUDED

#include <map>

/*!
 * Linear interpolation in one dimension.
 * Computational complexity is approx. O(log N).
 * warning: returns border values for evaluations outside domain.
 */
double interpolate(const std::map<double,double>& xy_map, double x);

/*!
 * Linear interpolation in two dimensions.
 * Computational complexity is approx. O(log N).
 * warning: returns border values for evaluations outside domain.
 */
double interpolate(const std::map<double,std::map<double,double>>& xyz_map, double x, double y);

#endif