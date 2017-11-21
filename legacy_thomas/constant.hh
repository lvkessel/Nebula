/*!
 * @file src/common/constant.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__CONSTANT__HEADER_INCLUDED
#define eSCATTER__COMMON__CONSTANT__HEADER_INCLUDED

namespace constant {

/*!
 * pi
 */
const double pi = 3.14159265359;
/*!
 * speed of light
 */
const double c = 299792458;
/*!
 * speed of light squared
 */
const double c_pow_2 = c*c;
/*!
 * elementary charge
 */
const double ec = 1.60217657e-19;
/*!
 * electric constant
 */
const double eps0 = 8.85418782e-12;
/*!
 * Planck's constant
 */
const double h = 6.62606957e-34;
/*!
 * reduced Planck's constant
 */
const double hbar = 1.05457173e-34;
/*!
 * Boltzmann's constant
 */
const double k = 1.3806488e-23;
/*!
 * electron mass
 */
const double me = 9.10938291e-31;
/*!
 * electron mass times speed of light squared
 */
const double me_c_pow_2 = me*c_pow_2;
/*!
 * magnetic constant
 */
const double mu0 = 1.25663706e-6;
/*!
 * Avogadro's number
 */
const double NA = 6.02214129e23;

}

#endif