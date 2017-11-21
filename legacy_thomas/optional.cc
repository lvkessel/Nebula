/*!
 * @file src/common/optional.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "optional.hh"

archive::ostream& operator<<(archive::ostream& oa, const optional<double>& obj) {
    oa.put_bool(obj.is_defined());
    if(obj.is_defined())
        oa.put_float64(obj());
    return oa;
}

archive::istream& operator>>(archive::istream& ia, optional<double>& obj) {
    bool is_defined;
    ia.get_bool(is_defined);
    if(is_defined) {
        double value;
        ia.get_float64(value);
        obj = value;
    } else
        obj.clear();
    return ia;
}
