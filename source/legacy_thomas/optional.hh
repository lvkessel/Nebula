/*!
 * @file src/common/optional.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__OPTIONAL__HEADER_INCLUDED
#define eSCATTER__COMMON__OPTIONAL__HEADER_INCLUDED

#include "archive.hh"

/*!
 *
 */
template<typename T>
class optional {
public:
    optional() = default;
    optional(const T&);
    optional(const optional&);
    ~optional();
    optional& operator=(const optional&);
    optional& operator=(const T&);
    /*!
      * @warning throws an exception if the value is not defined.
      */
    const T& operator()() const;
    /*!
      * @return
      */
    bool is_defined() const;
    /*!
      * @return
      */
    optional& clear();
private:
    T* _value_p = nullptr;
};

archive::ostream& operator<<(archive::ostream&, const optional<double>&);
archive::istream& operator>>(archive::istream&, optional<double>&);

#include "optional.inl"

#endif
