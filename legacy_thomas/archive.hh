/*!
 * @file src/common/archive.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__ARCHIVE__HEADER_INCLUDED
#define eSCATTER__COMMON__ARCHIVE__HEADER_INCLUDED

#include <cinttypes>
#include <iostream>
#include <string>

namespace archive {

/*!
 *
 */
enum data_type_enum : int8_t {
    BOOL,
    INT8, UINT8,
    INT16, UINT16,
    INT32, UINT32,
    INT64, UINT64,
    FLOAT32, FLOAT64,
    STRING, BLOB
};

/*!
 *
 */
class istream {
public:
    /*!
      * @param[in,out] is
      */
    istream(std::istream& is);
    istream(const istream&) = delete;
    istream& operator=(const istream&) = delete;
    data_type_enum peek();
    istream& get_bool(bool& x);
    istream& get_int8(int8_t& x);
    istream& get_uint8(uint8_t& x);
    istream& get_int16(int16_t& x);
    istream& get_uint16(uint16_t& x);
    istream& get_int32(int32_t& x);
    istream& get_uint32(uint32_t& x);
    istream& get_int64(int64_t& x);
    istream& get_uint64(uint64_t& x);
    istream& get_float32(float& x);
    istream& get_float64(double& x);
    istream& get_string(std::string& x);
    istream& get_blob(void* data, uint64_t size);
private:
    std::streambuf* _streambuf_p;
};

/*!
 *
 */
class ostream {
public:
    /*!
      * @param[in,out] os
      */
    ostream(std::ostream& os);
    ostream(const ostream&) = delete;
    ostream& operator=(const ostream&) = delete;
    ostream& put_bool(bool x);
    ostream& put_int8(int8_t x);
    ostream& put_uint8(uint8_t x);
    ostream& put_int16(int16_t x);
    ostream& put_uint16(uint16_t x);
    ostream& put_int32(int32_t x);
    ostream& put_uint32(uint32_t x);
    ostream& put_int64(int64_t x);
    ostream& put_uint64(uint64_t x);
    ostream& put_float32(float x);
    ostream& put_float64(double x);
    ostream& put_string(const std::string& x);
    ostream& put_blob(const void* data, uint64_t size);
private:
    std::streambuf* _streambuf_p;
};

}

#endif