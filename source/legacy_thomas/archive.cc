/**
 * @file src/common/archive.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "archive.hh"
#include <stdexcept>

namespace archive {

istream::istream(std::istream& is) {
    _streambuf_p = is.rdbuf();
}

data_type_enum istream::peek() {
    const data_type_enum data_type = static_cast<data_type_enum>(_streambuf_p->sgetc());
    switch(data_type) {
        case BOOL:
        case INT8:
        case UINT8:
        case INT16:
        case UINT16:
        case INT32:
        case UINT32:
        case INT64:
        case UINT64:
        case FLOAT32:
        case FLOAT64:
        case STRING:
        case BLOB:
            return data_type;
    }
    throw std::runtime_error("unknown data type");
}

istream& istream::get_bool(bool& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != BOOL)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_int8(int8_t& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != INT8)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_uint8(uint8_t& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != UINT8)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_int16(int16_t& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != INT16)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_uint16(uint16_t& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != UINT16)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_int32(int32_t& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != INT32)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_uint32(uint32_t& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != UINT32)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_int64(int64_t& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != INT64)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_uint64(uint64_t& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != UINT64)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_float32(float& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != FLOAT32)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_float64(double& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != FLOAT64)
        throw std::runtime_error("unexpected data type");
    _streambuf_p->sgetn(reinterpret_cast<char*>(&x), sizeof(x));
    return *this;
}

istream& istream::get_string(std::string& x) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != STRING)
        throw std::runtime_error("unexpected data type");
    uint64_t byte_count;
    _streambuf_p->sgetn(reinterpret_cast<char*>(&byte_count), sizeof(byte_count));
    x.assign(byte_count, 0);
    for(uint32_t i = 0; i < byte_count; i++)
        x[i] = _streambuf_p->sbumpc();
    return *this;
}

istream& istream::get_blob(void* data, uint64_t size) {
    if(static_cast<data_type_enum>(_streambuf_p->sbumpc()) != BLOB)
        throw std::runtime_error("unexpected data type");
    uint64_t byte_count;
    _streambuf_p->sgetn(reinterpret_cast<char*>(&byte_count), sizeof(byte_count));
    if(byte_count != size)
        throw std::runtime_error("blob size mismatch");
    _streambuf_p->sgetn(reinterpret_cast<char*>(data), byte_count);
    return *this;
}

ostream::ostream(std::ostream& os) {
    _streambuf_p = os.rdbuf();
}

ostream& ostream::put_bool(bool x) {
    _streambuf_p->sputc(BOOL);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_int8(int8_t x) {
    _streambuf_p->sputc(INT8);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_uint8(uint8_t x) {
    _streambuf_p->sputc(UINT8);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_int16(int16_t x) {
    _streambuf_p->sputc(INT16);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_uint16(uint16_t x) {
    _streambuf_p->sputc(UINT16);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_int32(int32_t x) {
    _streambuf_p->sputc(INT32);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_uint32(uint32_t x) {
    _streambuf_p->sputc(UINT32);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_int64(int64_t x) {
    _streambuf_p->sputc(INT64);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_uint64(uint64_t x) {
    _streambuf_p->sputc(UINT64);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_float32(float x) {
    _streambuf_p->sputc(FLOAT32);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_float64(double x) {
    _streambuf_p->sputc(FLOAT64);
    _streambuf_p->sputn(reinterpret_cast<const char*>(&x), sizeof(x));
    return *this;
}

ostream& ostream::put_string(const std::string& x) {
    _streambuf_p->sputc(STRING);
    const uint64_t byte_count = x.size();
    _streambuf_p->sputn(reinterpret_cast<const char*>(&byte_count), sizeof(byte_count));
    _streambuf_p->sputn(x.data(), byte_count);
    return *this;
}

ostream& ostream::put_blob(const void* data, uint64_t size) {
    _streambuf_p->sputc(BLOB);
    const uint64_t byte_count = size;
    _streambuf_p->sputn(reinterpret_cast<const char*>(&byte_count), sizeof(byte_count));
    _streambuf_p->sputn(reinterpret_cast<const char*>(data), byte_count);
    return *this;
}

}
