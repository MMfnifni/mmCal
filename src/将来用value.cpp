#pragma once
#include <cassert>
#include <complex>
#include <new>
#include <string>
#include <utility>

namespace mm::cal {

 struct InvalidValue {};

 struct Value {
   enum class Kind : unsigned char { Real, Complex, String, Invalid };

   Kind kind = Kind::Invalid;

   union Storage {
     double real;
     std::complex<double> complex;
     std::string str;

     Storage() {}
     ~Storage() {}
   } data;

   /* ============================
      ctor / dtor
      ============================ */

   Value() noexcept : kind(Kind::Invalid) {}

   Value(InvalidValue) noexcept : kind(Kind::Invalid) {}

   Value(double x) noexcept : kind(Kind::Real) { data.real = x; }

   Value(const std::complex<double> &c) : kind(Kind::Complex) { new (&data.complex) std::complex<double>(c); }

   Value(std::complex<double> &&c) noexcept : kind(Kind::Complex) { new (&data.complex) std::complex<double>(std::move(c)); }

   Value(const std::string &s) : kind(Kind::String) { new (&data.str) std::string(s); }

   Value(std::string &&s) noexcept : kind(Kind::String) { new (&data.str) std::string(std::move(s)); }

   Value(const char *s) : kind(Kind::String) { new (&data.str) std::string(s); }

   ~Value() { destroy(); }

   /* ============================
      copy / move
      ============================ */

   Value(const Value &other) { copyFrom(other); }

   Value(Value &&other) noexcept { moveFrom(std::move(other)); }

   Value &operator=(const Value &other) {
    if (this == &other) return *this;
    destroy();
    copyFrom(other);
    return *this;
   }

   Value &operator=(Value &&other) noexcept {
    if (this == &other) return *this;
    destroy();
    moveFrom(std::move(other));
    return *this;
   }

   /* ============================
      helpers
      ============================ */

   bool isReal() const noexcept { return kind == Kind::Real; }
   bool isComplex() const noexcept { return kind == Kind::Complex; }
   bool isString() const noexcept { return kind == Kind::String; }
   bool isInvalid() const noexcept { return kind == Kind::Invalid; }

   double &asReal() {
    assert(isReal());
    return data.real;
   }
   const double &asReal() const {
    assert(isReal());
    return data.real;
   }

   std::complex<double> &asComplex() {
    assert(isComplex());
    return data.complex;
   }
   const std::complex<double> &asComplex() const {
    assert(isComplex());
    return data.complex;
   }

   std::string &asString() {
    assert(isString());
    return data.str;
   }
   const std::string &asString() const {
    assert(isString());
    return data.str;
   }

  private:
   void destroy() {
    switch (kind) {
     case Kind::Complex: data.complex.~complex(); break;
     case Kind::String: data.str.~basic_string(); break;
     default: break;
    }
    kind = Kind::Invalid;
   }

   void copyFrom(const Value &other) {
    kind = other.kind;
    switch (other.kind) {
     case Kind::Real: data.real = other.data.real; break;
     case Kind::Complex: new (&data.complex) std::complex<double>(other.data.complex); break;
     case Kind::String: new (&data.str) std::string(other.data.str); break;
     case Kind::Invalid: break;
    }
   }

   void moveFrom(Value &&other) {
    kind = other.kind;
    switch (other.kind) {
     case Kind::Real: data.real = other.data.real; break;
     case Kind::Complex: new (&data.complex) std::complex<double>(std::move(other.data.complex)); break;
     case Kind::String: new (&data.str) std::string(std::move(other.data.str)); break;
     case Kind::Invalid: break;
    }
    other.destroy();
   }
 };

} // namespace mm::cal
