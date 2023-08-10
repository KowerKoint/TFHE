#pragma once

#include <cmath>
#include <complex>
#include <cstdint>
#include <iterator>
#include <memory>

namespace TFHE {
int p = 0;
template <int N>
struct is_valid_int_bits {
    constexpr static bool value = N == 8 || N == 16 || N == 32;
};
template <int N>
constexpr bool is_valid_int_bits_v = is_valid_int_bits<N>::value;
template <int N, std::enable_if_t<is_valid_int_bits_v<N>>* = nullptr>
using SignedInt = std::conditional_t<N == 8, int8_t,
    std::conditional_t<N == 16, int16_t, int32_t>>;

class TorusValue {
    uint32_t _val;

public:
    constexpr static uint32_t MU = 1U << 29;  // 1/8

    constexpr TorusValue() : _val(0U) {}
    explicit TorusValue(double d) {
        d = fmod(d, 1.);
        if (d < 0.) d += 1.;
        d *= (1LL << 32);
        _val = d;
    }
    constexpr explicit TorusValue(bool b) : _val(b ? MU : -MU) {}
    constexpr static TorusValue from_raw_value(uint32_t val) {
        TorusValue ret;
        ret._val = val;
        return ret;
    }

    constexpr explicit operator double() const {
        double d = _val;
        d /= (1LL << 32);
        return d;
    }
    constexpr explicit operator bool() const { return !(_val >> 31 & 1U); }
    constexpr const uint32_t& get_raw_value() const { return _val; }
    constexpr uint32_t& get_raw_value() { return _val; }

    constexpr TorusValue operator+() const { return from_raw_value(_val); }
    constexpr TorusValue operator-() const { return from_raw_value(-_val); }

    constexpr TorusValue& operator+=(const TorusValue& rhs) {
        _val += rhs._val;
        return *this;
    }
    constexpr TorusValue operator+(const TorusValue& rhs) const {
        return TorusValue(*this) += rhs;
    }
    constexpr TorusValue& operator-=(const TorusValue& rhs) {
        _val -= rhs._val;
        return *this;
    }
    constexpr TorusValue operator-(const TorusValue& rhs) const {
        return TorusValue(*this) -= rhs;
    }

    template <typename Scalar,
        std::enable_if_t<std::is_scalar_v<Scalar>>* = nullptr>
    constexpr TorusValue& operator*=(Scalar rhs) {
        _val *= rhs;
        return *this;
    }
    template <typename Scalar,
        std::enable_if_t<std::is_scalar_v<Scalar>>* = nullptr>
    constexpr TorusValue operator*(Scalar rhs) const {
        return TorusValue(*this) *= rhs;
    }

    template <typename Scalar,
        std::enable_if_t<std::is_scalar_v<Scalar>>* = nullptr>
    friend constexpr TorusValue operator*(Scalar lhs, const TorusValue& rhs) {
        return TorusValue::from_raw_value(rhs._val * lhs);
    }

    constexpr bool operator==(const TorusValue& rhs) const {
        return _val == rhs._val;
    }
    constexpr bool operator!=(const TorusValue& rhs) const {
        return _val != rhs._val;
    }
};

template <typename Value, int N, typename VectorOrPolynomial>
class _VectorOrPolynomial {
protected:
    std::unique_ptr<Value[]> _val;

public:
    _VectorOrPolynomial() : _val{std::make_unique<Value[]>(N)} {}
    _VectorOrPolynomial(const std::unique_ptr<Value[]>& val)
        : _val{std::make_unique<Value[]>(N)} {
        std::copy(val.get(), val.get() + N, _val.get());
    }
    _VectorOrPolynomial(std::unique_ptr<Value[]>&& val)
        : _val{std::move(val)} {}

    Value* begin() { return _val.get(); }
    const Value* begin() const { return _val.get(); }
    Value* end() { return _val.get() + N; }
    const Value* end() const { return _val.get() + N; }

    const Value& operator[](int idx) const { return _val[idx]; }
    Value& operator[](int idx) { return _val[idx]; }

    VectorOrPolynomial operator+() const { return VectorOrPolynomial(*this); }
    VectorOrPolynomial operator-() const {
        VectorOrPolynomial ret(*this);
        for (int i = 0; i < N; i++) ret[i] = -_val[i];
        return ret;
    }
    VectorOrPolynomial& operator+=(const VectorOrPolynomial& rhs) {
        for (int i = 0; i < N; i++) _val[i] += rhs[i];
        return *(static_cast<VectorOrPolynomial*>(this));
    }
    VectorOrPolynomial operator+(const VectorOrPolynomial& rhs) const {
        return VectorOrPolynomial(_val) += rhs;
    }
    VectorOrPolynomial& operator-=(const VectorOrPolynomial& rhs) {
        for (int i = 0; i < N; i++) _val[i] -= rhs[i];
        return *(static_cast<VectorOrPolynomial*>(this));
    }
    VectorOrPolynomial operator-(const VectorOrPolynomial& rhs) const {
        return VectorOrPolynomial(_val) -= rhs;
    }

    bool operator==(const VectorOrPolynomial& rhs) const {
        for (int i = 0; i < N; i++)
            if (_val[i] != rhs[i]) return false;
        return true;
    }
    bool operator!=(const VectorOrPolynomial& rhs) const {
        for (int i = 0; i < N; i++)
            if (_val[i] != rhs[i]) return true;
        return false;
    }
};

template <typename Value, int N>
class Polynomial;

template <int N>
Polynomial<std::complex<double>, N> fft(
    const Polynomial<std::complex<double>, N>& a,
    const std::array<std::complex<double>, (size_t)N>& omega_pow) {
    auto get_log_2 = [](int n) {
        int n_bit = 0;
        while ((1 << n_bit) < n) n_bit++;
        return n_bit;
    };
    constexpr int N_BIT = get_log_2(N);
    static_assert((1 << N_BIT) == N);

    auto bit_rev = [&](int i) {
        int ret = 0;
        for (int j = 0; j < N_BIT; j++) {
            if (i & (1 << j)) ret |= 1 << (N_BIT - j - 1);
        }
        return ret;
    };

    Polynomial<std::complex<double>, N> ret = a;
    for (int i = 0; i < N_BIT; i++) {
        Polynomial<std::complex<double>, N> nret;
        for (int j = 0; j < 1 << i; j++) {
            int rev = bit_rev(j << 1);
            for (int k = 0; k < 1 << (N_BIT - i - 1); k++) {
                std::complex<double> even = ret[(j << (N_BIT - i)) + k];
                std::complex<double> odd =
                    ret[(j << (N_BIT - i)) + k + (1 << (N_BIT - i - 1))];
                odd *= omega_pow[rev];
                nret[(j << (N_BIT - i)) + k] = even + odd;
                nret[(j << (N_BIT - i)) + k + (1 << (N_BIT - i - 1))] =
                    even - odd;
            }
        }
        ret = std::move(nret);
    }
    {
        Polynomial<std::complex<double>, N> nret;
        for (int i = 0; i < N; i++) {
            nret[bit_rev(i)] = ret[i];
        }
        ret = std::move(nret);
    }
    return ret;
}

template <typename Value, int N>
class Polynomial : public _VectorOrPolynomial<Value, N, Polynomial<Value, N>> {
    using Base = _VectorOrPolynomial<Value, N, Polynomial>;
    friend Base;

public:
    Polynomial() : Base() {}
    Polynomial(const Polynomial& p) : Base(p._val){};
    Polynomial(Polynomial&& p) : Base(std::move(p._val)){};
    Polynomial(const std::unique_ptr<Value[]>& val) : Base{val} {}
    Polynomial(std::unique_ptr<Value[]>&& val) : Base{std::move(val)} {}

    Polynomial& operator=(const Polynomial& rhs) {
        std::copy(rhs.begin(), rhs.end(), this->begin());
        return *this;
    }
    Polynomial& operator=(Polynomial&& rhs) {
        if (this != &rhs) this->_val = std::move(rhs._val);
        return *this;
    }

    template <typename RHSValue>
    Polynomial<decltype(std::declval<Value>() * std::declval<RHSValue>()), N>
    operator*(const Polynomial<RHSValue, N>& rhs) const {
        using Result =
            decltype(std::declval<Value>() * std::declval<RHSValue>());
        if constexpr (N < 32) {
            Polynomial<Result, N> ret;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (i + j < N)
                        ret[i + j] += (*this)[i] * rhs[j];
                    else
                        ret[i + j - N] -= (*this)[i] * rhs[j];
                }
            }
            return ret;
        }
        auto expanded_n = [](int n) {
            int ret = 1;
            while (ret < n) ret <<= 1;
            return ret;
        };
        constexpr int M = expanded_n(N) >> 1;
        Polynomial<std::complex<double>, M> lhs_c;
        for (int i = 0; i < N; i++) {
            if (i < M)
                lhs_c[i].real((double)((*this)[i]));
            else
                lhs_c[i - M].imag((double)((*this)[i]));
        }
        Polynomial<std::complex<double>, M> rhs_c;
        for (int i = 0; i < N; i++) {
            if (i < M)
                rhs_c[i].real((double)(rhs[i]));
            else
                rhs_c[i - M].imag((double)(rhs[i]));
        }
        std::complex<double> w = std::polar(1., M_PI / (M * 2));
        std::complex<double> w_pow = 1.;
        for (int i = 0; i < M; i++) {
            lhs_c[i] *= w_pow;
            rhs_c[i] *= w_pow;
            w_pow *= w;
        }

        std::complex<double> omega = std::polar(1., -M_PI * 2 / M);
        std::array<std::complex<double>, M> omega_pow;
        omega_pow[0] = 1.;
        for (int i = 0; i < M - 1; i++) omega_pow[i + 1] = omega_pow[i] * omega;

        auto lhs_fft = fft(lhs_c, omega_pow);
        auto rhs_fft = fft(rhs_c, omega_pow);
        Polynomial<std::complex<double>, M> ret_fft;
        for (int i = 0; i < M; i++) {
            ret_fft[i] = lhs_fft[i] * rhs_fft[i];
        }
        std::reverse(omega_pow.begin() + 1, omega_pow.end());
        auto ret_c = fft(ret_fft, omega_pow);
        std::complex<double> w_inv = std::polar(1., -M_PI / (M * 2));
        std::complex<double> w_inv_pow = 1.;
        for (int i = 0; i < M; i++) {
            ret_c[i] *= 1. / M;
            ret_c[i] *= w_inv_pow;
            w_inv_pow *= w_inv;
        }
        Polynomial<Result, N> ret;
        for (int i = 0; i < N; i++) {
            if (i < M)
                ret[i] = (Result)ret_c[i].real();
            else
                ret[i] = (Result)ret_c[i - M].imag();
        }
        return ret;
    }
    template <typename RHSValue>
    Polynomial& operator*=(const Polynomial<RHSValue, N>& rhs) {
        swap((*this)._val, (*this * rhs)._val);
        return *this;
    }
    Polynomial multiply_x_exp(int n) const {
        n %= N * 2;
        if (n < 0) n += N * 2;
        Polynomial ret;
        for (int i = 0; i < N; i++) {
            if (i + n >= N * 2)
                ret[i + n - N * 2] = (*this)[i];
            else if (i + n >= N)
                ret[i + n - N] = -(*this)[i];
            else
                ret[i + n] = (*this)[i];
        }
        return ret;
    }
};

template <typename Value, int M, int N>
class Matrix;

template <typename Value, int N>
class Vector : public _VectorOrPolynomial<Value, N, Vector<Value, N>> {
    using Base = _VectorOrPolynomial<Value, N, Vector>;
    friend Base;

    friend Matrix<Value, N, N>;

public:
    Vector() : Base{} {}
    Vector(const Vector& v) : Base{v._val} {}
    Vector(Vector&& v) : Base{std::move(v._val)} {}
    Vector(const std::unique_ptr<Value[]>& val) : Base{val} {}
    Vector(std::unique_ptr<Value[]>&& val) : Base{std::move(val)} {}

    Vector& operator=(const Vector& rhs) {
        std::copy(rhs.begin(), rhs.end(), this->begin());
        return *this;
    }
    Vector& operator=(Vector&& rhs) {
        if (this != &rhs) this->_val = std::move(rhs._val);
        return *this;
    }

    template <typename RHSValue>
    decltype(std::declval<Value>() * std::declval<RHSValue>()) dot(
        const Vector<RHSValue, N>& rhs) const {
        using Result =
            decltype(std::declval<Value>() * std::declval<RHSValue>());
        Result ret{};
        for (int i = 0; i < N; i++) {
            ret += (*this)[i] * rhs[i];
        }
        return ret;
    }

    Vector operator*=(const Matrix<Value, N, N>& rhs) {
        *this = *this * rhs;
        return *this;
    }
};

template <typename Value, int M, int N>
class Matrix {
    friend Vector<Value, M>;
    friend Vector<Value, N>;

    Vector<Vector<Value, N>, M> _val;

public:
    Matrix() : _val{} {}
    Matrix(const Matrix& v) : _val{v._val} {}
    Matrix(Matrix&& v) : _val{std::move(v._val)} {}

    Matrix& operator=(const Matrix& rhs) {
        this->_val = rhs._val;
        return *this;
    }
    Matrix& operator=(Matrix&& rhs) {
        if (this != &rhs) this->_val = std::move(rhs._val);
        return *this;
    }

    auto begin() { return _val.begin(); }
    auto begin() const { return _val.begin(); }
    auto end() { return _val.end(); }
    auto end() const { return _val.end(); }

    const Vector<Value, N>& operator[](int idx) const { return _val[idx]; }
    Vector<Value, N>& operator[](int idx) { return _val[idx]; }

    Matrix operator+() const { return Matrix(*this); }
    Matrix operator-() const {
        Matrix ret(*this);
        ret._val = -ret._val;
        return ret;
    }
    Matrix operator+=(const Matrix& rhs) {
        _val += rhs._val;
        return *this;
    }
    Matrix operator+(const Matrix& rhs) const { return Matrix(*this) += rhs; }
    Matrix operator-=(const Matrix& rhs) {
        _val -= rhs._val;
        return *this;
    }
    Matrix operator-(const Matrix& rhs) const { return Matrix(*this) -= rhs; }

    template <typename RHSValue, int O>
    Matrix<decltype(std::declval<Value>() * std::declval<RHSValue>()), M, O>
    operator*(const Matrix<RHSValue, N, O>& rhs) const {
        using Result =
            decltype(std::declval<Value>() * std::declval<RHSValue>());
        Matrix<Result, M, O> ret;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < O; k++) {
                    ret[i][k] += (*this)[i][j] * rhs[j][k];
                }
            }
        }
        return ret;
    }
    template <typename RHSValue>
    Vector<decltype(std::declval<Value>() * std::declval<RHSValue>()), M>
    operator*(const Vector<RHSValue, N>& rhs) const {
        using Result =
            decltype(std::declval<Value>() * std::declval<RHSValue>());
        const Vector<Result, M> ret;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                ret[i] += (*this)[i][j] * rhs[j];
            }
        }
        return ret;
    }
    template <typename LHSValue>
    friend Vector<decltype(std::declval<LHSValue>() * std::declval<Value>()), N>
    operator*(const Vector<LHSValue, M>& lhs, const Matrix<Value, M, N>& rhs) {
        using Result =
            decltype(std::declval<LHSValue>() * std::declval<Value>());
        Vector<Result, N> ret;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                ret[j] += lhs[i] * rhs[i][j];
            }
        }
        return ret;
    }
    template <typename RHSValue>
    Matrix& operator*=(const Matrix<RHSValue, N, N>& rhs) {
        swap((*this)._val, (*this * rhs)._val);
        return *this;
    }

    bool operator==(const Matrix& rhs) const { return _val == rhs._val; }
    bool operator!=(const Matrix& rhs) const { return _val != rhs._val; }
};
}  // namespace TFHE