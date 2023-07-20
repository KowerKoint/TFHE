#pragma once

#include <cmath>
#include <cstdint>

namespace TFHE {
class TorusValue {
    uint32_t _val;

public:
    constexpr static uint32_t MU = 1U << 29;  // 1/8

    TorusValue() = default;
    constexpr TorusValue(uint32_t u) : _val(u) {}
    TorusValue(double d) : _val(0.) {
        d = fmod(d, 1.);
        if (d < 0.) d += 1.;
        d *= (1LL << 32);
        _val = d;
    }
    constexpr TorusValue(bool b) : _val(0.) {
        if (b) {
            _val = MU;
        } else {
            _val = -MU;
        }
    }

    constexpr operator double() const {
        double d = _val;
        d /= (1LL << 32);
        return d;
    }
    constexpr operator bool() const { return !(_val >> 31 & 1U); }

    constexpr TorusValue operator+() const { return TorusValue(_val); }
    constexpr TorusValue operator-() { return TorusValue(-_val); }

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

    constexpr TorusValue& operator*=(bool rhs) {
        if (!rhs) _val = 0;
        return *this;
    }
    constexpr TorusValue operator*(bool rhs) const {
        return TorusValue(*this) *= rhs;
    }
};

template <int N>
class BitVector {
    constexpr static int M = (N + 63) / 64;
    std::array<uint64_t, M> _val;

public:
    struct Reference {
        BitVector& container;
        int idx;

        constexpr operator bool() { return container.get(idx); }
        constexpr bool operator!() { return !container.get(idx); }
        constexpr Reference operator=(bool b) {
            container.set(idx, b);
            return *this;
        }
        constexpr Reference operator&=(bool b) {
            bool nb = container.get(idx) & b;
            container.set(idx, nb);
            return *this;
        }
        constexpr Reference operator|=(bool b) {
            bool nb = container.get(idx) | b;
            container.set(idx, nb);
            return *this;
        }
        constexpr Reference operator^=(bool b) {
            bool nb = container.get(idx) ^ b;
            container.set(idx, nb);
            return *this;
        }
    };

    BitVector() = default;
    constexpr BitVector(const std::array<uint64_t, M>& val) : _val{val} {
        constexpr int NUM_UNUSE = M * 64 - N;
        _val.back() &= (1ULL << (64 - NUM_UNUSE)) - 1;
    }
    constexpr BitVector(std::array<uint64_t, N>&& val) : _val{std::move(val)} {
        constexpr int NUM_UNUSE = M * 64 - N;
        _val.back() &= (1ULL << (64 - NUM_UNUSE)) - 1;
    }

    constexpr bool get(int idx) const {
        return (_val[idx >> 6] >> (idx & 0x3f) & 1);
    }
    constexpr void set(int idx, bool b) {
        if (b) {
            _val[idx >> 6] |= 1ULL << (idx & 0x3f);
        } else {
            _val[idx >> 6] &= ~(1ULL << (idx & 0x3f));
        }
    }
    constexpr bool operator[](int idx) const { return get(idx); }
    constexpr Reference operator[](int idx) { return Reference{*this, idx}; }
};

// TODO: use SIMD
template <int N>
class TorusVector {
    std::array<TorusValue, N> _val;

public:
    TorusVector() = default;
    template <typename... Args>
    constexpr TorusVector(Args... args) : _val{args...} {}
    constexpr TorusVector(const std::array<TorusValue, N>& val) : _val{val} {}
    constexpr TorusVector(std::array<TorusValue, N>&& val)
        : _val{std::move(val)} {}

    constexpr TorusValue operator[](int idx) const { return _val[idx]; }
    constexpr TorusValue& operator[](int idx) { return _val[idx]; }

    using Iterator = typename std::array<TorusValue, N>::iterator;
    using ConstIterator = typename std::array<TorusValue, N>::const_iterator;
    using ReverseIterator =
        typename std::array<TorusValue, N>::reverse_iterator;
    using ConstReverseIterator =
        typename std::array<TorusValue, N>::const_reverse_iterator;
    constexpr Iterator begin() { return _val.begin(); }
    constexpr ConstIterator begin() const { return _val.begin(); }
    constexpr Iterator end() { return _val.end(); }
    constexpr ConstIterator end() const { return _val.end(); }
    constexpr ConstIterator cbegin() const { return _val.cbegin(); }
    constexpr ConstIterator cend() const { return _val.cend(); }
    constexpr Iterator rbegin() { return _val.rbegin(); }
    constexpr ConstIterator rbegin() const { return _val.rbegin(); }
    constexpr Iterator rend() { return _val.rend(); }
    constexpr ConstIterator rend() const { return _val.rend(); }
    constexpr ConstIterator crbegin() const { return _val.crbegin(); }
    constexpr ConstIterator crend() const { return _val.crend(); }

    constexpr TorusVector operator+() const { return TorusVector(*this); }
    constexpr TorusVector operator-() const {
        TorusVector ret(*this);
        for (int i = 0; i < N; i++) ret[i] = -_val[i];
        return ret;
    }
    constexpr TorusVector& operator+=(const TorusVector& rhs) {
        for (int i = 0; i < N; i++) _val[i] += rhs[i];
    }
    constexpr TorusVector operator+(const TorusVector& rhs) const {
        return TorusVector(*this) += rhs;
    }
    constexpr TorusVector& operator-=(const TorusVector& rhs) {
        for (int i = 0; i < N; i++) _val[i] -= rhs[i];
    }
    constexpr TorusVector operator-(const TorusVector& rhs) const {
        return TorusVector(*this) -= rhs;
    }

    constexpr TorusValue dot(const BitVector<N>& rhs) const {
        TorusValue ret = 0.;
        for (int i = 0; i < N; i++) {
            if (rhs[i]) ret += _val[i];
        }
        return ret;
    }
};

}  // namespace TFHE