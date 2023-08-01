#pragma once

#include <cmath>
#include <cstdint>
#include <iterator>

namespace TFHE {
class TorusValue {
    uint32_t _val;

public:
    constexpr static uint32_t MU = 1U << 29;  // 1/8

    TorusValue() : _val(0U) {}
    constexpr TorusValue(uint32_t u) : _val(u) {}
    TorusValue(double d) {
        d = fmod(d, 1.);
        if (d < 0.) d += 1.;
        d *= (1LL << 32);
        _val = d;
    }
    constexpr TorusValue(bool b) : _val(b ? MU : -MU) {}

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
        if (!rhs) _val = 0U;
        return *this;
    }
    constexpr TorusValue operator*(bool rhs) const {
        return TorusValue(*this) *= rhs;
    }

    constexpr bool operator==(const TorusValue& rhs) const {
        return _val == rhs._val;
    }
    constexpr bool operator!=(const TorusValue& rhs) const {
        return _val != rhs._val;
    }
};

template <int N, class BitVector>
class _BitVector {
protected:
    constexpr static int M = (N + 63) / 64;
    std::array<uint64_t, M> _val;

public:
    class Reference {
        _BitVector& container;
        int idx;

    public:
        constexpr Reference(_BitVector& container, int idx)
            : container(container), idx(idx) {}

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

    template <typename ContainerReference, bool Reverse>
    class IteratorBase
        : public std::iterator<std::random_access_iterator_tag, bool> {
        ContainerReference container;
        int idx;

    public:
        constexpr IteratorBase(ContainerReference container, int idx)
            : container(container), idx(idx) {}

        constexpr auto operator*() const { return container[idx]; }

        constexpr IteratorBase& operator++() {
            if constexpr (!Reverse)
                ++idx;
            else
                --idx;
            return *this;
        }
        constexpr IteratorBase operator++(int) {
            IteratorBase tmp{*this};
            if constexpr (!Reverse)
                ++idx;
            else
                --idx;
            return tmp;
        }
        constexpr IteratorBase& operator--() {
            if constexpr (!Reverse)
                --idx;
            else
                ++idx;
            return *this;
        }
        constexpr IteratorBase operator--(int) {
            IteratorBase tmp{*this};
            if constexpr (!Reverse)
                --idx;
            else
                ++idx;
            return tmp;
        }

        constexpr IteratorBase& operator+=(int rhs) {
            if constexpr (!Reverse)
                idx += rhs;
            else
                idx -= rhs;
            return *this;
        }
        constexpr IteratorBase operator+(int rhs) const {
            return IteratorBase(container, idx + rhs);
        }
        constexpr IteratorBase& operator-=(int rhs) {
            if constexpr (!Reverse)
                idx -= rhs;
            else
                idx += rhs;
            return *this;
        }
        constexpr IteratorBase operator-(int rhs) const {
            return IteratorBase{container, idx - rhs};
        }

        constexpr bool operator==(const IteratorBase& rhs) const {
            return idx == rhs.idx;
        }
        constexpr bool operator!=(const IteratorBase& rhs) const {
            return idx != rhs.idx;
        }
        constexpr bool operator<(const IteratorBase& rhs) const {
            if constexpr (!Reverse)
                return idx < rhs.idx;
            else
                return idx > rhs.idx;
        }
        constexpr bool operator<=(const IteratorBase& rhs) const {
            if constexpr (!Reverse)
                return idx <= rhs.idx;
            else
                return idx >= rhs.idx;
        }
        constexpr bool operator>(const IteratorBase& rhs) const {
            if constexpr (!Reverse)
                return idx > rhs.idx;
            else
                return idx < rhs.idx;
        }
        constexpr bool operator>=(const IteratorBase& rhs) const {
            if constexpr (!Reverse)
                return idx >= rhs.idx;
            else
                return idx <= rhs.idx;
        }
    };

    using Iterator = IteratorBase<_BitVector&, false>;
    using ConstIterator = IteratorBase<const _BitVector&, false>;
    using ReverseIterator = IteratorBase<BitVector&, true>;
    using ConstReverseIterator = IteratorBase<const BitVector&, true>;

    _BitVector() = default;
    constexpr _BitVector(const std::array<uint64_t, M>& val) : _val{val} {
        constexpr int NUM_UNUSE = M * 64 - N;
        if (NUM_UNUSE) _val.back() &= (1ULL << (64 - NUM_UNUSE)) - 1;
    }
    constexpr _BitVector(std::array<uint64_t, N>&& val) : _val{std::move(val)} {
        constexpr int NUM_UNUSE = M * 64 - N;
        if (NUM_UNUSE) _val.back() &= (1ULL << (64 - NUM_UNUSE)) - 1;
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
    constexpr Iterator begin() { return Iterator{*this, 0}; }
    constexpr Iterator end() { return Iterator{*this, N - 1}; }
    constexpr ConstIterator begin() const { return ConstIterator{*this, 0}; }
    constexpr ConstIterator end() const { return ConstIterator{*this, N}; }
    constexpr ConstIterator cbegin() const { return ConstIterator{*this, 0}; }
    constexpr ConstIterator cend() const { return ConstIterator{*this, N}; }
    constexpr Iterator rbegin() { return ReverseIterator{*this, N - 1}; }
    constexpr Iterator rend() { return ReverseIterator{*this, -1}; }
    constexpr ConstReverseIterator rbegin() const {
        return ConstReverseIterator{*this, N - 1};
    }
    constexpr ConstReverseIterator rend() const {
        return ConstReverseIterator{*this, -1};
    }
    constexpr ConstReverseIterator crbegin() const {
        return ConstReverseIterator{*this, N - 1};
    }
    constexpr ConstReverseIterator crend() const {
        return ConstReverseIterator{*this, -1};
    }

    constexpr bool operator[](int idx) const { return get(idx); }
    constexpr Reference operator[](int idx) { return Reference{*this, idx}; }

    constexpr bool operator==(const BitVector& rhs) const {
        return _val == rhs._val;
    }
    constexpr bool operator!=(const BitVector& rhs) const {
        return _val != rhs._val;
    }
};

template <typename Value, int N, typename VectorOrPolynomial>
class _VectorOrPolynomial {
protected:
    std::array<Value, N> _val;

public:
    _VectorOrPolynomial() = default;
    template <typename... Args>
    constexpr _VectorOrPolynomial(Args... args) : _val{args...} {}
    constexpr _VectorOrPolynomial(const std::array<Value, N>& val)
        : _val{val} {}
    constexpr _VectorOrPolynomial(std::array<Value, N>&& val)
        : _val{std::move(val)} {}

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

    constexpr Value operator[](int idx) const { return _val[idx]; }
    constexpr Value& operator[](int idx) { return _val[idx]; }

    constexpr VectorOrPolynomial operator+() const {
        return VectorOrPolynomial(*this);
    }
    constexpr VectorOrPolynomial operator-() const {
        VectorOrPolynomial ret(*this);
        for (int i = 0; i < N; i++) ret[i] = -_val[i];
        return ret;
    }
    constexpr VectorOrPolynomial operator+=(const VectorOrPolynomial& rhs) {
        for (int i = 0; i < N; i++) _val[i] += rhs[i];
        return *this;
    }
    constexpr VectorOrPolynomial operator+(
        const VectorOrPolynomial& rhs) const {
        return VectorOrPolynomial(*this) += rhs;
    }
    constexpr VectorOrPolynomial operator-=(const VectorOrPolynomial& rhs) {
        for (int i = 0; i < N; i++) _val[i] -= rhs[i];
        return *this;
    }
    constexpr VectorOrPolynomial operator-(
        const VectorOrPolynomial& rhs) const {
        return VectorOrPolynomial(*this) -= rhs;
    }

    constexpr bool operator==(const VectorOrPolynomial& rhs) const {
        for (int i = 0; i < 0; i++)
            if (_val[i] != rhs[i]) return false;
        return true;
    }
    constexpr bool operator!=(const VectorOrPolynomial& rhs) const {
        for (int i = 0; i < 0; i++)
            if (_val[i] != rhs[i]) return true;
        return false;
    }
};

template <typename Value, int N>
class Polynomial : public _VectorOrPolynomial<Value, N, Polynomial<Value, N>> {
    using Base = _VectorOrPolynomial<Value, N, Polynomial>;
    friend Base;

public:
    Polynomial() = default;
    constexpr Polynomial(const Polynomial& p) : Base(p._val){};
    constexpr Polynomial(Polynomial&& p) : Base(std::move(p._val)){};
    template <typename... Args>
    constexpr Polynomial(Args... args) : Base{args...} {}
    constexpr Polynomial(const std::array<Value, N>& val) : Base{val} {}
    constexpr Polynomial(std::array<Value, N>&& val) : Base{std::move(val)} {}

    constexpr Polynomial& operator=(const Polynomial& rhs) {
        this->_val = rhs._val;
        return *this;
    }
    constexpr Polynomial& operator=(Polynomial&& rhs) {
        this->_val = std::move(rhs._val);
        return *this;
    }

    template <typename RHSValue, typename FAKE = void>
    constexpr Polynomial& operator*=(const Polynomial<RHSValue, N>& rhs) {
        std::array<Value, N> ret;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i + j < N)
                    ret[i + j] += (*this)[i] * rhs[j];
                else
                    ret[i + j - N] -= (*this)[i] * rhs[j];
            }
        }
        this->_val.swap(ret);
        return *this;
    }
    template <typename RHSValue>
    constexpr Polynomial operator*(const Polynomial<RHSValue, N>& rhs) const {
        return Polynomial(*this) *= rhs;
    }
};

template <int N>
class Polynomial<bool, N> : public _BitVector<N, Polynomial<bool, N>> {
    using Base = _BitVector<N, Polynomial>;
    friend Base;

public:
    Polynomial() = default;
    constexpr Polynomial(const Polynomial& p) : Base{p._val} {}
    constexpr Polynomial(Polynomial&& p) : Base{std::move(p._val)} {}
    constexpr Polynomial(const std::array<uint64_t, Base::M>& val)
        : Base{val} {}
    constexpr Polynomial(std::array<uint64_t, Base::M>&& val)
        : Base{std::move(val)} {}

    constexpr Polynomial& operator=(const Polynomial& rhs) {
        this->_val = rhs._val;
        return *this;
    }
    constexpr Polynomial& operator=(Polynomial&& rhs) {
        this->_val = std::move(rhs._val);
        return *this;
    }
};

template <typename Value, int N>
class Vector : public _VectorOrPolynomial<Value, N, Vector<Value, N>> {
    using Base = _VectorOrPolynomial<Value, N, Vector>;
    friend Base;

public:
    Vector() = default;
    constexpr Vector(const Vector& v) : Base{v._val} {}
    constexpr Vector(Vector&& v) : Base{std::move(v._val)} {}
    template <typename... Args>
    constexpr Vector(Args... args) : Base{args...} {}
    constexpr Vector(const std::array<Value, N>& val) : Base{val} {}
    constexpr Vector(std::array<Value, N>&& val) : Base{std::move(val)} {}

    template <typename RHSValue, typename FAKE = void>
    constexpr Value dot(const Vector<RHSValue, N>& rhs) const {
        Value ret{};
        for (int i = 0; i < N; i++) {
            ret += (*this)[i] * rhs[i];
        }
        return ret;
    }
};

template <int N>
class Vector<bool, N> : public _BitVector<N, Vector<bool, N>> {
    using Base = _BitVector<N, Vector>;
    friend Base;

public:
    Vector() = default;
    constexpr Vector(const Vector& v) : Base{v._val} {}
    constexpr Vector(Vector&& v) : Base{std::move(v._val)} {}
    constexpr Vector(const std::array<uint64_t, Base::M>& val)
        : _BitVector<N, Vector<bool, N>>{val} {}
    constexpr Vector(std::array<uint64_t, Base::M>&& val)
        : _BitVector<N, Vector<bool, N>>{std::move(val)} {}
};
}  // namespace TFHE