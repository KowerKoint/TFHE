#pragma once

#include <array>
#include <bitset>
#include <tuple>

#include "random.hpp"
#include "types.hpp"

namespace TFHE {

struct TLWEParameter128BitSecurity {
    constexpr static int N = 636;
    constexpr static double ALPHA = 0.0000925119974676756;
};

template <typename TLWEParameter = TLWEParameter128BitSecurity>
class TLWE {
public:
    constexpr static int N = TLWEParameter::N;
    constexpr static double ALPHA = TLWEParameter::ALPHA;

private:
    Random random;
    Vector<TorusValue, N> generate_a() {
        Vector<TorusValue, N> a;
        std::generate(
            a.begin(), a.end(), [this] { return random.uniform_torus(); });
        return a;
    }
    TorusValue generate_e() { return TorusValue{random.normal() * ALPHA}; }

    template <int N>
    constexpr static Vector<TorusValue, N + 1> concat_ba(
        const TorusValue& b, const Vector<TorusValue, N>& a) {
        Vector<TorusValue, N + 1> ba;
        ba[0] = b;
        for (int i = 0; i < N; i++) ba[1 + i] = a[i];
        return ba;
    }

    template <int N>
    constexpr static std::pair<TorusValue, Vector<TorusValue, N>> decompose_ba(
        const Vector<TorusValue, N + 1>& ba) {
        TorusValue b = ba[0];
        Vector<TorusValue, N> a;
        for (int i = 0; i < N; i++) a[i] = ba[1 + i];
        return std::make_pair(b, a);
    }

public:
    TLWE() : random{} {}

    Vector<bool, N> generate_s() { return random.bit_vector<N>(); }

    Vector<TorusValue, N + 1> encrypt(
        const TorusValue& m, const Vector<bool, N>& s) {
        Vector<TorusValue, N> a = generate_a();
        TorusValue e = generate_e();
        TorusValue b = a.dot(s) + m + e;
        return concat_ba(b, a);
    }
    Vector<TorusValue, N + 1> encrypt_single_binary(
        bool m, const Vector<bool, N>& s) {
        return encrypt(TorusValue(m), s);
    }

    TorusValue decrypt(
        const Vector<TorusValue, N + 1>& ba, const Vector<bool, N>& s) {
        auto [b, a] = decompose_ba<N>(ba);
        return (b - a.dot(s));
    }
    bool decrypt_single_binary(
        const Vector<TorusValue, N + 1>& ba, const Vector<bool, N>& s) {
        return (bool)decrypt(ba, s);
    }
};
}  // namespace TFHE