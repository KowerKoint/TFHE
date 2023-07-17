#pragma once

#include <array>
#include <bitset>
#include <tuple>

#include "random.hpp"
#include "types.hpp"

namespace TFHE {
class TLWE {
    constexpr static int N = 636;
    constexpr static double ALPHA = 0.0000925119974676756;

    Random random;

    TorusVector<N> generate_a() {
        TorusVector<N> a;
        for (auto& ai : a) {
            ai = random.uniform_torus();
        }
        return a;
    }
    TorusValue generate_e() { return random.normal() * ALPHA; }

public:
    TLWE() : random{} {}

    BitVector<N> generate_s() { return random.bit_vector<N>(); }

    TorusVector<N + 1> encrypt_single_binary(bool m, const BitVector<N>& s) {
        TorusVector<N> a = generate_a();
        TorusValue e = generate_e();
        TorusValue b = a.dot(s) + TorusValue(m) + e;
        return concat_ba(b, a);
    }

    bool decrypt_single_binary(
        const TorusVector<N + 1>& ba, const BitVector<N>& s) {
        TorusValue b;
        TorusVector<N> a;
        std::tie(b, a) = decompose_ba<N>(ba);
        return b - a.dot(s);
    }
};
}  // namespace TFHE