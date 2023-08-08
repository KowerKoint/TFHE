#pragma once

#include <algorithm>
#include <random>

#include "types.hpp"

namespace TFHE {
class Random {
    std::random_device rd;
    std::normal_distribution<> normal_dist;
    std::uniform_real_distribution<> uniform_real_dist;
    std::uniform_int_distribution<uint32_t> uniform_int_dist_32;
    std::uniform_int_distribution<uint64_t> uniform_int_dist_64;

    template <int N>
    std::array<bool, N> _bit_array() {
        constexpr int M = (N + 63) / 64;
        std::array<uint64_t, M> ret_bits;
        std::generate(ret_bits.begin(), ret_bits.end(),
            [this] { return uniform_int_dist_64(rd); });
        std::array<bool, N> ret;
        for (int i = 0; i < N; i++) ret[i] = ret_bits[i >> 6] >> (i & 0x3F) & 1;
        return ret;
    }

public:
    Random()
        : rd{},
          normal_dist{},
          uniform_real_dist{},
          uniform_int_dist_32{},
          uniform_int_dist_64{} {}

    double normal() { return normal_dist(rd); }

    double uniform_real() { return uniform_real_dist(rd); }

    TorusValue uniform_torus() {
        return TorusValue::from_raw_value(uniform_int_dist_32(rd));
    }

    template <int N>
    Polynomial<bool, N> bit_polynomial() {
        return Polynomial<bool, N>{_bit_array<N>()};
    }

    template <int N>
    Vector<bool, N> bit_vector() {
        return Vector<bool, N>{_bit_array<N>()};
    }
};
}  // namespace TFHE