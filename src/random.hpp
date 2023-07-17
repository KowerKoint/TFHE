#pragma once

#include <random>

#include "types.hpp"

namespace TFHE {
class Random {
    std::random_device rd;
    std::normal_distribution<> normal_dist;
    std::uniform_real_distribution<> uniform_real_dist;
    std::uniform_int_distribution<uint32_t> uniform_int_dist_32;
    std::uniform_int_distribution<uint64_t> uniform_int_dist_64;

public:
    Random()
        : rd{},
          normal_dist{},
          uniform_real_dist{},
          uniform_int_dist_32{},
          uniform_int_dist_64{} {}

    double normal() { return normal_dist(rd); }

    double uniform_real() { return uniform_real_dist(rd); }

    TorusValue uniform_torus() { return uniform_int_dist_32(rd); }

    template <int N>
    BitVector<N> bit_vector() {
        constexpr int M = (N + 63) / 64;
        std::array<uint64_t, M> ret_array;
        for (int i = 0; i < M; i++) {
            ret_array[i] = uniform_int_dist_64(rd);
        }
        return BitVector<N>(std::move(ret_array));
    }
};
}  // namespace TFHE