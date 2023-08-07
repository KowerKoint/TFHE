#pragma once

#include "trlwe.hpp"

namespace TFHE {
struct TRGSWParameterDefault {
    constexpr static int K = 2;
    constexpr static int BG_BIT = 8;
    constexpr static int L = 2;
};

template <typename TRGSWParameter = TRGSWParameterDefault,
    typename TRLWEParameter = TRLWEParameterDefault,
    std::enable_if_t<TRGSWParameter::K == TRLWEParameter::K>* = nullptr>
class TRGSW {
    constexpr static int BG_BIT = TRGSWParameter::BG_BIT;
    constexpr static int N = TRLWEParameter::N;
    constexpr static int K = TRGSWParameter::K;
    constexpr static int L = TRGSWParameter::L;
    using Int = SignedInt<BG_BIT>;

    TRLWE<TRLWEParameter>& trlwe;

    Vector<Polynomial<Int, N>, L> decomposition(
        const Polynomial<TorusValue, N>& a) {
        constexpr uint32_t round_offset =
            1 << (32 - L * BG_BIT - 1);  // bg^(-l)/2 四捨五入のために足す
        Vector<Polynomial<int, N>, L> a_hat;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < N; j++) {
                a_hat[i][j] = (a[j].get_raw_value() + round_offset) >>
                                  (32 - BG_BIT * (i + 1)) &
                              ((1 << BG_BIT) - 1);
            }
        }
        Vector<Polynomial<Int, N>, L> a_bar;
        for (int i = L - 1; i >= 0; i--) {
            for (int j = 0; j < N; j++) {
                if (a_hat[i][j] >= (1 << (BG_BIT - 1))) {
                    a_bar[i][j] = a_hat[i][j] - (1 << BG_BIT);
                    if (i) a_hat[i - 1][j]++;
                } else {
                    a_bar[i][j] = a_hat[i][j];
                }
            }
        }
        return a_bar;
    }

    Vector<Polynomial<TorusValue, N>, K + 1> encrypt_zero_by_trlwe(
        const Vector<Polynomial<bool, N>, K>& s) {
        return trlwe.encrypt_binary_polynomial({}, s);
    }

public:
    TRGSW(TRLWE<TRLWEParameter>& trlwe) : trlwe(trlwe) {}

    Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1>
    encrypt_integer_polynomial(
        const Polynomial<Int, N>& mu, const Vector<Polynomial<bool, N>, K>& s) {
        Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1> ret;
        for (int i = 0; i < K + 1; i++) {
            for (int j = 0; j < L; j++) {
                std::transform(mu.begin(), mu.end(), ret[i * L + j][i].begin(),
                    [&](Int val) {
                        // 値としてはval/2^(bgbit*(j+1))が正しいが、
                        // TorusValue型のuint32_t表現を直接出す
                        return TorusValue::from_raw_value(
                            (uint32_t)val << (32 - BG_BIT * (j + 1)));
                    });
            }
        }
        for (int i = 0; i < K + 1; i++) {
            for (int j = 0; j < L; j++) {
                ret[i * L + j] += encrypt_zero_by_trlwe(s);
            }
        }
        return ret;
    }

    Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1>
    encrypt_binary_polynomial(const Polynomial<bool, N>& mu,
        const Vector<Polynomial<bool, N>, K>& s) {
        Polynomial<Int, N> mu_integer;
        std::transform(mu.begin(), mu.end(), mu_integer.begin(),
            [](bool b) { return (Int)b; });
        return encrypt_integer_polynomial(mu_integer, s);
    }

    Vector<Polynomial<TorusValue, N>, K + 1> external_product(
        const Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1>& c,
        const Vector<Polynomial<TorusValue, N>, K + 1>& ba) {
        Vector<Polynomial<Int, N>, (K + 1) * L> decomposed_ba;
        for (int i = 0; i < K + 1; i++) {
            Vector<Polynomial<Int, N>, L> decomposed = decomposition(ba[i]);
            std::copy(decomposed.begin(), decomposed.end(),
                decomposed_ba.begin() + i * (K + 1));
        }
        return decomposed_ba * c;
    }

    Vector<Polynomial<TorusValue, N>, K + 1> cmux(
        const Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1>& c,
        const Vector<Polynomial<TorusValue, N>, K + 1>& ba_0,
        const Vector<Polynomial<TorusValue, N>, K + 1>& ba_1) {
        return external_product(c, ba_1 - ba_0) + ba_0;
    }
};
}  // namespace TFHE