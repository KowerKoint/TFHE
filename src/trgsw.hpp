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
public:
    constexpr static int BG_BIT = TRGSWParameter::BG_BIT;
    constexpr static int N = TRLWEParameter::N;
    constexpr static int K = TRGSWParameter::K;
    constexpr static int L = TRGSWParameter::L;
    using Int = SignedInt<BG_BIT>;

private:
    TRLWE<TRLWEParameter>& trlwe;

    constexpr static Vector<Polynomial<TorusValue, N>, K + 1>
    gen_test_vector() {
        Vector<Polynomial<TorusValue, N>, K + 1> ret;
        for (int i = 0; i < N; i++) {
            ret[0][i] = TorusValue(true);
        }
        return ret;
    }

    Vector<Polynomial<Int, N>, L> decomposition(
        const Polynomial<TorusValue, N>& a) {
        constexpr uint32_t round_offset =
            1U << (32 - L * BG_BIT - 1);  // bg^(-l)/2 四捨五入のために足す
        Vector<Polynomial<int, N>, L> a_hat;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < N; j++) {
                a_hat[i][j] = (a[j].get_raw_value() + round_offset) >>
                                  (32 - BG_BIT * (i + 1)) &
                              ((1U << BG_BIT) - 1);
            }
        }
        Vector<Polynomial<Int, N>, L> a_bar;
        for (int i = L - 1; i >= 0; i--) {
            for (int j = 0; j < N; j++) {
                if (a_hat[i][j] >= (1U << (BG_BIT - 1))) {
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
        constexpr Polynomial<TorusValue, N> zero_polynomial;
        return trlwe.encrypt(zero_polynomial, s);
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
                decomposed_ba.begin() + i * L);
        }
        return decomposed_ba * c;
    }

    Vector<Polynomial<TorusValue, N>, K + 1> cmux(
        const Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1>& c,
        const Vector<Polynomial<TorusValue, N>, K + 1>& ba_0,
        const Vector<Polynomial<TorusValue, N>, K + 1>& ba_1) {
        return external_product(c, ba_1 - ba_0) + ba_0;
    }

    template <int TLWE_N>
    Vector<Polynomial<TorusValue, N>, K + 1> blind_rotate(
        const Vector<TorusValue, TLWE_N + 1>& tlwe_lv0,
        const Vector<Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1>,
            TLWE_N>& bk,
        const Vector<Polynomial<TorusValue, N>, K + 1>& test_vector) {
        auto get_log_2 = [](int n) {
            int n_bit = 0;
            while ((1 << n_bit) < n) n_bit++;
            return n_bit;
        };
        constexpr int N_bit = get_log_2(N);
        static_assert((1 << N_bit) == N);  // N = 2^N_bit

        auto trlwe_multiply_x_exp =
            [](const Vector<Polynomial<TorusValue, N>, K + 1>& trlwe, int n) {
                // trlweにX^nを乗算
                Vector<Polynomial<TorusValue, N>, K + 1> ret;
                std::transform(trlwe.begin(), trlwe.end(), ret.begin(),
                    [&](const Polynomial<TorusValue, N>& p) {
                        return p.multiply_x_exp(n);
                    });
                return ret;
            };
        int b_2n = tlwe_lv0[0].get_raw_value() >>
                   (32 - (N_bit + 1));  // floor(b * 2N);
        Vector<Polynomial<TorusValue, N>, K + 1> ret =
            trlwe_multiply_x_exp(test_vector, -b_2n);
        for (int i = 0; i < TLWE_N; i++) {
            int a_2n = (tlwe_lv0[1 + i].get_raw_value() +
                           (1U << (32 - (N_bit + 1) - 1))) >>
                       (32 - (N_bit + 1));  // round(a_i * 2N)
            Vector<Polynomial<TorusValue, N>, K + 1> rotated =
                trlwe_multiply_x_exp(ret, a_2n);
            Vector<Polynomial<TorusValue, N>, K + 1> selected =
                cmux(bk[i], ret, rotated);  // s==1のときret←ret*x^(a_2n)
            ret = selected;
        }
        return ret;
    }

    template <int TLWE_N>
    Vector<TorusValue, N * K + 1> gate_bootstrapping_tlwe_to_tlwe(
        const Vector<TorusValue, TLWE_N + 1>& tlwe_lv0,
        const Vector<Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1>,
            TLWE_N>& bk) {
        constexpr Vector<Polynomial<TorusValue, N>, K + 1> test_vector =
            gen_test_vector();
        return trlwe.sample_extract_index(
            blind_rotate(tlwe_lv0, bk, test_vector), 0);
    }
};
}  // namespace TFHE