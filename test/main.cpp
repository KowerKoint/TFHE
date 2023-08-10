#include <cassert>
#include <chrono>
#include <iostream>
#include <random>

#include "hom_nand.hpp"
#include "tlwe.hpp"
#include "tlwe_key_switching.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"
#include "types.hpp"

void test_TLWE() {
    TFHE::TLWE<> tlwe;
    std::mt19937 mt(0);
    for (int i = 0; i < 100; i++) {
        auto s = tlwe.generate_s();
        bool msg = mt() & 1;
        auto c = tlwe.encrypt_single_binary(msg, s);
        assert(tlwe.decrypt_single_binary(c, s) == msg);
    }
}

void test_TRLWE() {
    TFHE::TRLWE<> trlwe;
    constexpr int N = TFHE::TRLWE<>::N;
    std::mt19937 mt(0);
    for (int i = 0; i < 100; i++) {
        auto s = trlwe.generate_s();
        TFHE::Polynomial<bool, N> msg;
        std::generate(msg.begin(), msg.end(), [&mt] { return mt() & 1; });
        auto c = trlwe.encrypt_binary_polynomial(msg, s);
        assert(trlwe.decrypt_binary_polynomial(c, s) == msg);
    }
}

void test_ExternalProduct() {
    TFHE::TRLWE<> trlwe;
    constexpr int N = TFHE::TRLWE<>::N;
    constexpr int K = TFHE::TRLWE<>::K;
    TFHE::TRGSW<> trgsw(trlwe);
    std::mt19937 mt(0);
    TFHE::Polynomial<bool, N> msg;
    for (int i = 0; i < N; i++) msg[i] = true;
    auto s = trlwe.generate_s();
    auto ba = trlwe.encrypt_binary_polynomial(msg, s);
    TFHE::Polynomial<int8_t, 512> one;
    one[0] = 1;
    auto c = trgsw.encrypt_integer_polynomial(one, s);
    auto e = trgsw.external_product(c, ba);
    for (int i = 0; i < K + 1; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << ' ' << (double)ba[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << "=====\n";
    for (int i = 0; i < K + 1; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << ' ' << (double)e[i][j];
        }
        std::cout << std::endl;
    }
}

void test_CMUX() {
    TFHE::TRLWE<> trlwe;
    constexpr int N = TFHE::TRLWE<>::N;
    constexpr int K = TFHE::TRLWE<>::K;
    TFHE::TRGSW<> trgsw(trlwe);
    std::mt19937 mt(0);
    for (int i = 0; i < 100; i++) {
        auto s = trlwe.generate_s();
        TFHE::Polynomial<bool, N> msgs[2];
        TFHE::Vector<TFHE::Polynomial<TFHE::TorusValue, N>, K + 1>
            msgs_encrypted[2];
        for (int j = 0; j < 2; j++) {
            std::generate(
                msgs[j].begin(), msgs[j].end(), [&mt] { return mt() & 1; });
            msgs_encrypted[j] = trlwe.encrypt_binary_polynomial(msgs[j], s);
        }
        int answer = mt() & 1;
        TFHE::Polynomial<bool, N> answer_poly;
        answer_poly[0] = answer;
        auto c = trgsw.encrypt_binary_polynomial(answer_poly, s);
        auto selected_msg = trgsw.cmux(c, msgs_encrypted[0], msgs_encrypted[1]);
        assert(
            msgs[answer] == trlwe.decrypt_binary_polynomial(selected_msg, s));
    }
}

void test_PolynomialXExp() {
    std::mt19937 mt(0);
    for (int _ = 0; _ < 100; _++) {
        TFHE::Polynomial<int, 512> poly;
        for (int i = 0; i < 512; i++) poly[i] = mt() % 10000;
        int n = mt() % 1024;
        TFHE::Polynomial<int, 512> x_n;
        if (n >= 512) {
            x_n[n - 512] = -1;
            n -= 1024;
        } else {
            x_n[n] = 1;
        }
        auto lhs = poly * x_n;
        auto rhs = poly.multiply_x_exp(n);
        assert(lhs == rhs);
    }
}

void test_Bootstrapping() {
    TFHE::TLWE<> tlwe_lv0;
    TFHE::TLWE<TFHE::TLWELv1ParameterDefault> tlwe_lv1;
    TFHE::TRLWE<> trlwe;
    constexpr int TLWE_N = TFHE::TLWE<>::N;
    constexpr int N = TFHE::TRLWE<>::N;
    constexpr int K = TFHE::TRLWE<>::K;
    constexpr int L = TFHE::TRGSW<>::L;
    TFHE::TRGSW<> trgsw(trlwe);
    std::mt19937 mt(0);
    std::uniform_real_distribution<> rd[2] = {
        std::uniform_real_distribution<>{0.55, 0.95},
        std::uniform_real_distribution<>{0.05, 0.45}};
    for (int _ = 0; _ < 10; _++) {
        // secrets
        auto tlwe_lv0_s = tlwe_lv0.generate_s();
        auto trlwe_s = trlwe.generate_s();
        auto tlwe_lv1_s = trlwe.extract_tlwe_lv0_key(trlwe_s);
        TFHE::Vector<TFHE::Matrix<TFHE::Polynomial<TFHE::TorusValue, N>,
                         (K + 1) * L, K + 1>,
            TLWE_N>
            bk;
        for (int i = 0; i < TLWE_N; i++) {
            TFHE::Polynomial<bool, N> tlwe_s_polynomial;
            tlwe_s_polynomial[0] = tlwe_lv0_s[i];
            bk[i] = trgsw.encrypt_binary_polynomial(tlwe_s_polynomial, trlwe_s);
        }

        for (int i = 0; i < 10; i++) {
            bool input = mt() & 1;
            TFHE::TorusValue m{
                rd[input](mt)};  // true: [0.05,0.45), false: [0.55,0.95)
            auto lv0 = tlwe_lv0.encrypt(m, tlwe_lv0_s);
            auto lv1 = trgsw.gate_bootstrapping_tlwe_to_tlwe(lv0, bk);
            assert(input == tlwe_lv1.decrypt_single_binary(lv1, tlwe_lv1_s));
        }
    }
}

void test_KeySwitching() {
    TFHE::TLWE<> tlwe_lv0;
    TFHE::TLWE<TFHE::TLWELv1ParameterDefault> tlwe_lv1;
    constexpr int N_LV0 = TFHE::TLWEParameter128BitSecurity::N;
    constexpr int N_LV1 = TFHE::TLWELv1ParameterDefault::N;
    TFHE::TLWEKeySwitching<decltype(tlwe_lv0), decltype(tlwe_lv1)>
        key_switching(tlwe_lv0);
    std::mt19937 mt(0);
    std::uniform_real_distribution<> rd[2] = {
        std::uniform_real_distribution<>{0.55, 0.95},
        std::uniform_real_distribution<>{0.05, 0.45}};
    for (int _ = 0; _ < 3; _++) {
        // secrets
        auto tlwe_lv0_s = tlwe_lv0.generate_s();
        auto tlwe_lv1_s = tlwe_lv1.generate_s();
        auto ks = key_switching.make_ks(tlwe_lv0_s, tlwe_lv1_s);
        for (int i = 0; i < 10; i++) {
            std::cout << "i=" << i << std::endl;
            bool input = mt() & 1;
            TFHE::TorusValue m{
                rd[input](mt)};  // true: [0.05,0.45), false: [0.55,0.95)
            auto lv1 = tlwe_lv1.encrypt(m, tlwe_lv1_s);
            auto lv1_lv0 = key_switching.identity_key_switch(lv1, ks);
            auto d = tlwe_lv0.decrypt(lv1_lv0, tlwe_lv0_s);
            std::cout << (double)m << ' ' << (double)d << std::endl;
            assert(input == (bool)d);
        }
    }
}

void test_HomNAND() {
    TFHE::HomNAND hom_nand;
    std::mt19937 mt;
    for (int _ = 0; _ < 3; _++) {
        auto secret_key = hom_nand.generate_secret_key();
        auto evaluate_key = hom_nand.make_evaluate_key(secret_key);
        bool a[10];
        TFHE::HomNAND<>::Cipher ca[10];
        for (int i = 0; i < 10; i++) {
            a[i] = mt() & 1;
            ca[i] = hom_nand.encrypt(a[i], secret_key);
        }
        int x[50], y[50], z[50];
        for (int i = 0; i < 50; i++) {
            x[i] = mt() % 10;
            y[i] = mt() % 9;
            if (y[i] == x[i]) y[i] = 9;
            z[i] = mt() % 8;
            if (z[i] == x[i]) z[i] = 8;
            if (z[i] == y[i]) z[i] = 9;
        }
        auto st = std::chrono::system_clock::now();
        for (int i = 0; i < 50; i++) {
            // a[z] <- a[x] nand a[y]
            ca[z[i]] = hom_nand.nand(ca[x[i]], ca[y[i]], evaluate_key);
        }
        auto ed = std::chrono::system_clock::now();
        std::cout << "#" << _ << ": "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(
                         ed - st)
                             .count() /
                         50.
                  << "ms / nand" << std::endl;
        for (int i = 0; i < 50; i++) {
            a[z[i]] = !(a[x[i]] && a[y[i]]);
        }
        for (int i = 0; i < 10; i++) {
            std::cout << ' ' << a[i];
        }
        std::cout << std::endl;
        for (int i = 0; i < 10; i++) {
            assert(a[i] == hom_nand.decrypt(ca[i], secret_key));
        }
    }
}

void test_FFT() {
    std::mt19937 mt(0);
    TFHE::Polynomial<std::complex<double>, 512> a;
    for (int i = 0; i < 512; i++) {
        a[i].real(mt());
        a[i].imag(mt());
    }
    std::complex<double> omega = std::polar(1., -M_PI * 2 / 512);
    std::array<std::complex<double>, 512> omega_pow;
    omega_pow[0] = 1.;
    for (int i = 0; i < 511; i++) omega_pow[i + 1] = omega_pow[i] * omega;
    auto a_fft = fft(a, omega_pow);
    std::reverse(omega_pow.begin() + 1, omega_pow.end());
    auto aa = fft(a_fft, omega_pow);
    for (int i = 0; i < 512; i++) {
        aa[i] /= 512;
        std::cout << a[i] << ' ' << aa[i] << std::endl;
    }
}

// TODO: Google Testとか使う
// CMake Targetでテストできるようにしたい
int main() {
    /* test_TLWE(); */
    /* test_TRLWE(); */
    /* test_ExternalProduct(); */
    /* test_CMUX(); */
    /* test_PolynomialXExp(); */
    /* test_Bootstrapping(); */
    /* test_KeySwitching(); */
    test_HomNAND();
    /* test_FFT(); */
}