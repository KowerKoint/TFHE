#include <cassert>
#include <random>

#include "tlwe.hpp"
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
            for (int j = 0; j < N; j++) assert(tlwe_s_polynomial[j] == 0);
            tlwe_s_polynomial[0] = tlwe_lv0_s[i];
            bk[i] = trgsw.encrypt_binary_polynomial(tlwe_s_polynomial, trlwe_s);
        }

        for (int i = 0; i < 10; i++) {
            bool input = mt() & 1;
            /*
            TFHE::TorusValue m{
                rd[input](mt)};  // true: [0.05,0.45), false: [0.55,0.95)
            */
            TFHE::TorusValue m(input);
            auto lv0 = tlwe_lv0.encrypt(m, tlwe_lv0_s);
            std::cout << "tlwelv0_key:\n";
            for (int j = 0; j < TLWE_N; j++)
                std::cout << tlwe_lv0_s[j] << (j == TLWE_N - 1 ? '\n' : ' ');
            auto lv1 = trgsw.gate_bootstrapping_tlwe_to_tlwe(lv0, bk);
            std::cout << (double)m << ' '
                      << tlwe_lv1.decrypt_single_binary(lv1, tlwe_lv1_s)
                      << std::endl;
        }
    }
}

// TODO: Google Testとか使う
// CMake Targetでテストできるようにしたい
int main() {
    test_TLWE();
    test_TRLWE();
    test_CMUX();
    test_Bootstrapping();
}