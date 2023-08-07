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

// TODO: Google Testとか使う
// CMake Targetでテストできるようにしたい
int main() {
    test_TLWE();
    test_TRLWE();
    test_CMUX();
}