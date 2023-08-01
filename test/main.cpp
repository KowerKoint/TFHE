#include <cassert>
#include <random>

#include "tlwe.hpp"
#include "trlwe.hpp"

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

// TODO: Google Testとか使う
// CMake Targetでテストできるようにしたい
int main() {
    test_TLWE();
    test_TRLWE();
}