#include <cassert>
#include <random>

#include "tlwe.hpp"

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

// TODO: Google Testとか使う
// CMake Targetでテストできるようにしたい
int main() { test_TLWE(); }