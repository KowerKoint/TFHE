#pragma once

#include "tlwe.hpp"

namespace TFHE {
struct TLWEKeySwitchingParameterDefault {
    constexpr static int T = 5;
    constexpr static int BASE_BIT = 2;
};

template <typename TLWELv0, typename TLWELv1,
    typename TLWEKeySwitchingParameter = TLWEKeySwitchingParameterDefault>
class TLWEKeySwitching {
public:
    constexpr static int N_LV0 = TLWELv0::N;
    constexpr static int N_LV1 = TLWELv1::N;
    constexpr static int T = TLWEKeySwitchingParameter::T;
    constexpr static int BASE_BIT = TLWEKeySwitchingParameter::BASE_BIT;

private:
    TLWELv0& tlwe_lv0;

public:
    TLWEKeySwitching(TLWELv0& tlwe_lv0) : tlwe_lv0(tlwe_lv0) {}

    Vector<Vector<Vector<Vector<TorusValue, N_LV0 + 1>, (1 << BASE_BIT)>, T>,
        N_LV1>
    make_ks(
        const Vector<bool, N_LV0>& s_lv0, const Vector<bool, N_LV1>& s_lv1) {
        Vector<
            Vector<Vector<Vector<TorusValue, N_LV0 + 1>, (1 << BASE_BIT)>, T>,
            N_LV1>
            ks;
        for (int i = 0; i < N_LV1; i++) {
            for (int m = 0; m < T; m++) {
                for (int o = 1; o < (1 << BASE_BIT); o++) {
                    TorusValue s_imo = TorusValue::from_raw_value(
                        s_lv1[i] ? (uint32_t)o << (32 - (m + 1) * BASE_BIT)
                                 : 0U);
                    ks[i][m][o] = tlwe_lv0.encrypt(s_imo, s_lv0);
                }
            }
        }
        return ks;
    }

    Vector<TorusValue, N_LV0 + 1> identity_key_switch(
        const Vector<TorusValue, N_LV1 + 1>& ba_lv1,
        const Vector<
            Vector<Vector<Vector<TorusValue, N_LV0 + 1>, (1 << BASE_BIT)>, T>,
            N_LV1>& ks) {
        constexpr uint32_t round_offset =
            1U << (32 - T * BASE_BIT - 1);  // base^(-t)/2 四捨五入のために足す
        Vector<TorusValue, N_LV0 + 1> ba_lv0;
        ba_lv0[0] = ba_lv1[0];
        for (int i = 0; i < N_LV1; i++) {
            for (int m = 0; m < T; m++) {
                int o = (ba_lv1[1 + i].get_raw_value() + round_offset) >>
                            (32 - (m + 1) * BASE_BIT) &
                        ((1 << BASE_BIT) - 1);
                if (o) ba_lv0 -= ks[i][m][o];
            }
        }
        return ba_lv0;
    }
};
}  // namespace TFHE