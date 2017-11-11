//
// Created by Game on 2017/11/11.
//

#include <cmath>
#include <vector>
#include <string>

#include "compressQuat.h"
#include "gtest/gtest.h"


const float mPi = 3.1415926535897932384626433832795f;


TEST(compressQuat, InitQuat) {
    cqQuaternion q{};
    cq_identity(&q);
}


TEST(compressQuat, Euler2Quaternion) {
    size_t test_size = 720;

    for (size_t x = 0; x < test_size;) {
        float x_rad = x * (mPi / 180);
        x += 2;

        for (int y = 0; y < test_size;) {
            float y_rad = y * (mPi / 180);
            y += 2;

            for (int z = 0; z < test_size;) {
                float z_rad = z * (mPi / 180);
                z += 2;

                cqQuaternion test_quat1 {0,0,0,0};
                cqQuaternion test_quat2 {0,0,0,0};

                cqFloat3 euler1 { x_rad, y_rad, z_rad };
                cqFloat3 euler2 { 0, 0, 0 };

                cq_fromEuler(&euler1, &test_quat1);
                cq_toEuler(&test_quat1, &euler2);
                cq_fromEuler(&euler2, &test_quat2);

                auto check = cq_compare(&test_quat1, &test_quat2, 0.5f) == cqTRUE;
                const float diff_w = test_quat1.x * test_quat2.x +
                                     test_quat1.y * test_quat1.y +
                                     test_quat1.z * test_quat1.z +
                                     test_quat1.w * test_quat1.w;
                // Converts w back to an angle.
                const float angle = 2.f * acosf(fminf(fabsf(diff_w), 1.f));
                EXPECT_TRUE(check) << std::string("angle : ") + std::to_string(angle);
            }
        }
    }
}

TEST(compressQuat, QuaternionCompress) {
    size_t test_size = 720;

    for (size_t x = 0; x < test_size;) {
        float x_rad = x * (mPi / 180);
        x += 2;

        for (int y = 0; y < test_size;) {
            float y_rad = y * (mPi / 180);
            y += 2;

            for (int z = 0; z < test_size;) {
                float z_rad = z * (mPi / 180);
                z += 2;

                cqQuaternion test_quat {0,0,0,0};

                cqFloat3 euler { x_rad, y_rad, z_rad };
                cq_fromEuler(&euler, &test_quat);

                CompressQuat cquat {};
                compress_pack(&test_quat, &cquat);

                cqQuaternion unquat {};
                cqQuaternion unquatN {};
                cq_identity(&unquat);
                cq_identity(&unquatN);
                uncompress_pack(&cquat, &unquat);
                uncompress_packN(&cquat, &unquatN);

                auto check = cq_compare(&unquat, &test_quat, 2.0f) == cqTRUE &&
                        cq_compare(&unquatN, &test_quat, 2.5f) == cqTRUE;

                const float diff_w = test_quat.x * unquat.x +
                        test_quat.y * unquat.y +
                        test_quat.z * unquat.z +
                        test_quat.w * unquat.w;
                const float normalize_diff_w = test_quat.x * unquatN.x +
                                               test_quat.y * unquatN.y +
                                               test_quat.z * unquatN.z +
                                               test_quat.w * unquatN.w;

                // Converts w back to an angle.
                const float angle = 2.f * acosf(fminf(fabsf(diff_w), 1.f));
                const float normalize_angle = 2.f * acosf(fminf(fabsf(normalize_diff_w), 1.f));
                EXPECT_TRUE(check) << " Compress Quat : x= " +  std::to_string(test_quat.x) +
                                    " y= " + std::to_string(test_quat.y) +
                                    " z= " + std::to_string(test_quat.z) +
                                    " w= " + std::to_string(test_quat.w) +
                                    " Uncompress Quat : x= "+  std::to_string(unquat.x) +
                                    " y= " + std::to_string(unquat.y) +
                                    " z= " + std::to_string(unquat.z) +
                                    " w= " + std::to_string(unquat.w) +
                                    " angle= " + std::to_string(angle) +
                                    " angle(normalize)= " + std::to_string(normalize_angle);
            }
        }
    }
}