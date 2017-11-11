#ifndef COMPRESS_QUAT_LIBRARY_H
#define COMPRESS_QUAT_LIBRARY_H
#include <stdint.h>

#define cqBOOL int
#define cqTRUE 1
#define cqFALSE 0
#define cqNullptr 0

#ifdef __cplusplus
extern "C" {
#endif

struct cqFloat3 {
    float x, y, z;
};


struct cqQuaternion {
    float x, y, z, w;
};

// Defines the compress quaternion, precision : 1/16384
// Rotation value is a quaternion. Quaternion are normalized, which means each
// component is in range [0:1]. This property allows to quantize the 3
// components to 3 signed integer 16 bits values. The 4th component is restored
// at runtime, using the knowledge that |w| = sqrt(1 - (a^2 + b^2 + c^2)).
// The sign of this 4th component is stored using 1 bit taken from the track
// member.
//
// In more details, compression algorithm stores the 3 smallest components of
// the quaternion and restores the largest. The 3 smallest can be pre-multiplied
// by sqrt(2) to gain some precision indeed.
struct CompressQuat {
    uint16_t largest : 2;  // The largest component of the quaternion.
    uint16_t sign    : 1;  // The sign of the largest component. 1 for negative.
    uint16_t sign_a  : 1;
    uint16_t a       : 14;
    uint16_t sign_b  : 1;
    uint16_t b       : 14;
    uint16_t sign_c  : 1;
    uint16_t c       : 14;
};

void cq_identity(struct cqQuaternion *_q);
void cq_mult(const struct cqQuaternion *_a, const struct cqQuaternion *_b, struct cqQuaternion *out);
void cq_conjugate(const struct cqQuaternion *_a, struct cqQuaternion *out);
void cq_normalize(struct cqQuaternion *_q);
cqBOOL cq_isNormalize(const struct cqQuaternion *_q);
cqBOOL cq_compare(const struct cqQuaternion *_a, const struct cqQuaternion *_b, float _tolerance);
void cq_fromEuler(const struct cqFloat3 *_euler, struct cqQuaternion *out);
void cq_toEuler(const struct cqQuaternion *q, struct cqFloat3 *_euler);

// compresses quaternion
// The 3 smallest components of the quaternion are quantized to 16 bits
// integers, while the largest is recomputed thanks to quaternion normalization
// property (x^2+y^2+z^2+w^2 = 1). Because the 3 components are the 3 smallest,
// their value cannot be greater than sqrt(2)/2. Thus quantization quality is
// improved by pre-multiplying each componenent by sqrt(2).
void compress_pack(const struct cqQuaternion *_src, struct CompressQuat *out);
void uncompress_pack(const struct CompressQuat *_src, struct cqQuaternion *out);

#ifdef __cplusplus
}
#endif

#endif