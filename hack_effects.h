#ifndef HACK_EFFECTS_H
#define HACK_EFFECTS_H

#include <stdint.h>
#define MAX_NUM_SAMPLES_OF_NOISE   (240000) /* 5 seconds at 48KHz, 8-bit PCM */
#define PHASE_ADJUST_IN_SAMPLES    (12)     /* 250 useconds */
#define FXP_Q31xQ31_MPY(_a, _b, scale) ((int32_t)((((int64_t)_a * (int64_t)_b) + ((int64_t)1 << (scale - 1))) >> scale))
#define GAIN_SCALAR                    (15)
#define MAX_SIGNAL_AMP                 (2147483647)
#define MIN_SIGNAL_AMP                 (-2147483647 - 1)
#define CRACKLE_AMP_DECAY_FACTOR       (15)
#define RANDOMNESS_SCALAR              (15)

/* Empirically determined */
#define RANDOMNESS_SEED                (322122547)
#define RANDOMNESS_ALPHA               (104729)
#define RANDOMNESS_INCREMENT_C         (41919)
#define RANDOMNESS_APPLICATION_HI_LIM  (134217728)
#define RANDOMNESS_APPLICATION_LO_LIM  (16777216)
#define RANDOMNESS_CADENCE_IN_SAMPLES  (4096)


void hack_ambient_noise_effect(int *inpBufPtr, int *opBufPtr, unsigned int numSamples);
void hack_vinyl_crackle_effect(int *inpBufPtr, int *opBufPtr, unsigned int numSamples);


#endif // !HACK_EFFECTS_H
