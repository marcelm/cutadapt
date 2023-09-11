#include <stdint.h>
#ifdef __SSE2__
#include "emmintrin.h"
#endif

static const double SCORE_TO_ERROR_RATE[94] = {
    1.0L,                     // 0
    0.7943282347242815L,      // 1
    0.6309573444801932L,      // 2
    0.5011872336272722L,      // 3
    0.3981071705534972L,      // 4
    0.31622776601683794L,     // 5
    0.251188643150958L,       // 6
    0.19952623149688797L,     // 7
    0.15848931924611134L,     // 8
    0.12589254117941673L,     // 9
    0.1L,                     // 10
    0.07943282347242814L,     // 11
    0.06309573444801933L,     // 12
    0.05011872336272722L,     // 13
    0.039810717055349734L,    // 14
    0.03162277660168379L,     // 15
    0.025118864315095794L,    // 16
    0.0199526231496888L,      // 17
    0.015848931924611134L,    // 18
    0.012589254117941675L,    // 19
    0.01L,                    // 20
    0.007943282347242814L,    // 21
    0.00630957344480193L,     // 22
    0.005011872336272725L,    // 23
    0.003981071705534973L,    // 24
    0.0031622776601683794L,   // 25
    0.0025118864315095794L,   // 26
    0.001995262314968879L,    // 27
    0.001584893192461114L,    // 28
    0.0012589254117941675L,   // 29
    0.001L,                   // 30
    0.0007943282347242813L,   // 31
    0.000630957344480193L,    // 32
    0.0005011872336272725L,   // 33
    0.00039810717055349735L,  // 34
    0.00031622776601683794L,  // 35
    0.00025118864315095795L,  // 36
    0.00019952623149688788L,  // 37
    0.00015848931924611142L,  // 38
    0.00012589254117941674L,  // 39
    0.0001L,                  // 40
    7.943282347242822E-05L,   // 41
    6.309573444801929E-05L,   // 42
    5.011872336272725E-05L,   // 43
    3.9810717055349695E-05L,  // 44
    3.1622776601683795E-05L,  // 45
    2.5118864315095822E-05L,  // 46
    1.9952623149688786E-05L,  // 47
    1.584893192461114E-05L,   // 48
    1.2589254117941661E-05L,  // 49
    1E-05L,                   // 50
    7.943282347242822E-06L,   // 51
    6.30957344480193E-06L,    // 52
    5.011872336272725E-06L,   // 53
    3.981071705534969E-06L,   // 54
    3.162277660168379E-06L,   // 55
    2.5118864315095823E-06L,  // 56
    1.9952623149688787E-06L,  // 57
    1.584893192461114E-06L,   // 58
    1.2589254117941661E-06L,  // 59
    1E-06L,                   // 60
    7.943282347242822E-07L,   // 61
    6.30957344480193E-07L,    // 62
    5.011872336272725E-07L,   // 63
    3.981071705534969E-07L,   // 64
    3.162277660168379E-07L,   // 65
    2.5118864315095823E-07L,  // 66
    1.9952623149688787E-07L,  // 67
    1.584893192461114E-07L,   // 68
    1.2589254117941662E-07L,  // 69
    1E-07L,                   // 70
    7.943282347242822E-08L,   // 71
    6.30957344480193E-08L,    // 72
    5.011872336272725E-08L,   // 73
    3.981071705534969E-08L,   // 74
    3.162277660168379E-08L,   // 75
    2.511886431509582E-08L,   // 76
    1.9952623149688786E-08L,  // 77
    1.5848931924611143E-08L,  // 78
    1.2589254117941661E-08L,  // 79
    1E-08L,                   // 80
    7.943282347242822E-09L,   // 81
    6.309573444801943E-09L,   // 82
    5.011872336272715E-09L,   // 83
    3.981071705534969E-09L,   // 84
    3.1622776601683795E-09L,  // 85
    2.511886431509582E-09L,   // 86
    1.9952623149688828E-09L,  // 87
    1.584893192461111E-09L,   // 88
    1.2589254117941663E-09L,  // 89
    1E-09L,                   // 90
    7.943282347242822E-10L,   // 91
    6.309573444801942E-10L,   // 92
    5.011872336272714E-10L,   // 93
};

static inline double
expected_errors_from_phreds(const uint8_t *phreds, size_t phreds_length, uint8_t base) {
    const uint8_t *end_ptr = phreds + phreds_length;
    const uint8_t *cursor = phreds;
    double expected_errors = 0.0;
    uint8_t max_phred = 126 - base;
    #ifdef __SSE2__ 
    const uint8_t *vec_end_ptr = end_ptr - sizeof(__m128i);
    __m128d accumulator = _mm_set1_pd(0.0);
    while (cursor < vec_end_ptr) {
        __m128i phred_array = _mm_loadu_si128((__m128i *)cursor);
        __m128i illegal_phreds = _mm_cmpgt_epi8(phred_array, _mm_set1_epi8(126));
        illegal_phreds = _mm_or_si128(
            illegal_phreds, _mm_cmplt_epi8(phred_array, _mm_set1_epi8(base)));
        if (_mm_movemask_epi8(illegal_phreds)) {
            return -1.0;
        }
        /* By explicitly setting multiple accumulators, the processor 
           can perform out of order execution for increased speed 
            See also: https://stackoverflow.com/a/36591776/16437839   
        */
        __m128d accumulator1 = _mm_add_pd(
            _mm_set_pd(
                SCORE_TO_ERROR_RATE[cursor[0] - base],
                SCORE_TO_ERROR_RATE[cursor[1] - base]
            ),
            _mm_set_pd(
                SCORE_TO_ERROR_RATE[cursor[2] - base],
                SCORE_TO_ERROR_RATE[cursor[3] - base]
            )
        );
        __m128d accumulator2 = _mm_add_pd(
            _mm_set_pd(
                SCORE_TO_ERROR_RATE[cursor[4] - base],
                SCORE_TO_ERROR_RATE[cursor[5] - base]
            ),
            _mm_set_pd(
                SCORE_TO_ERROR_RATE[cursor[6] - base],
                SCORE_TO_ERROR_RATE[cursor[7] - base]
            )
        );
        __m128d accumulator3 = _mm_add_pd(
            _mm_set_pd(
                SCORE_TO_ERROR_RATE[cursor[8] - base],
                SCORE_TO_ERROR_RATE[cursor[9] - base]
            ),
            _mm_set_pd(
                SCORE_TO_ERROR_RATE[cursor[10] - base],
                SCORE_TO_ERROR_RATE[cursor[11] - base]
            )
        );
        __m128d accumulator4 = _mm_add_pd(
            _mm_set_pd(
                SCORE_TO_ERROR_RATE[cursor[12] - base],
                SCORE_TO_ERROR_RATE[cursor[13] - base]
            ),
            _mm_set_pd(
                SCORE_TO_ERROR_RATE[cursor[14] - base],
                SCORE_TO_ERROR_RATE[cursor[15] - base]
            )
        ); 
        accumulator = _mm_add_pd(accumulator, accumulator1);
        accumulator = _mm_add_pd(accumulator, accumulator2);
        accumulator = _mm_add_pd(accumulator, accumulator3);
        accumulator = _mm_add_pd(accumulator, accumulator4);
        cursor += sizeof(__m128i);
    }
    double double_store[2];
    _mm_store_pd(double_store, accumulator);
    expected_errors = double_store[0] + double_store[1];
    #endif
    while (cursor < end_ptr) {
        uint8_t phred = *cursor - base;
        if (phred > max_phred) {
            return -1.0;
        }
        expected_errors += SCORE_TO_ERROR_RATE[phred];
        cursor += 1;
    }
    return expected_errors;
}
