/*
 * Copyright (c) 2020 Muhammad Faiz <mfcc64@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/* Audio visualization based on showcqt mpv/ffmpeg audio visualization */
/* See https://github.com/FFmpeg/FFmpeg/blob/master/libavfilter/avf_showcqt.c */

// Modified by Eknous for unscope - https://github.com/Eknous-P/unscope

#ifndef SHOWCQT_H_INCLUDED
#define SHOWCQT_H_INCLUDED 1

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdint.h>

#define DECLARE_ALIGNED(n) __attribute__((__aligned__(n)))
#define ALWAYS_INLINE __inline__ __attribute__((__always_inline__))

#define CQT_MAX_FFT_SIZE 32768
#define CQT_MAX_KERNEL_SIZE (6*256*1024)
#define CQT_MAX_WIDTH 7680

// uncomment this define to use builtin color array rendering
// #define SHOWCQT_RENDERER

#ifdef SHOWCQT_RENDERER
#define CQT_MAX_HEIGHT 4320
#define CQT_MIN_VOL 1.0f
#define CQT_MAX_VOL 100.0f
#endif


typedef struct Complex {
    float re, im;
} Complex;

#ifdef SHOWCQT_RENDERER
typedef struct ColorF {
    float r, g, b, h;
} ColorF;
#endif

typedef union Kernel {
    float f;
    int   i;
} Kernel;

typedef struct KernelIndex {
    int len;
    int start;
} KernelIndex;

typedef struct ShowCQT {
    /* args */
    float        input[2][CQT_MAX_FFT_SIZE];
#ifdef SHOWCQT_RENDERER
    unsigned int output[CQT_MAX_WIDTH];
#endif
    float        magOutput[CQT_MAX_WIDTH];

    /* tables */
    Complex     exp_tbl[CQT_MAX_FFT_SIZE+CQT_MAX_FFT_SIZE/4];
    int16_t     perm_tbl[CQT_MAX_FFT_SIZE/4];
    float       attack_tbl[CQT_MAX_FFT_SIZE/8];
#ifndef SHOWCQT_RENDERER
    uint8_t     padding[1024];
#endif

    /* buffers */
    Complex     fft_buf[CQT_MAX_FFT_SIZE+128];
#ifdef SHOWCQT_RENDERER
    ColorF      color_buf[CQT_MAX_WIDTH*2];
#endif
    float       rcp_h_buf[CQT_MAX_WIDTH];

    /* kernel */
    KernelIndex kernel_index[CQT_MAX_WIDTH*2];
    float       kernel[CQT_MAX_KERNEL_SIZE];

    /* props */
#ifdef SHOWCQT_RENDERER
    int         height;
#endif
    int         width;
    int         aligned_width;
    int         fft_size;
    int         t_size;
    int         attack_size;
#ifdef SHOWCQT_RENDERER
    float       sono_v;
    float       bar_v;
    int         prerender;
#endif
} ShowCQT;

float* CQT_get_input_array(ShowCQT* cqt, int index);
#ifdef SHOWCQT_RENDERER
unsigned int *CQT_get_output_array(ShowCQT* cqt);
ColorF *CQT_get_color_array(ShowCQT* cqt);
#endif
int CQT_init(ShowCQT* cqt, int rate, int width, int height, float bar_v, float sono_v, int super);
void CQT_calc(ShowCQT* cqt);
#ifdef SHOWCQT_RENDERER
void CQT_render_line_alpha(ShowCQT* cqt, int y, uint8_t alpha);
void CQT_render_line_opaque(ShowCQT* cqt, int y);
#endif

#ifdef __cplusplus
}
#endif

#endif
