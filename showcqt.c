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

#include "showcqt.h"

float *CQT_get_input_array(ShowCQT* cqt, int index) {
    return cqt->input[!!index];
}

#ifdef SHOWCQT_RENDERER
unsigned int *CQT_get_output_array(ShowCQT* cqt) {
    return cqt->output;
}

ColorF *CQT_get_color_array(ShowCQT* cqt) {
    return cqt->color_buf;
}
#endif

unsigned int CQT_revbin(unsigned x, int bits) {
    unsigned int m = 0x55555555;
    x = ((x & m) << 1) | ((x & ~m) >> 1);
    m = 0x33333333;
    x = ((x & m) << 2) | ((x & ~m) >> 2);
    m = 0x0F0F0F0F;
    x = ((x & m) << 4) | ((x & ~m) >> 4);
    m = 0x00FF00FF;
    x = ((x & m) << 8) | ((x & ~m) >> 8);
    m = 0x0000FFFF;
    x = ((x & m) << 16) | ((x & ~m) >> 16);
    return (x >> (32 - bits)) & ((1 << bits) - 1);
}

void CQT_gen_perm_tbl(ShowCQT* cqt, int bits) {
    int n = 1 << bits;
    for (int x = 0; x < n; x++)
        cqt->perm_tbl[x] = CQT_revbin(x, bits);
}

void CQT_gen_exp_tbl(ShowCQT* cqt) {
    int n=cqt->fft_size;
    double mul;
    for (int k = 2; k < n; k *= 2) {
        mul = 2.0 * M_PI / k;
        for (int x = 0; x < k/2; x++)
            cqt->exp_tbl[k+x] = (Complex){ cosf(mul*x), -sinf(mul*x) };
        mul = 3.0 * M_PI / k;
        for (int x = 0; x < k/2; x++)
            cqt->exp_tbl[k+k/2+x] = (Complex){ cosf(mul*x), -sinf(mul*x) };
    }
    mul = 2.0 * M_PI / n;
    for (int x = 0; x < n/4; x++)
        cqt->exp_tbl[n+x] = (Complex){ cosf(mul*x), -sinf(mul*x) };
}

#define C_ADD(a, b) (Complex){ (a).re + (b).re, (a).im + (b).im }
#define C_SUB(a, b) (Complex){ (a).re - (b).re, (a).im - (b).im }
#define C_MUL(a, b) (Complex){ (a).re * (b).re - (a).im * (b).im, (a).re * (b).im + (a).im * (b).re }
#define C_AIM(a, b) (Complex){ (a).re - (b).im, (a).im + (b).re }
#define C_SIM(a, b) (Complex){ (a).re + (b).im, (a).im - (b).re }

#define FFT_CALC_FUNC(n, q)                                                     \
void CQT_fft_calc_ ## n(ShowCQT* cqt, Complex* v) {                             \
    const Complex* e2 = cqt->exp_tbl + 2*q;                                     \
    const Complex* e3 = cqt->exp_tbl + 3*q;                                     \
    const Complex* e1 = cqt->exp_tbl + 4*q;                                     \
    Complex v0, v1, v2, v3;                                                     \
    Complex a02, a13, s02, s13;                                                 \
                                                                                \
    CQT_fft_calc_ ## q(cqt, v);                                                 \
    CQT_fft_calc_ ## q(cqt, q+v);                                               \
    CQT_fft_calc_ ## q(cqt, 2*q+v);                                             \
    CQT_fft_calc_ ## q(cqt, 3*q+v);                                             \
                                                                                \
    v0 = v[0];                                                                  \
    v2 = v[q]; /* bit reversed */                                               \
    v1 = v[2*q];                                                                \
    v3 = v[3*q];                                                                \
    a02 = C_ADD(v0, v2);                                                        \
    s02 = C_SUB(v0, v2);                                                        \
    a13 = C_ADD(v1, v3);                                                        \
    s13 = C_SUB(v1, v3);                                                        \
    v[0] = C_ADD(a02, a13);                                                     \
    v[q] = C_SIM(s02, s13);                                                     \
    v[2*q] = C_SUB(a02, a13);                                                   \
    v[3*q] = C_AIM(s02, s13);                                                   \
                                                                                \
    for (int x = 1; x < q; x++) {                                               \
        v0 = v[x];                                                              \
        v2 = C_MUL(e2[x], v[q+x]); /* bit reversed */                           \
        v1 = C_MUL(e1[x], v[2*q+x]);                                            \
        v3 = C_MUL(e3[x], v[3*q+x]);                                            \
        a02 = C_ADD(v0, v2);                                                    \
        s02 = C_SUB(v0, v2);                                                    \
        a13 = C_ADD(v1, v3);                                                    \
        s13 = C_SUB(v1, v3);                                                    \
        v[x] = C_ADD(a02, a13);                                                 \
        v[q+x] = C_SIM(s02, s13);                                               \
        v[2*q+x] = C_SUB(a02, a13);                                             \
        v[3*q+x] = C_AIM(s02, s13);                                             \
    }                                                                           \
}

ALWAYS_INLINE void CQT_fft_calc_1(ShowCQT* cqt, Complex* v) { }
ALWAYS_INLINE void CQT_fft_calc_2(ShowCQT* cqt, Complex* v) {
    Complex v0 = v[0], v1 = v[1];
    v[0] = C_ADD(v0, v1);
    v[1] = C_SUB(v0, v1);
}

FFT_CALC_FUNC(4, 1)
FFT_CALC_FUNC(8, 2)
FFT_CALC_FUNC(16, 4)
FFT_CALC_FUNC(32, 8)
FFT_CALC_FUNC(64, 16)
FFT_CALC_FUNC(128, 32)
FFT_CALC_FUNC(256, 64)
FFT_CALC_FUNC(512, 128)
FFT_CALC_FUNC(1024, 256)
FFT_CALC_FUNC(2048, 512)
FFT_CALC_FUNC(4096, 1024)
FFT_CALC_FUNC(8192, 2048)
FFT_CALC_FUNC(16384, 4096)
FFT_CALC_FUNC(32768, 8192)

void CQT_fft_calc(ShowCQT* cqt, Complex* v, int n) {
    switch (n) {
        case  1024: CQT_fft_calc_1024(cqt,v); break;
        case  2048: CQT_fft_calc_2048(cqt,v); break;
        case  4096: CQT_fft_calc_4096(cqt,v); break;
        case  8192: CQT_fft_calc_8192(cqt,v); break;
        case 16384: CQT_fft_calc_16384(cqt,v); break;
        case 32768: CQT_fft_calc_32768(cqt,v); break;
    }
}

int CQT_init(ShowCQT* cqt, int rate, int width, int height, float bar_v, float sono_v, int super) {
#ifdef SHOWCQT_RENDERER
    if (height <= 0 || height > CQT_MAX_HEIGHT)
        return 0;
#endif
    if (width <= 0 || width > CQT_MAX_WIDTH)
        return 0;

    cqt->width = width;
    cqt->aligned_width = width;
#ifdef SHOWCQT_RENDERER
    cqt->height = height;

    cqt->bar_v = (bar_v > CQT_MAX_VOL) ? CQT_MAX_VOL : (bar_v > CQT_MIN_VOL) ? bar_v : CQT_MIN_VOL;
    cqt->sono_v = (sono_v > CQT_MAX_VOL) ? CQT_MAX_VOL : (sono_v > CQT_MIN_VOL) ? sono_v : CQT_MIN_VOL;
#endif
    if (rate < 8000 || rate > 100000)
        return 0;

    int bits = ceil(log(rate * 0.33)/ M_LN2);
    if (bits > 20 || bits < 10)
        return 0;
    cqt->fft_size = 1 << bits;
    if (cqt->fft_size > CQT_MAX_FFT_SIZE)
        return 0;

    CQT_gen_perm_tbl(cqt, bits - 2);
    CQT_gen_exp_tbl(cqt);

    cqt->attack_size = ceil(rate * 0.033);
    for (int x = 0; x < cqt->attack_size; x++) {
        double y = M_PI * x / (rate * 0.033);
        cqt->attack_tbl[x] = 0.355768 + 0.487396 * cos(y) + 0.144232 * cos(2*y) + 0.012604 * cos(3*y);
    }

    cqt->t_size = cqt->width * (1 + !!super);
    double log_base = log(20.01523126408007475);
    double log_end = log(20495.59681441799654);
    for (int f = 0, idx = 0; f < cqt->t_size; f++) {
        double freq = exp(log_base + (f + 0.5) * (log_end - log_base) * (1.0/cqt->t_size));

        if (freq >= 0.5 * rate) {
            cqt->kernel_index[f].len = 0;
            cqt->kernel_index[f].start = 0;
            continue;
        }

        double tlen = 384*0.33 / (384/0.17 + 0.33*freq/(1-0.17)) + 384*0.33 / (0.33*freq/0.17 + 384/(1-0.17));
        double flen = 8.0 * cqt->fft_size / (tlen * rate);
        double center = freq * cqt->fft_size / rate;
        int start = ceil(center - 0.5*flen);
        int end = floor(center + 0.5*flen);
        int len = end - start + 1;

        if (idx + len + 1000 > CQT_MAX_KERNEL_SIZE)
            return 0;
        cqt->kernel_index[f].len = len;
        cqt->kernel_index[f].start = start;

        for (int x = start; x < start + len; x++) {
            if (x > end) {
                cqt->kernel[idx+x-start] = 0;
                continue;
            }
            int sign = (x & 1) ? (-1) : 1;
            double y = 2.0 * M_PI * (x - center) * (1.0 / flen);
            double w = 0.355768 + 0.487396 * cos(y) + 0.144232 * cos(2*y) + 0.012604 * cos(3*y);
            w *= sign * (1.0/cqt->fft_size);
            cqt->kernel[idx+x-start] = w;
        }

        idx += len;
    }
    return cqt->fft_size;
}

Complex CQT_cqt_calc(ShowCQT* cqt, const float *kernel, int start, int len) {
    Complex a = { 0, 0 }, b = { 0, 0 };

    for (int m = 0, i = start, j = cqt->fft_size - start; m < len; m++, i++, j--) {
        float u = kernel[m];
        a.re += u * cqt->fft_buf[i].re;
        a.im += u * cqt->fft_buf[i].im;
        b.re += u * cqt->fft_buf[j].re;
        b.im += u * cqt->fft_buf[j].im;
    }

    Complex v0 = { a.re + b.re, a.im - b.im };
    Complex v1 = { b.im + a.im, b.re - a.re };
    float r0 = v0.re*v0.re + v0.im*v0.im;
    float r1 = v1.re*v1.re + v1.im*v1.im;
    return (Complex){ r0, r1 };
}

void CQT_calc(ShowCQT* cqt) {
    int fft_size_h = cqt->fft_size >> 1;
    int fft_size_q = cqt->fft_size >> 2;
    int shift = fft_size_h - cqt->attack_size;

    for (int x = 0; x < cqt->attack_size; x++) {
        int i = 4 * cqt->perm_tbl[x];
        cqt->fft_buf[i] = (Complex){ cqt->input[0][shift+x], cqt->input[1][shift+x] };
        cqt->fft_buf[i+1].re = cqt->attack_tbl[x] * cqt->input[0][fft_size_h+shift+x];
        cqt->fft_buf[i+1].im = cqt->attack_tbl[x] * cqt->input[1][fft_size_h+shift+x];
        cqt->fft_buf[i+2] = (Complex){ cqt->input[0][fft_size_q+shift+x], cqt->input[1][fft_size_q+shift+x] };
        cqt->fft_buf[i+3] = (Complex){0,0};
    }

    for (int x = cqt->attack_size; x < fft_size_q; x++) {
        int i = 4 * cqt->perm_tbl[x];
        cqt->fft_buf[i] = (Complex){ cqt->input[0][shift+x], cqt->input[1][shift+x] };
        cqt->fft_buf[i+1] = (Complex){0,0};
        cqt->fft_buf[i+2] = (Complex){ cqt->input[0][fft_size_q+shift+x], cqt->input[1][fft_size_q+shift+x] };
        cqt->fft_buf[i+3] = (Complex){0,0};
    }

    CQT_fft_calc(cqt, cqt->fft_buf, cqt->fft_size);

    for (int x = 0, m = 0; x < cqt->t_size; x++) {
        int len = cqt->kernel_index[x].len;
        int start = cqt->kernel_index[x].start;
        if (!len) {
#ifdef SHOWCQT_RENDERER
            cqt->color_buf[x] = (ColorF){0,0,0,0};
#endif
            continue;
        }

        Complex r = CQT_cqt_calc(cqt, cqt->kernel + m, start, len);

        cqt->magOutput[x] = sqrtf(r.re*r.re + r.im*r.im);

#ifdef SHOWCQT_RENDERER
        cqt->color_buf[x].r = sqrtf(cqt->sono_v * sqrtf(r.re));
        cqt->color_buf[x].g = sqrtf(cqt->sono_v * sqrtf(0.5f * (r.re + r.im)));
        cqt->color_buf[x].b = sqrtf(cqt->sono_v * sqrtf(r.im));
        cqt->color_buf[x].h = cqt->bar_v * sqrtf(0.5f * (r.re + r.im));
#endif

        m += len;
    }

#ifdef SHOWCQT_RENDERER
    if (cqt->t_size != cqt->width) {
        for (int x = 0; x < cqt->width; x++) {
            cqt->color_buf[x].r = 0.5f * (cqt->color_buf[2*x].r + cqt->color_buf[2*x+1].r);
            cqt->color_buf[x].g = 0.5f * (cqt->color_buf[2*x].g + cqt->color_buf[2*x+1].g);
            cqt->color_buf[x].b = 0.5f * (cqt->color_buf[2*x].b + cqt->color_buf[2*x+1].b);
            cqt->color_buf[x].h = 0.5f * (cqt->color_buf[2*x].h + cqt->color_buf[2*x+1].h);
        }
    }

    cqt->prerender = 1;
#endif
}

#ifdef SHOWCQT_RENDERER
void CQT_prerender(ShowCQT* cqt) {
    for (int x = 0; x < cqt->width; x++) {
        ColorF *c = cqt->color_buf;
        c[x].r = 255.5f * (c[x].r >= 0.0f ? (c[x].r <= 1.0f ? c[x].r : 1.0f) : 0.0f);
        c[x].g = 255.5f * (c[x].g >= 0.0f ? (c[x].g <= 1.0f ? c[x].g : 1.0f) : 0.0f);
        c[x].b = 255.5f * (c[x].b >= 0.0f ? (c[x].b <= 1.0f ? c[x].b : 1.0f) : 0.0f);
        c[x].h = c[x].h >= 0.0f ? c[x].h : 0.0f;
    }

    for (int x = 0; x < cqt->aligned_width; x++)
        cqt->rcp_h_buf[x] = 1.0f / (cqt->color_buf[x].h + 0.0001f);

    cqt->prerender = 0;
}

void CQT_render_line_alpha(ShowCQT* cqt, int y, uint8_t alpha) {
    if (cqt->prerender)
        CQT_prerender(cqt);

    int a = ((int) alpha) << 24;

    if (y >= 0 && y < cqt->height) {
        float ht = (cqt->height - y) / (float) cqt->height;
        for (int x = 0; x < cqt->width; x++) {
            if (cqt->color_buf[x].h <= ht) {
                cqt->output[x] = a;
            } else {
                float mul = (cqt->color_buf[x].h - ht) * cqt->rcp_h_buf[x];
                int r = mul * cqt->color_buf[x].r;
                int g = mul * cqt->color_buf[x].g;
                int b = mul * cqt->color_buf[x].b;
                g = g << 8;
                b = b << 16;
                cqt->output[x] = (r | g) | (b | a);
            }
        }
    } else {
        for (int x = 0; x < cqt->width; x++) {
            int r = cqt->color_buf[x].r;
            int g = cqt->color_buf[x].g;
            int b = cqt->color_buf[x].b;
            g = g << 8;
            b = b << 16;
            cqt->output[x] = (r | g) | (b | a);
        }
    }
}

void CQT_render_line_opaque(ShowCQT* cqt, int y) {
    CQT_render_line_alpha(cqt, y, 255);
}

void CQT_set_volume(ShowCQT* cqt, float bar_v, float sono_v) {
    cqt->bar_v = (bar_v > CQT_MAX_VOL) ? CQT_MAX_VOL : (bar_v > CQT_MIN_VOL) ? bar_v : CQT_MIN_VOL;
    cqt->sono_v = (sono_v > CQT_MAX_VOL) ? CQT_MAX_VOL : (sono_v > CQT_MIN_VOL) ? sono_v : CQT_MIN_VOL;
}

void CQT_set_height(ShowCQT* cqt, int height) {
    cqt->height = (height > CQT_MAX_HEIGHT) ? CQT_MAX_HEIGHT : (height > 1) ? height : 1;
}
#endif
