/*
 *  Copyright 2017 SALMON developers
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/* Hand-Code Vector processing for Knights Corner */

#include <complex.h>
#include <assert.h>
#include "./glue.h"

void stencil_c_tuned_seq_imp( const int PNLx
                            , const int PNLy
                            , const int PNLz
                            , const int NLx
                            , const int NLy
                            , const int NLz
                            , const int NHx
                            , const int NHy
                            , const int NHz
                            , int    const * restrict modx
                            , int    const * restrict mody
                            , int    const * restrict modz
                            , int    const            igs[restrict 3]
                            , int    const            ige[restrict 3]
                            , double const * restrict A
                            , double const            B[restrict NLz][NLy][NLx]
                            , double const            C[restrict 12]
                            , double const            D[restrict 12]
                            , double complex const    E[restrict PNLz][PNLy][PNLx]
                            , double complex          F[restrict PNLz][PNLy][PNLx]
)
{
  const int isx = igs[0], isy = igs[1], isz = igs[2];
  const int iex = ige[0], iey = ige[1], iez = ige[2];
  const int NSx = NHx * 2, NSy = NHy * 2, NSz = NHz * 2;

  __m512d at   = _mm512_set1_pd(*A);
  __m512d HALF = _mm512_set1_pd(-0.5);
  __m512i INV  = _mm512_set4_epi64(1LL << 63, 0, 1LL << 63, 0);

  __declspec(align(64)) double G[12];
  for(int n = 0 ; n < 12 ; ++n)
    G[n] = C[n] * -0.5;

  __m512i nly = _mm512_set1_epi32(PNLy);
  __m512i nlx = _mm512_set1_epi32(PNLx);
  __m512i nyz = _mm512_mask_blend_epi32(0xFF00, _mm512_set1_epi32(PNLy), _mm512_set1_epi32(PNLz));
  __m512i nlyx = _mm512_mask_mullo_epi32(nlx, 0xFF00, nlx, nly);

  __declspec(align(64)) int yz_table[16];
  __m512i *yz = (__m512i*) yz_table;


#pragma noprefetch
#pragma novector
  for(int iz = isz ; iz < iez ; ++iz)
  {
    __m512i tiz = _mm512_set1_epi32(iz);
    __m512i mzm = _mm512_loadu_prefetch_epi32(modz + (iz - 4 + NLz + NSz));
    __m512i mzp = _mm512_alignr_epi32(mzm, mzm, 1);
    __m512i zmp = _mm512_mask_blend_epi32(0xF0F0, mzm, mzp);
            zmp = _mm512_permute4f128_epi32(zmp, _MM_PERM_BADC);
  for(int iy = isy ; iy < iey ; ++iy)
  {
    __m512i tiy = _mm512_set1_epi32(iy);
    __m512i mym = _mm512_loadu_prefetch_epi32(mody + (iy - 4 + NLy + NSy));
    __m512i myp = _mm512_alignr_epi32(mym, mym, 1);
    __m512i ymp = _mm512_mask_blend_epi32(0xF0F0, mym, myp);
    __m512i uyz = _mm512_mask_blend_epi32(0xFF00, ymp, zmp);
    __m512i tyz = _mm512_mask_blend_epi32(0xFF00, tiy, tiz);

    double         const* b = &B[iz-NHz][iy-NHy][0];
    double complex const* e = &E[iz    ][iy    ][0];
    double complex      * f = &F[iz    ][iy    ][0];

    for(int ix = isx ; ix < iex ; ix += 4)
    {
      __m512i tix = _mm512_set1_epi32(ix);
      __m512i mm  = _mm512_sub_epi32(tyz, uyz);
      *yz = _mm512_sub_epi32(tix, _mm512_mullo_epi32(mm, nlyx));

      __m512d ex = _mm512_load_prefetch_pd(e + ix);
      __m512d tt = _mm512_setzero_pd();
      __m512d ut = _mm512_setzero_pd();

      __m512d m, p, bt, v0, v1, v2, v3, v4;

      /* x-dimension (unit stride) */
      {
        __m512i x0, x2;
        x0 = _mm512_load_prefetch_epi64(e + modx[ix - 4 + NLx + NSx]);
        x2 = _mm512_load_prefetch_epi64(e + modx[ix + 4 + NLx + NSx]);
        {
          m  = (__m512d) _mm512_alignr_epi32((__m512i) ex, x0, 12);
          p  = (__m512d) _mm512_alignr_epi32(x2, (__m512i) ex,  4);
          v4 = _mm512_sub_pd(p, m);
          v3 = _mm512_add_pd(p, m);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[0]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[0]), v3, tt);
        }
        {
          m  = (__m512d) _mm512_alignr_epi32((__m512i) ex, x0,  8);
          p  = (__m512d) _mm512_alignr_epi32(x2, (__m512i) ex,  8);
          v4 = _mm512_sub_pd(p, m);
          v3 = _mm512_add_pd(p, m);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[1]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[1]), v3, tt);
        }
        {
          m  = (__m512d) _mm512_alignr_epi32((__m512i) ex, x0,  4);
          p  = (__m512d) _mm512_alignr_epi32(x2, (__m512i) ex, 12);
          v4 = _mm512_sub_pd(p, m);
          v3 = _mm512_add_pd(p, m);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[2]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[2]), v3, tt);
        }
        {
          m  = (__m512d) x0;
          p  = (__m512d) x2;
          v4 = _mm512_sub_pd(p, m);
          v3 = _mm512_add_pd(p, m);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[3]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[3]), v3, tt);
        }
      }

      /* y-dimension (NLx stride) */
      {
#pragma unroll(4)
        for(int n = 0 ; n < 4 ; ++n) {
          m  = _mm512_load_prefetch_pd(e + yz_table[3-n]);
          p  = _mm512_load_prefetch_pd(e + yz_table[n+4]);
          v4 = _mm512_sub_pd(p, m);
          v3 = _mm512_add_pd(p, m);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[n+4]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[n+4]), v3, tt);
        }
      }

      /* z-dimension (NLy*NLx stride)  */
      {
#pragma unroll(4)
        for(int n = 0 ; n < 4 ; ++n) {
          m  = _mm512_load_prefetch_pd(e + yz_table[11-n]);
          p  = _mm512_load_prefetch_pd(e + yz_table[n+12]);
          v4 = _mm512_sub_pd(p, m);
          v3 = _mm512_add_pd(p, m);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[n+8]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[n+8]), v3, tt);
        }
      }

      bt = dcast_to_dcmplx(b + ix - NHx);
      v2 = _mm512_fmadd_pd(at, ex, tt);
      v4 = (__m512d) _mm512_shuffle_epi32((__m512i) ut, _MM_PERM_BADC);
      v3 = (__m512d) _mm512_xor_si512((__m512i) v4, INV);
      v1 = _mm512_add_pd(v2, v3);
      v0 = _mm512_fmadd_pd(bt, ex, v1);

      _mm512_storenrngo_pd(&f[ix], v0);
    } /* NLx */
  } /* NLy */
  } /* NLz */
}

/*
 * is_array: one origin
 * ie_array: one origin
 * is      : one origin
 * ie      : one origin
 */
void stencil_c_tuned_seq_( int            const            is_array[restrict 3]
                         , int            const            ie_array[restrict 3]
                         , int            const            is[restrict 3]
                         , int            const            ie[restrict 3]
                         , int            const * restrict modx
                         , int            const * restrict mody
                         , int            const * restrict modz
                         , int            const            igs_[restrict 3]
                         , int            const            ige_[restrict 3]
                         , double complex const * restrict E
                         , double complex       * restrict F
                         , double         const * restrict B
                         , double         const * restrict A_
                         , double         const            C[restrict 12]
                         , double         const            D[restrict 12]
)
{
#define INT_ABS(X) (X) < 0 ? -(X) : (X)
  const int PNLx = ie_array[0] - is_array[0] + 1;
  const int PNLy = ie_array[1] - is_array[1] + 1;
  const int PNLz = ie_array[2] - is_array[2] + 1;

  const int NLx  = ie[0] - is[0] + 1; // lattice
  const int NLy  = ie[1] - is[1] + 1;
  const int NLz  = ie[2] - is[2] + 1;

  const int NHx  = INT_ABS(is_array[0] - is[0]); // shadow
  const int NHy  = INT_ABS(is_array[1] - is[1]);
  const int NHz  = INT_ABS(is_array[2] - is[2]);

  int igs[3] = { NHx,          NHy,          NHz };
  int ige[3] = { igs[0] + NLx, igs[1] + NLy, igs[2] + NLz };

  for (int i = 0 ; i < 3 ; ++i)
  {
    if (is[i] < igs_[i])
      igs[i] = igs[i] + igs_[i] - is[i];
    if (ige_[i] < ie[i])
      ige[i] = igs[i] + ige_[i] - igs_[i] + 1;
  }

  if (ige[0] > NHx + NLx) ige[0] = NHx + NLx;
  if (ige[1] > NHy + NLy) ige[1] = NHy + NLy;
  if (ige[2] > NHz + NLz) ige[2] = NHz + NLz;

  assert(NLx % 4 == 0);
  assert(NHx == 4 || NHx == 0);
  assert(NHy == 4 || NHy == 0);
  assert(NHz == 4 || NHz == 0);
  assert(igs[0] % 4 == 0);
  assert(ige[0] % 4 == 0);
#undef INT_ABS

  stencil_c_tuned_seq_imp(PNLx, PNLy, PNLz, NLx, NLy, NLz, NHx, NHy, NHz, modx, mody, modz, igs, ige
                         , A_
                         , (double         const (* restrict)[NLy][NLx])(B)
                         , C
                         , D
                         , (double complex const (* restrict)[PNLy][PNLx])(E)
                         , (double complex       (* restrict)[PNLy][PNLx])(F)
                         );
}
