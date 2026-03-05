#ifndef __SALMON2_ARRAY_INDEX_H__
#define __SALMON2_ARRAY_INDEX_H__

#define ARRAY_INDEX_1D(C_i, Fort_i_start) ((C_i) - (Fort_i_start))
#define ARRAY_INDEX_2D(C_i, C_j, Fort_i_start, Fort_i_end, Fort_j_start) \
  (((C_i) - (Fort_i_start)) \
   + ((C_j) - (Fort_j_start)) * ((Fort_i_end) - (Fort_i_start) + 1))
#define ARRAY_INDEX_3D(C_i, C_j, C_k, Fort_i_start, Fort_i_end, Fort_j_start, Fort_j_end, Fort_k_start) \
  (((C_i) - (Fort_i_start)) \
   + ((C_j) - (Fort_j_start)) * ((Fort_i_end) - (Fort_i_start) + 1) \
   + ((C_k) - (Fort_k_start)) * ((Fort_i_end) - (Fort_i_start) + 1) * ((Fort_j_end) - (Fort_j_start) + 1))
#define ARRAY_INDEX_7D(C0, C1, C2, C3, C4, C5, C6, F0s, F0e, F1s, F1e, F2s, F2e, F3s, F3e, F4s, F4e, F5s, F5e, F6s) \
  (((C0) - (F0s)) \
  + ((C1) - (F1s)) * ((F0e) - (F0s) + 1) \
  + ((C2) - (F2s)) * ((F0e) - (F0s) + 1) * ((F1e) - (F1s) + 1) \
  + ((C3) - (F3s)) * ((F0e) - (F0s) + 1) * ((F1e) - (F1s) + 1)  * ((F2e) - (F2s) + 1) \
  + ((C4) - (F4s)) * ((F0e) - (F0s) + 1) * ((F1e) - (F1s) + 1)  * ((F2e) - (F2s) + 1) * ((F3e) - (F3s) + 1) \
  + ((C5) - (F5s)) * ((F0e) - (F0s) + 1) * ((F1e) - (F1s) + 1)  * ((F2e) - (F2s) + 1) * ((F3e) - (F3s) + 1) * ((F4e) - (F4s) + 1) \
  + ((C6) - (F6s)) * ((F0e) - (F0s) + 1) * ((F1e) - (F1s) + 1)  * ((F2e) - (F2s) + 1) * ((F3e) - (F3s) + 1) * ((F4e) - (F4s) + 1) * ((F5e) - (F5s) + 1) \
   )

#endif
