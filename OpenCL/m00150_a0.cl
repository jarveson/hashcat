/**
 * Author......: See docs/credits.txt
 * License.....: MIT
 */

#define NEW_SIMD_CODE

#include "inc_vendor.cl"
#include "inc_hash_constants.h"
#include "inc_hash_functions.cl"
#include "inc_types.cl"
#include "inc_common.cl"
#include "inc_rp.h"
#include "inc_rp.cl"
#include "inc_simd.cl"
#include "inc_hash_sha1.cl"

__kernel void m00150_m04 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  /**
   * modifier
   */

  const u32 lid = get_local_id (0);

  /**
   * base
   */

  const u32 gid = get_global_id (0);

  /**
   * salt
   */

  __local salt_t s_salt_buf[1];

  if (lid == 0)
  {
    s_salt_buf[0] = salt_bufs[salt_pos];
  }

  barrier (CLK_LOCAL_MEM_FENCE);

  if (gid >= gid_max) return;

  u32 pw_buf0[4];
  u32 pw_buf1[4];

  pw_buf0[0] = pws[gid].i[0];
  pw_buf0[1] = pws[gid].i[1];
  pw_buf0[2] = pws[gid].i[2];
  pw_buf0[3] = pws[gid].i[3];
  pw_buf1[0] = pws[gid].i[4];
  pw_buf1[1] = pws[gid].i[5];
  pw_buf1[2] = pws[gid].i[6];
  pw_buf1[3] = pws[gid].i[7];

  const u32 pw_len = pws[gid].pw_len;

  const u32 salt_len = salt_bufs[salt_pos].salt_len;

  /**
   * loop
   */

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    u32x w0[4] = { 0 };
    u32x w1[4] = { 0 };

    apply_rules_vect (pw_buf0, pw_buf1, pw_len, rules_buf, il_pos, w0, w1);

    sha1_ctx_t sha1_ctx;
    sha1_init(&sha1_ctx);

    //printf("key?: 0x%x, 0x%x, 0x%x, 0x%x\n", swap32(w0[0]), swap32(w0[1]), w0[2], w0[3]);
    //printf("saltLen: 0x%x\n", salt_len);
    //printf("salt?: 0x%x, 0x%x, 0x%x, 0x%x\n", s_salt_buf[0].salt_buf[0], s_salt_buf[0].salt_buf[1], s_salt_buf[0].salt_buf[2], s_salt_buf[0].salt_buf[3]);
    //printf("salt?: 0x%x, 0x%x, 0x%x, 0x%x\n", salt_bufs[salt_pos].salt_buf[0], salt_bufs[salt_pos].salt_buf[1], salt_bufs[salt_pos].salt_buf[2], salt_bufs[salt_pos].salt_buf[3]);

    u32 ipadData[16];
    ipadData[0] = swap32 (w0[0]) ^ 0x36363636;
    ipadData[1] = swap32 (w0[1]) ^ 0x36363636;
    ipadData[2] = swap32 (w0[2]) ^ 0x36363636;
    ipadData[3] = swap32 (w0[3]) ^ 0x36363636;
    ipadData[4] = swap32 (w1[0]) ^ 0x36363636;
    ipadData[5] = swap32 (w1[1]) ^ 0x36363636;
    ipadData[6] = swap32 (w1[2]) ^ 0x36363636;
    ipadData[7] = swap32 (w1[3]) ^ 0x36363636;
    ipadData[8] = 0x36363636;
    ipadData[9] = 0x36363636;
    ipadData[10] = 0x36363636;
    ipadData[11] = 0x36363636;
    ipadData[12] = 0x36363636;
    ipadData[13] = 0x36363636;
    ipadData[14] = 0x36363636;
    ipadData[15] = 0x36363636;

    u32 opadData[16];
    opadData[0] = swap32 (w0[0]) ^ 0x5c5c5c5c;
    opadData[1] = swap32 (w0[1]) ^ 0x5c5c5c5c;
    opadData[2] = swap32 (w0[2]) ^ 0x5c5c5c5c;
    opadData[3] = swap32 (w0[3]) ^ 0x5c5c5c5c;
    opadData[4] = swap32 (w1[0]) ^ 0x5c5c5c5c;
    opadData[5] = swap32 (w1[1]) ^ 0x5c5c5c5c;
    opadData[6] = swap32 (w1[2]) ^ 0x5c5c5c5c;
    opadData[7] = swap32 (w1[3]) ^ 0x5c5c5c5c;
    opadData[8] = 0x5c5c5c5c;
    opadData[9] = 0x5c5c5c5c;
    opadData[10] = 0x5c5c5c5c;
    opadData[11] = 0x5c5c5c5c;
    opadData[12] = 0x5c5c5c5c;
    opadData[13] = 0x5c5c5c5c;
    opadData[14] = 0x5c5c5c5c;
    opadData[15] = 0x5c5c5c5c;

    sha1_update_64(&sha1_ctx, &ipadData[0], &ipadData[4], &ipadData[8], &ipadData[12], 64);
    //printf("statesha1: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.h[0], sha1_ctx.h[1], sha1_ctx.h[2], sha1_ctx.h[3], sha1_ctx.h[4]);

    //sha1_update(&sha1_ctx, salt_buf0, salt_len);
    sha1_update_local(&sha1_ctx, s_salt_buf[0].salt_buf, salt_len);

    //printf("statesha1: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.h[0], sha1_ctx.h[1], sha1_ctx.h[2], sha1_ctx.h[3], sha1_ctx.h[4]);

    sha1_final(&sha1_ctx);
    //printf("statesha1: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.h[0], sha1_ctx.h[1], sha1_ctx.h[2], sha1_ctx.h[3], sha1_ctx.h[4]);

    /*printf("wstatesha1: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.w0[0], sha1_ctx.w0[1], sha1_ctx.w0[2], sha1_ctx.w0[3], sha1_ctx.w0[4]);
    printf("wstatesha2: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.w1[0], sha1_ctx.w1[1], sha1_ctx.w1[2], sha1_ctx.w1[3], sha1_ctx.w1[4]);
    printf("wstatesha3: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.w2[0], sha1_ctx.w2[1], sha1_ctx.w2[2], sha1_ctx.w2[3], sha1_ctx.w2[4]);
    printf("wstatesha4: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.w3[0], sha1_ctx.w3[1], sha1_ctx.w3[2], sha1_ctx.w3[3], sha1_ctx.w3[4]);*/


    sha1_ctx_t opad_ctx;
    sha1_init(&opad_ctx);

    sha1_update_64(&opad_ctx, &opadData[0], &opadData[4], &opadData[8], &opadData[12], 64);
    //printf("stateopad2: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", opad_ctx.h[0], opad_ctx.h[1], opad_ctx.h[2], opad_ctx.h[3], opad_ctx.h[4]);

    u32 tmpData[16];
    tmpData[0] = sha1_ctx.h[0];
    tmpData[1] = sha1_ctx.h[1];
    tmpData[2] = sha1_ctx.h[2];
    tmpData[3] = sha1_ctx.h[3];
    tmpData[4] = sha1_ctx.h[4];
    tmpData[5] = 0;
    tmpData[6] = 0;
    tmpData[7] = 0;
    tmpData[8] = 0;
    tmpData[9] = 0;
    tmpData[10] = 0;
    tmpData[11] = 0;
    tmpData[12] = 0;
    tmpData[13] = 0;
    tmpData[14] = 0;
    tmpData[15] = 0;

    sha1_update(&opad_ctx, tmpData, 20);
    //printf("stateopadupdate: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", opad_ctx.h[0], opad_ctx.h[1], opad_ctx.h[2], opad_ctx.h[3], opad_ctx.h[4]);

    sha1_final(&opad_ctx);

    COMPARE_M_SIMD (opad_ctx.h[3], opad_ctx.h[4], opad_ctx.h[2], opad_ctx.h[1]);
  }
}

__kernel void m00150_m08 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m00150_m16 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m00150_s04 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  /**
   * modifier
   */

  const u32 lid = get_local_id (0);

  /**
   * base
   */

  const u32 gid = get_global_id (0);

  /**
   * salt
   */

  __local salt_t s_salt_buf[1];

  if (lid == 0)
  {
    s_salt_buf[0] = salt_bufs[salt_pos];
  }

  barrier (CLK_LOCAL_MEM_FENCE);

  if (gid >= gid_max) return;

  u32 pw_buf0[4];
  u32 pw_buf1[4];

  pw_buf0[0] = pws[gid].i[0];
  pw_buf0[1] = pws[gid].i[1];
  pw_buf0[2] = pws[gid].i[2];
  pw_buf0[3] = pws[gid].i[3];
  pw_buf1[0] = pws[gid].i[4];
  pw_buf1[1] = pws[gid].i[5];
  pw_buf1[2] = pws[gid].i[6];
  pw_buf1[3] = pws[gid].i[7];

  const u32 pw_len = pws[gid].pw_len;

  /**
   * salt
   */

  /*u32 salt_buf0[0xAB0];

  for (int i=0; i<0xAB0; i+=16) {
    salt_buf0[i+0] = salt_bufs[salt_pos].salt_buf[i+ 0];
    salt_buf0[i+1] = salt_bufs[salt_pos].salt_buf[ i+1];
    salt_buf0[i+2] = salt_bufs[salt_pos].salt_buf[ i+2];
    salt_buf0[i+3] = salt_bufs[salt_pos].salt_buf[ i+3];
    salt_buf0[i+4] = salt_bufs[salt_pos].salt_buf[ i+4];
    salt_buf0[i+5] = salt_bufs[salt_pos].salt_buf[ i+5];
    salt_buf0[i+6] = salt_bufs[salt_pos].salt_buf[ i+6];
    salt_buf0[i+7] = salt_bufs[salt_pos].salt_buf[ i+7];
    salt_buf0[i+8] = salt_bufs[salt_pos].salt_buf[ i+8];
    salt_buf0[i+9] = salt_bufs[salt_pos].salt_buf[ i+9];
    salt_buf0[i+10] = salt_bufs[salt_pos].salt_buf[i+10];
    salt_buf0[i+11] = salt_bufs[salt_pos].salt_buf[i+11];
    salt_buf0[i+12] = salt_bufs[salt_pos].salt_buf[i+12];
    salt_buf0[i+13] = salt_bufs[salt_pos].salt_buf[i+13];
    salt_buf0[i+14] = salt_bufs[salt_pos].salt_buf[i+14];
    salt_buf0[i+15] = salt_bufs[salt_pos].salt_buf[i+15];
  }*/


  const u32 salt_len = salt_bufs[salt_pos].salt_len;

  /**
   * digest
   */

  const u32 search[4] =
  {
    digests_buf[digests_offset].digest_buf[DGST_R0],
    digests_buf[digests_offset].digest_buf[DGST_R1],
    digests_buf[digests_offset].digest_buf[DGST_R2],
    digests_buf[digests_offset].digest_buf[DGST_R3]
  };

  /**
   * loop
   */

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    u32x w0[4] = { 0 };
    u32x w1[4] = { 0 };

    apply_rules_vect (pw_buf0, pw_buf1, pw_len, rules_buf, il_pos, w0, w1);

    sha1_ctx_t sha1_ctx;
    sha1_init(&sha1_ctx);

    //printf("key?: 0x%x, 0x%x, 0x%x, 0x%x\n", swap32(w0[0]), swap32(w0[1]), w0[2], w0[3]);
    //printf("saltLen: 0x%x\n", salt_len);
    //printf("salt?: 0x%x, 0x%x, 0x%x, 0x%x\n", s_salt_buf[0].salt_buf[0], s_salt_buf[0].salt_buf[1], s_salt_buf[0].salt_buf[2], s_salt_buf[0].salt_buf[3]);
    //printf("salt?: 0x%x, 0x%x, 0x%x, 0x%x\n", salt_bufs[salt_pos].salt_buf[0], salt_bufs[salt_pos].salt_buf[1], salt_bufs[salt_pos].salt_buf[2], salt_bufs[salt_pos].salt_buf[3]);

    u32 ipadData[16];
    ipadData[0] = swap32 (w0[0]) ^ 0x36363636;
    ipadData[1] = swap32 (w0[1]) ^ 0x36363636;
    ipadData[2] = swap32 (w0[2]) ^ 0x36363636;
    ipadData[3] = swap32 (w0[3]) ^ 0x36363636;
    ipadData[4] = swap32 (w1[0]) ^ 0x36363636;
    ipadData[5] = swap32 (w1[1]) ^ 0x36363636;
    ipadData[6] = swap32 (w1[2]) ^ 0x36363636;
    ipadData[7] = swap32 (w1[3]) ^ 0x36363636;
    ipadData[8] = 0x36363636;
    ipadData[9] = 0x36363636;
    ipadData[10] = 0x36363636;
    ipadData[11] = 0x36363636;
    ipadData[12] = 0x36363636;
    ipadData[13] = 0x36363636;
    ipadData[14] = 0x36363636;
    ipadData[15] = 0x36363636;

    u32 opadData[16];
    opadData[0] = swap32 (w0[0]) ^ 0x5c5c5c5c;
    opadData[1] = swap32 (w0[1]) ^ 0x5c5c5c5c;
    opadData[2] = swap32 (w0[2]) ^ 0x5c5c5c5c;
    opadData[3] = swap32 (w0[3]) ^ 0x5c5c5c5c;
    opadData[4] = swap32 (w1[0]) ^ 0x5c5c5c5c;
    opadData[5] = swap32 (w1[1]) ^ 0x5c5c5c5c;
    opadData[6] = swap32 (w1[2]) ^ 0x5c5c5c5c;
    opadData[7] = swap32 (w1[3]) ^ 0x5c5c5c5c;
    opadData[8] = 0x5c5c5c5c;
    opadData[9] = 0x5c5c5c5c;
    opadData[10] = 0x5c5c5c5c;
    opadData[11] = 0x5c5c5c5c;
    opadData[12] = 0x5c5c5c5c;
    opadData[13] = 0x5c5c5c5c;
    opadData[14] = 0x5c5c5c5c;
    opadData[15] = 0x5c5c5c5c;

    sha1_update_64(&sha1_ctx, &ipadData[0], &ipadData[4], &ipadData[8], &ipadData[12], 64);
    //printf("statesha1: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.h[0], sha1_ctx.h[1], sha1_ctx.h[2], sha1_ctx.h[3], sha1_ctx.h[4]);

    //sha1_update(&sha1_ctx, salt_buf0, salt_len);
    sha1_update_local(&sha1_ctx, s_salt_buf[0].salt_buf, salt_len);

    //printf("statesha1: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.h[0], sha1_ctx.h[1], sha1_ctx.h[2], sha1_ctx.h[3], sha1_ctx.h[4]);

    sha1_final(&sha1_ctx);
    //printf("statesha1: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.h[0], sha1_ctx.h[1], sha1_ctx.h[2], sha1_ctx.h[3], sha1_ctx.h[4]);

    /*printf("wstatesha1: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.w0[0], sha1_ctx.w0[1], sha1_ctx.w0[2], sha1_ctx.w0[3], sha1_ctx.w0[4]);
    printf("wstatesha2: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.w1[0], sha1_ctx.w1[1], sha1_ctx.w1[2], sha1_ctx.w1[3], sha1_ctx.w1[4]);
    printf("wstatesha3: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.w2[0], sha1_ctx.w2[1], sha1_ctx.w2[2], sha1_ctx.w2[3], sha1_ctx.w2[4]);
    printf("wstatesha4: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", sha1_ctx.w3[0], sha1_ctx.w3[1], sha1_ctx.w3[2], sha1_ctx.w3[3], sha1_ctx.w3[4]);*/

    sha1_ctx_t opad_ctx;
    sha1_init(&opad_ctx);

    sha1_update_64(&opad_ctx, &opadData[0], &opadData[4], &opadData[8], &opadData[12], 64);
    //printf("stateopad2: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", opad_ctx.h[0], opad_ctx.h[1], opad_ctx.h[2], opad_ctx.h[3], opad_ctx.h[4]);

    u32 tmpData[16];
    tmpData[0] = sha1_ctx.h[0];
    tmpData[1] = sha1_ctx.h[1];
    tmpData[2] = sha1_ctx.h[2];
    tmpData[3] = sha1_ctx.h[3];
    tmpData[4] = sha1_ctx.h[4];
    tmpData[5] = 0;
    tmpData[6] = 0;
    tmpData[7] = 0;
    tmpData[8] = 0;
    tmpData[9] = 0;
    tmpData[10] = 0;
    tmpData[11] = 0;
    tmpData[12] = 0;
    tmpData[13] = 0;
    tmpData[14] = 0;
    tmpData[15] = 0;

    sha1_update(&opad_ctx, tmpData, 20);
    //printf("stateopadupdate: 0x%x, 0x%x, 0x%x, 0x%x, 0x%x\n", opad_ctx.h[0], opad_ctx.h[1], opad_ctx.h[2], opad_ctx.h[3], opad_ctx.h[4]);

    sha1_final(&opad_ctx);

    //printf("outputhash: 0x%x, 0x%x, 0x%x, 0x%x\n", opad_ctx.h[3], opad_ctx.h[4], opad_ctx.h[2], opad_ctx.h[1]);
    COMPARE_S_SIMD (opad_ctx.h[3], opad_ctx.h[4], opad_ctx.h[2], opad_ctx.h[1]);
  }
}

__kernel void m00150_s08 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m00150_s16 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}
