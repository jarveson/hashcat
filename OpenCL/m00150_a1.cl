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

  const u32 pw_l_len = pws[gid].pw_len;

  const u32 salt_len = salt_bufs[salt_pos].salt_len;

  /**
   * loop
   */

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    const u32x pw_r_len = pwlenx_create_combt (combs_buf, il_pos);

    const u32x pw_len = pw_l_len + pw_r_len;

    /**
     * concat password candidate
     */

    u32x wordl0[4] = { 0 };
    u32x wordl1[4] = { 0 };
    u32x wordl2[4] = { 0 };
    u32x wordl3[4] = { 0 };

    wordl0[0] = pw_buf0[0];
    wordl0[1] = pw_buf0[1];
    wordl0[2] = pw_buf0[2];
    wordl0[3] = pw_buf0[3];
    wordl1[0] = pw_buf1[0];
    wordl1[1] = pw_buf1[1];
    wordl1[2] = pw_buf1[2];
    wordl1[3] = pw_buf1[3];

    u32x wordr0[4] = { 0 };
    u32x wordr1[4] = { 0 };
    u32x wordr2[4] = { 0 };
    u32x wordr3[4] = { 0 };

    wordr0[0] = ix_create_combt (combs_buf, il_pos, 0);
    wordr0[1] = ix_create_combt (combs_buf, il_pos, 1);
    wordr0[2] = ix_create_combt (combs_buf, il_pos, 2);
    wordr0[3] = ix_create_combt (combs_buf, il_pos, 3);
    wordr1[0] = ix_create_combt (combs_buf, il_pos, 4);
    wordr1[1] = ix_create_combt (combs_buf, il_pos, 5);
    wordr1[2] = ix_create_combt (combs_buf, il_pos, 6);
    wordr1[3] = ix_create_combt (combs_buf, il_pos, 7);

    if (combs_mode == COMBINATOR_MODE_BASE_LEFT)
    {
      switch_buffer_by_offset_le_VV (wordr0, wordr1, wordr2, wordr3, pw_l_len);
    }
    else
    {
      switch_buffer_by_offset_le_VV (wordl0, wordl1, wordl2, wordl3, pw_r_len);
    }

    u32x w0[4];
    u32x w1[4];
    u32x w2[4];
    u32x w3[4];

    w0[0] = wordl0[0] | wordr0[0];
    w0[1] = wordl0[1] | wordr0[1];
    w0[2] = wordl0[2] | wordr0[2];
    w0[3] = wordl0[3] | wordr0[3];
    w1[0] = wordl1[0] | wordr1[0];
    w1[1] = wordl1[1] | wordr1[1];
    w1[2] = wordl1[2] | wordr1[2];
    w1[3] = wordl1[3] | wordr1[3];
    w2[0] = wordl2[0] | wordr2[0];
    w2[1] = wordl2[1] | wordr2[1];
    w2[2] = wordl2[2] | wordr2[2];
    w2[3] = wordl2[3] | wordr2[3];
    w3[0] = wordl3[0] | wordr3[0];
    w3[1] = wordl3[1] | wordr3[1];
    w3[2] = wordl3[2] | wordr3[2];
    w3[3] = wordl3[3] | wordr3[3];

    w0[0] = swap32 (w0[0]);
    w0[1] = swap32 (w0[1]);
    w0[2] = swap32 (w0[2]);
    w0[3] = swap32 (w0[3]);
    w1[0] = swap32 (w1[0]);
    w1[1] = swap32 (w1[1]);
    w1[2] = swap32 (w1[2]);
    w1[3] = swap32 (w1[3]);
    w2[0] = swap32 (w2[0]);
    w2[1] = swap32 (w2[1]);
    w2[2] = swap32 (w2[2]);
    w2[3] = swap32 (w2[3]);
    w3[0] = swap32 (w3[0]);
    w3[1] = swap32 (w3[1]);
    w3[2] = swap32 (w3[2]);
    w3[3] = swap32 (w3[3]);

    
    sha1_ctx_t sha1_ctx;
    sha1_init(&sha1_ctx);

    //printf("key?: 0x%x, 0x%x, 0x%x, 0x%x\n", swap32(w0[0]), swap32(w0[1]), w0[2], w0[3]);
    //printf("saltLen: 0x%x\n", salt_len);
    //printf("salt?: 0x%x, 0x%x, 0x%x, 0x%x\n", s_salt_buf[0].salt_buf[0], s_salt_buf[0].salt_buf[1], s_salt_buf[0].salt_buf[2], s_salt_buf[0].salt_buf[3]);
    //printf("salt?: 0x%x, 0x%x, 0x%x, 0x%x\n", salt_bufs[salt_pos].salt_buf[0], salt_bufs[salt_pos].salt_buf[1], salt_bufs[salt_pos].salt_buf[2], salt_bufs[salt_pos].salt_buf[3]);

    u32x ipadData[16];
    ipadData[0] = w0[0] ^ 0x36363636;
    ipadData[1] = w0[1] ^ 0x36363636;
    ipadData[2] = w0[2] ^ 0x36363636;
    ipadData[3] = w0[3] ^ 0x36363636;
    ipadData[4] = w1[0] ^ 0x36363636;
    ipadData[5] = w1[1] ^ 0x36363636;
    ipadData[6] = w1[2] ^ 0x36363636;
    ipadData[7] = w1[3] ^ 0x36363636;
    ipadData[8] = w2[0] ^ 0x36363636;
    ipadData[9] = w2[1] ^ 0x36363636;
    ipadData[10] = w2[2] ^ 0x36363636;
    ipadData[11] = w2[3] ^ 0x36363636;
    ipadData[12] = w3[0] ^ 0x36363636;
    ipadData[13] = w3[1] ^ 0x36363636;
    ipadData[14] = w3[2] ^ 0x36363636;
    ipadData[15] = w3[3] ^ 0x36363636;

    u32x opadData[16];
    opadData[0] = w0[0] ^ 0x5c5c5c5c;
    opadData[1] = w0[1] ^ 0x5c5c5c5c;
    opadData[2] = w0[2] ^ 0x5c5c5c5c;
    opadData[3] = w0[3] ^ 0x5c5c5c5c;
    opadData[4] = w1[0] ^ 0x5c5c5c5c;
    opadData[5] = w1[1] ^ 0x5c5c5c5c;
    opadData[6] = w1[2] ^ 0x5c5c5c5c;
    opadData[7] = w1[3] ^ 0x5c5c5c5c;
    opadData[8] = w2[0] ^ 0x5c5c5c5c;
    opadData[9] = w2[1] ^ 0x5c5c5c5c;
    opadData[10] = w2[2] ^ 0x5c5c5c5c;
    opadData[11] = w2[3] ^ 0x5c5c5c5c;
    opadData[12] = w3[0] ^ 0x5c5c5c5c;
    opadData[13] = w3[1] ^ 0x5c5c5c5c;
    opadData[14] = w3[2] ^ 0x5c5c5c5c;
    opadData[15] = w2[3] ^ 0x5c5c5c5c;

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

    u32x tmpData[16];
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

  const u32 pw_l_len = pws[gid].pw_len;

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
    const u32x pw_r_len = pwlenx_create_combt (combs_buf, il_pos);

    const u32x pw_len = pw_l_len + pw_r_len;

    /**
     * concat password candidate
     */

    u32x wordl0[4] = { 0 };
    u32x wordl1[4] = { 0 };
    u32x wordl2[4] = { 0 };
    u32x wordl3[4] = { 0 };

    wordl0[0] = pw_buf0[0];
    wordl0[1] = pw_buf0[1];
    wordl0[2] = pw_buf0[2];
    wordl0[3] = pw_buf0[3];
    wordl1[0] = pw_buf1[0];
    wordl1[1] = pw_buf1[1];
    wordl1[2] = pw_buf1[2];
    wordl1[3] = pw_buf1[3];

    u32x wordr0[4] = { 0 };
    u32x wordr1[4] = { 0 };
    u32x wordr2[4] = { 0 };
    u32x wordr3[4] = { 0 };

    wordr0[0] = ix_create_combt (combs_buf, il_pos, 0);
    wordr0[1] = ix_create_combt (combs_buf, il_pos, 1);
    wordr0[2] = ix_create_combt (combs_buf, il_pos, 2);
    wordr0[3] = ix_create_combt (combs_buf, il_pos, 3);
    wordr1[0] = ix_create_combt (combs_buf, il_pos, 4);
    wordr1[1] = ix_create_combt (combs_buf, il_pos, 5);
    wordr1[2] = ix_create_combt (combs_buf, il_pos, 6);
    wordr1[3] = ix_create_combt (combs_buf, il_pos, 7);

    if (combs_mode == COMBINATOR_MODE_BASE_LEFT)
    {
      switch_buffer_by_offset_le_VV (wordr0, wordr1, wordr2, wordr3, pw_l_len);
    }
    else
    {
      switch_buffer_by_offset_le_VV (wordl0, wordl1, wordl2, wordl3, pw_r_len);
    }

    u32x w0[4];
    u32x w1[4];
    u32x w2[4];
    u32x w3[4];

    w0[0] = wordl0[0] | wordr0[0];
    w0[1] = wordl0[1] | wordr0[1];
    w0[2] = wordl0[2] | wordr0[2];
    w0[3] = wordl0[3] | wordr0[3];
    w1[0] = wordl1[0] | wordr1[0];
    w1[1] = wordl1[1] | wordr1[1];
    w1[2] = wordl1[2] | wordr1[2];
    w1[3] = wordl1[3] | wordr1[3];
    w2[0] = wordl2[0] | wordr2[0];
    w2[1] = wordl2[1] | wordr2[1];
    w2[2] = wordl2[2] | wordr2[2];
    w2[3] = wordl2[3] | wordr2[3];
    w3[0] = wordl3[0] | wordr3[0];
    w3[1] = wordl3[1] | wordr3[1];
    w3[2] = wordl3[2] | wordr3[2];
    w3[3] = wordl3[3] | wordr3[3];

    w0[0] = swap32 (w0[0]);
    w0[1] = swap32 (w0[1]);
    w0[2] = swap32 (w0[2]);
    w0[3] = swap32 (w0[3]);
    w1[0] = swap32 (w1[0]);
    w1[1] = swap32 (w1[1]);
    w1[2] = swap32 (w1[2]);
    w1[3] = swap32 (w1[3]);
    w2[0] = swap32 (w2[0]);
    w2[1] = swap32 (w2[1]);
    w2[2] = swap32 (w2[2]);
    w2[3] = swap32 (w2[3]);
    w3[0] = swap32 (w3[0]);
    w3[1] = swap32 (w3[1]);
    w3[2] = swap32 (w3[2]);
    w3[3] = swap32 (w3[3]);


    sha1_ctx_t sha1_ctx;
    sha1_init(&sha1_ctx);

    //printf("key?: 0x%x, 0x%x, 0x%x, 0x%x\n", swap32(w0[0]), swap32(w0[1]), w0[2], w0[3]);
    //printf("saltLen: 0x%x\n", salt_len);
    //printf("salt?: 0x%x, 0x%x, 0x%x, 0x%x\n", s_salt_buf[0].salt_buf[0], s_salt_buf[0].salt_buf[1], s_salt_buf[0].salt_buf[2], s_salt_buf[0].salt_buf[3]);
    //printf("salt?: 0x%x, 0x%x, 0x%x, 0x%x\n", salt_bufs[salt_pos].salt_buf[0], salt_bufs[salt_pos].salt_buf[1], salt_bufs[salt_pos].salt_buf[2], salt_bufs[salt_pos].salt_buf[3]);

    u32x ipadData[16];
    ipadData[0] = w0[0] ^ 0x36363636;
    ipadData[1] = w0[1] ^ 0x36363636;
    ipadData[2] = w0[2] ^ 0x36363636;
    ipadData[3] = w0[3] ^ 0x36363636;
    ipadData[4] = w1[0] ^ 0x36363636;
    ipadData[5] = w1[1] ^ 0x36363636;
    ipadData[6] = w1[2] ^ 0x36363636;
    ipadData[7] = w1[3] ^ 0x36363636;
    ipadData[8] = w2[0] ^ 0x36363636;
    ipadData[9] = w2[1] ^ 0x36363636;
    ipadData[10] = w2[2] ^ 0x36363636;
    ipadData[11] = w2[3] ^ 0x36363636;
    ipadData[12] = w3[0] ^ 0x36363636;
    ipadData[13] = w3[1] ^ 0x36363636;
    ipadData[14] = w3[2] ^ 0x36363636;
    ipadData[15] = w3[3] ^ 0x36363636;

    u32x opadData[16];
    opadData[0] = w0[0] ^ 0x5c5c5c5c;
    opadData[1] = w0[1] ^ 0x5c5c5c5c;
    opadData[2] = w0[2] ^ 0x5c5c5c5c;
    opadData[3] = w0[3] ^ 0x5c5c5c5c;
    opadData[4] = w1[0] ^ 0x5c5c5c5c;
    opadData[5] = w1[1] ^ 0x5c5c5c5c;
    opadData[6] = w1[2] ^ 0x5c5c5c5c;
    opadData[7] = w1[3] ^ 0x5c5c5c5c;
    opadData[8] = w2[0] ^ 0x5c5c5c5c;
    opadData[9] = w2[1] ^ 0x5c5c5c5c;
    opadData[10] = w2[2] ^ 0x5c5c5c5c;
    opadData[11] = w2[3] ^ 0x5c5c5c5c;
    opadData[12] = w3[0] ^ 0x5c5c5c5c;
    opadData[13] = w3[1] ^ 0x5c5c5c5c;
    opadData[14] = w3[2] ^ 0x5c5c5c5c;
    opadData[15] = w2[3] ^ 0x5c5c5c5c;

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

    u32x tmpData[16];
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


    COMPARE_S_SIMD (opad_ctx.h[3], opad_ctx.h[4], opad_ctx.h[2], opad_ctx.h[1]);
  }
}

__kernel void m00150_s08 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m00150_s16 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}
