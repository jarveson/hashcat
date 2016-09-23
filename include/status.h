/**
 * Author......: See docs/credits.txt
 * License.....: MIT
 */

#ifndef _STATUS_H
#define _STATUS_H

#include <stdio.h>
#include <inttypes.h>

double get_avg_exec_time (hc_device_param_t *device_param, const int last_num_entries);

void status_display_machine_readable (opencl_ctx_t *opencl_ctx, const hashes_t *hashes, const user_options_t *user_options, const user_options_extra_t *user_options_extra);
void status_display (opencl_ctx_t *opencl_ctx, const hashconfig_t *hashconfig, const hashes_t *hashes, const user_options_t *user_options, const user_options_extra_t *user_options_extra);
void status_benchmark_automate (opencl_ctx_t *opencl_ctx, const hashconfig_t *hashconfig);
void status_benchmark (opencl_ctx_t *opencl_ctx, const hashconfig_t *hashconfig, const user_options_t *user_options);

#endif // _STATUS_H