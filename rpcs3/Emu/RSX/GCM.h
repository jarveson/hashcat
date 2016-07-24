#pragma once

#include "Emu/Memory/vm.h"

enum
{
	CELL_GCM_DISPLAY_HSYNC            = 1,
	CELL_GCM_DISPLAY_VSYNC            = 2,
	CELL_GCM_DISPLAY_HSYNC_WITH_NOISE = 3,
};

enum
{
	CELL_GCM_DEBUG_LEVEL0 = 1,
	CELL_GCM_DEBUG_LEVEL1 = 2,
	CELL_GCM_DEBUG_LEVEL2 = 3,
};

enum
{
	CELL_GCM_DISPLAY_FREQUENCY_59_94HZ = 1,
	CELL_GCM_DISPLAY_FREQUENCY_SCANOUT = 2,
	CELL_GCM_DISPLAY_FREQUENCY_DISABLE = 3,
};

namespace rsx
{
	enum class vertex_base_type : u8
	{
		s1, ///< signed byte
		f, ///< float
		sf, ///< half float
		ub, ///< unsigned byte interpreted as 0.f and 1.f
		s32k, ///< signed 32bits int
		cmp, ///< compressed aka X11G11Z10 and always 1. W.
		ub256, ///< unsigned byte interpreted as between 0 and 255.
	};

	vertex_base_type to_vertex_base_type(u8 in);

	enum class index_array_type : u8
	{
		u32,
		u16,
	};

	index_array_type to_index_array_type(u8 in);

	enum class primitive_type : u8
	{
		invalid,
		points,
		lines,
		line_loop, // line strip with last end being joined with first end.
		line_strip,
		triangles,
		triangle_strip,
		triangle_fan, // like strip except that every triangle share the first vertex and one instead of 2 from previous triangle.
		quads,
		quad_strip,
		polygon, // convex polygon
	};

	primitive_type to_primitive_type(u8 in);

	enum class surface_target : u8
	{
		none,
		surface_a,
		surface_b,
		surfaces_a_b,
		surfaces_a_b_c,
		surfaces_a_b_c_d,
	};

	surface_target to_surface_target(u8 in);

	enum class surface_depth_format : u8
	{
		z16, // unsigned 16 bits depth
		z24s8, // unsigned 24 bits depth + 8 bits stencil
	};

	surface_depth_format to_surface_depth_format(u8 in);

	enum class surface_antialiasing : u8
	{
		center_1_sample,
		diagonal_centered_2_samples,
		square_centered_4_samples,
		square_rotated_4_samples,
	};

	surface_antialiasing to_surface_antialiasing(u8 in);

	enum class surface_color_format : u8
	{
		x1r5g5b5_z1r5g5b5,
		x1r5g5b5_o1r5g5b5,
		r5g6b5,
		x8r8g8b8_z8r8g8b8,
		x8r8g8b8_o8r8g8b8,
		a8r8g8b8,
		b8,
		g8b8,
		w16z16y16x16,
		w32z32y32x32,
		x32,
		x8b8g8r8_z8b8g8r8,
		x8b8g8r8_o8b8g8r8,
		a8b8g8r8,
	};

	surface_color_format to_surface_color_format(u8 in);

	enum class window_origin : u8
	{
		top,
		bottom
	};

	window_origin to_window_origin(u8 in);

	enum class window_pixel_center : u8
	{
		half,
		integer
	};

	window_pixel_center to_window_pixel_center(u8 in);

	enum class comparison_function : u8
	{
		never,
		less,
		equal,
		less_or_equal,
		greater,
		not_equal,
		greater_or_equal,
		always
	};

	comparison_function to_comparison_function(u16 in);

	enum class fog_mode : u8
	{
		linear,
		exponential,
		exponential2,
		exponential_abs,
		exponential2_abs,
		linear_abs
	};

	fog_mode to_fog_mode(u32 in);

	enum class texture_dimension : u8
	{
		dimension1d,
		dimension2d,
		dimension3d,
	};

	texture_dimension to_texture_dimension(u8 in);

	enum class texture_wrap_mode : u8
	{
		wrap,
		mirror,
		clamp_to_edge,
		border,
		clamp,
		mirror_once_clamp_to_edge,
		mirror_once_border,
		mirror_once_clamp,
	};

	texture_wrap_mode to_texture_wrap_mode(u8 in);

	enum class texture_max_anisotropy : u8
	{
		x1,
		x2,
		x4,
		x6,
		x8,
		x10,
		x12,
		x16,
	};

	texture_max_anisotropy to_texture_max_anisotropy(u8 in);

	enum class texture_minify_filter : u8
	{
		nearest, ///< no filtering, mipmap base level
		linear, ///< linear filtering, mipmap base level
		nearest_nearest, ///< no filtering, closest mipmap level
		linear_nearest, ///< linear filtering, closest mipmap level
		nearest_linear, ///< no filtering, linear mix between closest mipmap levels
		linear_linear, ///< linear filtering, linear mix between closest mipmap levels
		convolution_min, ///< Unknown mode but looks close to linear_linear
	};

	texture_minify_filter to_texture_minify_filter(u8 in);

	enum class texture_magnify_filter : u8
	{
		nearest, ///< no filtering
		linear, ///< linear filtering
		convolution_mag, ///< Unknown mode but looks close to linear
	};

	texture_magnify_filter to_texture_magnify_filter(u8 in);

	enum class stencil_op : u8
	{
		keep,
		zero,
		replace,
		incr,
		decr,
		invert,
		incr_wrap,
		decr_wrap,
	};

	stencil_op to_stencil_op(u16 in);

	enum class blend_equation : u8
	{
		add,
		min,
		max,
		substract,
		reverse_substract,
		reverse_substract_signed,
		add_signed,
		reverse_add_signed,
	};

	blend_equation to_blend_equation(u16 in);

	enum class blend_factor : u8
	{
		zero,
		one,
		src_color,
		one_minus_src_color,
		dst_color,
		one_minus_dst_color,
		src_alpha,
		one_minus_src_alpha,
		dst_alpha,
		one_minus_dst_alpha,
		src_alpha_saturate,
		constant_color,
		one_minus_constant_color,
		constant_alpha,
		one_minus_constant_alpha,
	};

	blend_factor to_blend_factor(u16 in);

	enum class logic_op : u8
	{
		logic_clear,
		logic_and,
		logic_and_reverse,
		logic_copy,
		logic_and_inverted,
		logic_noop,
		logic_xor,
		logic_or,
		logic_nor,
		logic_equiv,
		logic_invert,
		logic_or_reverse,
		logic_copy_inverted,
		logic_or_inverted,
		logic_nand,
		logic_set,
	};

	logic_op to_logic_op(u16 in);

	enum class front_face : u8
	{
		cw, /// clockwise
		ccw /// counter clockwise
	};

	front_face to_front_face(u16 in);

	enum class cull_face : u8
	{
		front,
		back,
		front_and_back,
	};

	cull_face to_cull_face(u16 in);

	enum class user_clip_plane_op : u8
	{
		disable,
		less_than,
		greather_or_equal,
	};

	user_clip_plane_op to_user_clip_plane_op(u8 in);

	enum class shading_mode : u8
	{
		smooth,
		flat,
	};

	shading_mode to_shading_mode(u32 in);

	enum class polygon_mode : u8
	{
		point,
		line,
		fill,
	};

	polygon_mode to_polygon_mode(u32 in);

	namespace blit_engine
	{
		enum class transfer_origin : u8
		{
			center,
			corner,
		};

		transfer_origin to_transfer_origin(u8 in);

		enum class transfer_interpolator : u8
		{
			zoh,
			foh,
		};

		transfer_interpolator to_transfer_interpolator(u8 in);

		enum class transfer_operation : u8
		{
			srccopy_and,
			rop_and,
			blend_and,
			srccopy,
			srccopy_premult,
			blend_premult,
		};

		transfer_operation to_transfer_operation(u8 in);

		enum class transfer_source_format : u8
		{
			a1r5g5b5,
			x1r5g5b5,
			a8r8g8b8,
			x8r8g8b8,
			cr8yb8cb8ya8,
			yb8cr8ya8cb8,
			r5g6b5,
			y8,
			ay8,
			eyb8ecr8eya8ecb8,
			ecr8eyb8ecb8eya8,
			a8b8g8r8,
			x8b8g8r8,
		};

		transfer_source_format to_transfer_source_format(u8 in);

		enum class transfer_destination_format : u8
		{
			r5g6b5,
			a8r8g8b8,
			y32,
		};

		transfer_destination_format to_transfer_destination_format(u8 in);

		enum class context_surface : u8
		{
			surface2d,
			swizzle2d,
		};

		context_surface to_context_surface(u32 in);

		enum class context_dma : u8
		{
			to_memory_get_report,
			report_location_main,
		};

		context_dma to_context_dma(u32 in);
	}
}

enum
{
	CELL_GCM_DISPLAY_FLIP_STATUS_        = 0,
	CELL_GCM_DISPLAY_FLIP_STATUS_WAITING = 1,
};

enum
{
	CELL_GCM_LOCATION_LOCAL = 0,
	CELL_GCM_LOCATION_MAIN  = 1,
};

enum
{
	CELL_GCM_FREQUENCY_MODULO = 1,
	CELL_GCM_FREQUENCY_DIVIDE = 0,
};

enum CellRescTableElement
{
	CELL_RESC_ELEMENT_HALF  = 0,
	CELL_RESC_ELEMENT_FLOAT = 1,
};

enum
{
	CELL_GCM_SYSTEM_MODE_IOMAP_512MB = 1,
};

// GCM Texture
enum
{
	// Color Flag
	CELL_GCM_TEXTURE_B8                     = 0x81,
	CELL_GCM_TEXTURE_A1R5G5B5               = 0x82,
	CELL_GCM_TEXTURE_A4R4G4B4               = 0x83,
	CELL_GCM_TEXTURE_R5G6B5                 = 0x84,
	CELL_GCM_TEXTURE_A8R8G8B8               = 0x85,
	CELL_GCM_TEXTURE_COMPRESSED_DXT1        = 0x86,
	CELL_GCM_TEXTURE_COMPRESSED_DXT23       = 0x87,
	CELL_GCM_TEXTURE_COMPRESSED_DXT45       = 0x88,
	CELL_GCM_TEXTURE_G8B8                   = 0x8B,
	CELL_GCM_TEXTURE_R6G5B5                 = 0x8F,
	CELL_GCM_TEXTURE_DEPTH24_D8             = 0x90,
	CELL_GCM_TEXTURE_DEPTH24_D8_FLOAT       = 0x91,
	CELL_GCM_TEXTURE_DEPTH16                = 0x92,
	CELL_GCM_TEXTURE_DEPTH16_FLOAT          = 0x93,
	CELL_GCM_TEXTURE_X16                    = 0x94,
	CELL_GCM_TEXTURE_Y16_X16                = 0x95,
	CELL_GCM_TEXTURE_R5G5B5A1               = 0x97,
	CELL_GCM_TEXTURE_COMPRESSED_HILO8       = 0x98,
	CELL_GCM_TEXTURE_COMPRESSED_HILO_S8     = 0x99,
	CELL_GCM_TEXTURE_W16_Z16_Y16_X16_FLOAT  = 0x9A,
	CELL_GCM_TEXTURE_W32_Z32_Y32_X32_FLOAT  = 0x9B,
	CELL_GCM_TEXTURE_X32_FLOAT              = 0x9C,
	CELL_GCM_TEXTURE_D1R5G5B5               = 0x9D,
	CELL_GCM_TEXTURE_D8R8G8B8               = 0x9E,
	CELL_GCM_TEXTURE_Y16_X16_FLOAT          = 0x9F,
	CELL_GCM_TEXTURE_COMPRESSED_B8R8_G8R8   = 0xAD,
	CELL_GCM_TEXTURE_COMPRESSED_R8B8_R8G8   = 0xAE,

	// Swizzle Flag
	CELL_GCM_TEXTURE_SZ = 0x00,
	CELL_GCM_TEXTURE_LN = 0x20,
	
	// Normalization Flag
	CELL_GCM_TEXTURE_NR = 0x00,
	CELL_GCM_TEXTURE_UN = 0x40,
};

// GCM Surface
enum
{
	// Surface type
	CELL_GCM_SURFACE_PITCH    = 1,
	CELL_GCM_SURFACE_SWIZZLE  = 2,
};

enum
{
	CELL_GCM_TEXTURE_UNSIGNED_REMAP_NORMAL = 0,
	CELL_GCM_TEXTURE_UNSIGNED_REMAP_BIASED = 1,

	CELL_GCM_TEXTURE_SIGNED_REMAP_NORMAL = 0x0,
	CELL_GCM_TEXTURE_SIGNED_REMAP_CLAMPED = 0x3,

	CELL_GCM_TEXTURE_ZFUNC_NEVER = 0,
	CELL_GCM_TEXTURE_ZFUNC_LESS = 1,
	CELL_GCM_TEXTURE_ZFUNC_EQUAL = 2,
	CELL_GCM_TEXTURE_ZFUNC_LEQUAL = 3,
	CELL_GCM_TEXTURE_ZFUNC_GREATER = 4,
	CELL_GCM_TEXTURE_ZFUNC_NOTEQUAL = 5,
	CELL_GCM_TEXTURE_ZFUNC_GEQUAL = 6,
	CELL_GCM_TEXTURE_ZFUNC_ALWAYS = 7,

	CELL_GCM_TEXTURE_GAMMA_R = 1 << 0,
	CELL_GCM_TEXTURE_GAMMA_G = 1 << 1,
	CELL_GCM_TEXTURE_GAMMA_B = 1 << 2,
	CELL_GCM_TEXTURE_GAMMA_A = 1 << 3,

	CELL_GCM_TEXTURE_ANISO_SPREAD_0_50_TEXEL = 0x0,
	CELL_GCM_TEXTURE_ANISO_SPREAD_1_00_TEXEL = 0x1,
	CELL_GCM_TEXTURE_ANISO_SPREAD_1_125_TEXEL = 0x2,
	CELL_GCM_TEXTURE_ANISO_SPREAD_1_25_TEXEL = 0x3,
	CELL_GCM_TEXTURE_ANISO_SPREAD_1_375_TEXEL = 0x4,
	CELL_GCM_TEXTURE_ANISO_SPREAD_1_50_TEXEL = 0x5,
	CELL_GCM_TEXTURE_ANISO_SPREAD_1_75_TEXEL = 0x6,
	CELL_GCM_TEXTURE_ANISO_SPREAD_2_00_TEXEL = 0x7,

	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX0_U = 1 << 0,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX0_V = 1 << 1,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX0_P = 1 << 2,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX0_Q = 1 << 3,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX1_U = 1 << 4,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX1_V = 1 << 5,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX1_P = 1 << 6,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX1_Q = 1 << 7,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX2_U = 1 << 8,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX2_V = 1 << 9,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX2_P = 1 << 10,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX2_Q = 1 << 11,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX3_U = 1 << 12,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX3_V = 1 << 13,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX3_P = 1 << 14,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX3_Q = 1 << 15,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX4_U = 1 << 16,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX4_V = 1 << 17,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX4_P = 1 << 18,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX4_Q = 1 << 19,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX5_U = 1 << 20,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX5_V = 1 << 21,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX5_P = 1 << 22,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX5_Q = 1 << 23,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX6_U = 1 << 24,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX6_V = 1 << 25,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX6_P = 1 << 26,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX6_Q = 1 << 27,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX7_U = 1 << 28,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX7_V = 1 << 29,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX7_P = 1 << 30,
	CELL_GCM_TEXTURE_CYLINDRICAL_WRAP_ENABLE_TEX7_Q = 1 << 31,

	CELL_GCM_COLOR_MASK_B = 1 << 0,
	CELL_GCM_COLOR_MASK_G = 1 << 8,
	CELL_GCM_COLOR_MASK_R = 1 << 16,
	CELL_GCM_COLOR_MASK_A = 1 << 24,

	CELL_GCM_COLOR_MASK_MRT1_A = 1 << 4,
	CELL_GCM_COLOR_MASK_MRT1_R = 1 << 5,
	CELL_GCM_COLOR_MASK_MRT1_G = 1 << 6,
	CELL_GCM_COLOR_MASK_MRT1_B = 1 << 7,
	CELL_GCM_COLOR_MASK_MRT2_A = 1 << 8,
	CELL_GCM_COLOR_MASK_MRT2_R = 1 << 9,
	CELL_GCM_COLOR_MASK_MRT2_G = 1 << 10,
	CELL_GCM_COLOR_MASK_MRT2_B = 1 << 11,
	CELL_GCM_COLOR_MASK_MRT3_A = 1 << 12,
	CELL_GCM_COLOR_MASK_MRT3_R = 1 << 13,
	CELL_GCM_COLOR_MASK_MRT3_G = 1 << 14,
	CELL_GCM_COLOR_MASK_MRT3_B = 1 << 15,

	CELL_GCM_NEVER = 0x0200,
	CELL_GCM_LESS = 0x0201,
	CELL_GCM_EQUAL = 0x0202,
	CELL_GCM_LEQUAL = 0x0203,
	CELL_GCM_GREATER = 0x0204,
	CELL_GCM_NOTEQUAL = 0x0205,
	CELL_GCM_GEQUAL = 0x0206,
	CELL_GCM_ALWAYS = 0x0207,

	CELL_GCM_ZERO = 0,
	CELL_GCM_ONE = 1,

	CELL_GCM_FRONT = 0x0404,
	CELL_GCM_BACK = 0x0405,
	CELL_GCM_FRONT_AND_BACK = 0x0408,

	CELL_GCM_CW = 0x0900,
	CELL_GCM_CCW = 0x0901,

	CELL_GCM_TRANSFER_LOCAL_TO_LOCAL = 0,
	CELL_GCM_TRANSFER_MAIN_TO_LOCAL = 1,
	CELL_GCM_TRANSFER_LOCAL_TO_MAIN = 2,
	CELL_GCM_TRANSFER_MAIN_TO_MAIN = 3,

	CELL_GCM_INVALIDATE_TEXTURE = 1,
	CELL_GCM_INVALIDATE_VERTEX_TEXTURE = 2,

	CELL_GCM_COMPMODE_DISABLED = 0,
	CELL_GCM_COMPMODE_C32_2X1 = 7,
	CELL_GCM_COMPMODE_C32_2X2 = 8,
	CELL_GCM_COMPMODE_Z32_SEPSTENCIL = 9,
	CELL_GCM_COMPMODE_Z32_SEPSTENCIL_REG = 10,
	CELL_GCM_COMPMODE_Z32_SEPSTENCIL_REGULAR = 10,
	CELL_GCM_COMPMODE_Z32_SEPSTENCIL_DIAGONAL = 11,
	CELL_GCM_COMPMODE_Z32_SEPSTENCIL_ROTATED = 12,

	CELL_GCM_ZCULL_Z16 = 1,
	CELL_GCM_ZCULL_Z24S8 = 2,
	CELL_GCM_ZCULL_MSB = 0,
	CELL_GCM_ZCULL_LONES = 1,
	CELL_GCM_ZCULL_LESS = 0,
	CELL_GCM_ZCULL_GREATER = 1,

	CELL_GCM_SCULL_SFUNC_NEVER = 0,
	CELL_GCM_SCULL_SFUNC_LESS = 1,
	CELL_GCM_SCULL_SFUNC_EQUAL = 2,
	CELL_GCM_SCULL_SFUNC_LEQUAL = 3,
	CELL_GCM_SCULL_SFUNC_GREATER = 4,
	CELL_GCM_SCULL_SFUNC_NOTEQUAL = 5,
	CELL_GCM_SCULL_SFUNC_GEQUAL = 6,
	CELL_GCM_SCULL_SFUNC_ALWAYS = 7,

	CELL_GCM_ATTRIB_OUTPUT_FRONTDIFFUSE = 0,
	CELL_GCM_ATTRIB_OUTPUT_FRONTSPECULAR = 1,
	CELL_GCM_ATTRIB_OUTPUT_BACKDIFFUSE = 2,
	CELL_GCM_ATTRIB_OUTPUT_BACKSPECULAR = 3,
	CELL_GCM_ATTRIB_OUTPUT_FOG = 4,
	CELL_GCM_ATTRIB_OUTPUT_POINTSIZE = 5,
	CELL_GCM_ATTRIB_OUTPUT_UC0 = 6,
	CELL_GCM_ATTRIB_OUTPUT_UC1 = 7,
	CELL_GCM_ATTRIB_OUTPUT_UC2 = 8,
	CELL_GCM_ATTRIB_OUTPUT_UC3 = 9,
	CELL_GCM_ATTRIB_OUTPUT_UC4 = 10,
	CELL_GCM_ATTRIB_OUTPUT_UC5 = 11,
	CELL_GCM_ATTRIB_OUTPUT_TEX8 = 12,
	CELL_GCM_ATTRIB_OUTPUT_TEX9 = 13,
	CELL_GCM_ATTRIB_OUTPUT_TEX0 = 14,
	CELL_GCM_ATTRIB_OUTPUT_TEX1 = 15,
	CELL_GCM_ATTRIB_OUTPUT_TEX2 = 16,
	CELL_GCM_ATTRIB_OUTPUT_TEX3 = 17,
	CELL_GCM_ATTRIB_OUTPUT_TEX4 = 18,
	CELL_GCM_ATTRIB_OUTPUT_TEX5 = 19,
	CELL_GCM_ATTRIB_OUTPUT_TEX6 = 20,
	CELL_GCM_ATTRIB_OUTPUT_TEX7 = 21,

	CELL_GCM_ATTRIB_OUTPUT_MASK_FRONTDIFFUSE = 1 << CELL_GCM_ATTRIB_OUTPUT_FRONTDIFFUSE,
	CELL_GCM_ATTRIB_OUTPUT_MASK_FRONTSPECULAR = 1 << CELL_GCM_ATTRIB_OUTPUT_FRONTSPECULAR,
	CELL_GCM_ATTRIB_OUTPUT_MASK_BACKDIFFUSE = 1 << CELL_GCM_ATTRIB_OUTPUT_BACKDIFFUSE,
	CELL_GCM_ATTRIB_OUTPUT_MASK_BACKSPECULAR = 1 << CELL_GCM_ATTRIB_OUTPUT_BACKSPECULAR,
	CELL_GCM_ATTRIB_OUTPUT_MASK_FOG = 1 << CELL_GCM_ATTRIB_OUTPUT_FOG,
	CELL_GCM_ATTRIB_OUTPUT_MASK_POINTSIZE = 1 << CELL_GCM_ATTRIB_OUTPUT_POINTSIZE,
	CELL_GCM_ATTRIB_OUTPUT_MASK_UC0 = 1 << CELL_GCM_ATTRIB_OUTPUT_UC0,
	CELL_GCM_ATTRIB_OUTPUT_MASK_UC1 = 1 << CELL_GCM_ATTRIB_OUTPUT_UC1,
	CELL_GCM_ATTRIB_OUTPUT_MASK_UC2 = 1 << CELL_GCM_ATTRIB_OUTPUT_UC2,
	CELL_GCM_ATTRIB_OUTPUT_MASK_UC3 = 1 << CELL_GCM_ATTRIB_OUTPUT_UC3,
	CELL_GCM_ATTRIB_OUTPUT_MASK_UC4 = 1 << CELL_GCM_ATTRIB_OUTPUT_UC4,
	CELL_GCM_ATTRIB_OUTPUT_MASK_UC5 = 1 << CELL_GCM_ATTRIB_OUTPUT_UC5,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX8 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX8,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX9 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX9,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX0 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX0,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX1 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX1,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX2 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX2,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX3 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX3,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX4 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX4,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX5 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX5,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX6 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX6,
	CELL_GCM_ATTRIB_OUTPUT_MASK_TEX7 = 1 << CELL_GCM_ATTRIB_OUTPUT_TEX7,

	CELL_GCM_TRUE = 1,
	CELL_GCM_FALSE = 0,
};

enum
{
	CELL_GCM_SHADER_CONTROL_DEPTH_EXPORT = 0xe, ///< shader program exports the depth of the shaded fragment
	CELL_GCM_SHADER_CONTROL_32_BITS_EXPORTS = 0x40 ///< shader program exports 32 bits registers values (instead of 16 bits ones)
};

// GCM Reports
enum
{
	CELL_GCM_ZPASS_PIXEL_CNT   = 1,
	CELL_GCM_ZCULL_STATS       = 2,
	CELL_GCM_ZCULL_STATS1      = 3,
	CELL_GCM_ZCULL_STATS2      = 4,
	CELL_GCM_ZCULL_STATS3      = 5,
};

// GPU Class Handles
enum
{
	CELL_GCM_CONTEXT_DMA_MEMORY_FRAME_BUFFER   = 0xFEED0000, // Local memory
	CELL_GCM_CONTEXT_DMA_MEMORY_HOST_BUFFER    = 0xFEED0001, // Main memory
	CELL_GCM_CONTEXT_DMA_TO_MEMORY_GET_REPORT = 0x66626660,
	CELL_GCM_CONTEXT_DMA_REPORT_LOCATION_MAIN = 0xBAD68000,
	CELL_GCM_CONTEXT_DMA_NOTIFY_MAIN_0         = 0x6660420F,
};

struct CellGcmControl
{
	atomic_be_t<u32> put;
	atomic_be_t<u32> get;
	atomic_be_t<u32> ref;
};

struct CellGcmConfig
{
	be_t<u32> localAddress;
	be_t<u32> ioAddress;
	be_t<u32> localSize;
	be_t<u32> ioSize;
	be_t<u32> memoryFrequency;
	be_t<u32> coreFrequency;
};

struct CellGcmContextData;

typedef s32(CellGcmContextCallback)(vm::ps3::ptr<CellGcmContextData>, u32);

struct CellGcmContextData
{
	vm::ps3::bptr<u32> begin;
	vm::ps3::bptr<u32> end;
	vm::ps3::bptr<u32> current;
	vm::ps3::bptr<CellGcmContextCallback> callback;
};

struct gcmInfo
{
	u32 config_addr;
	u32 context_addr;
	u32 control_addr;
	u32 label_addr;
};

struct CellGcmSurface
{
	u8 type;
	u8 antialias;
	u8 colorFormat;
	u8 colorTarget;
	u8 colorLocation[4];
	be_t<u32> colorOffset[4];
	be_t<u32> colorPitch[4];
	u8 depthFormat;
	u8 depthLocation;
	u8 _padding[2];
	be_t<u32> depthOffset;
	be_t<u32> depthPitch;
	be_t<u16> width;
	be_t<u16> height;
	be_t<u16> x;
	be_t<u16> y;
};

struct CellGcmReportData
{
	be_t<u64> timer;
	be_t<u32> value;
	be_t<u32> padding;
};

struct CellGcmZcullInfo
{
	u32 region;
	u32 size;
	u32 start;
	u32 offset;
	u32 status0;
	u32 status1;
};

struct CellGcmDisplayInfo
{
	be_t<u32> offset;
	be_t<u32> pitch;
	be_t<u32> width;
	be_t<u32> height;
};

struct CellGcmTileInfo
{
	be_t<u32> tile;
	be_t<u32> limit;
	be_t<u32> pitch;
	be_t<u32> format;
};

struct GcmZcullInfo
{
	u32 offset;
	u32 width;
	u32 height;
	u32 cullStart;
	u32 zFormat;
	u32 aaFormat;
	u32 zcullDir;
	u32 zcullFormat;
	u32 sFunc;
	u32 sRef;
	u32 sMask;
	bool binded;

	GcmZcullInfo()
	{
		memset(this, 0, sizeof(*this));
	}

	CellGcmZcullInfo pack() const
	{
		CellGcmZcullInfo ret;

		ret.region = (1<<0) | (zFormat<<4) | (aaFormat<<8);
		ret.size = ((width>>6)<<22) | ((height>>6)<<6);
		ret.start = cullStart&(~0xFFF);
		ret.offset = offset;
		ret.status0 = (zcullDir<<1) | (zcullFormat<<2) | ((sFunc&0xF)<<12) | (sRef<<16) | (sMask<<24);
		ret.status1 = (0x2000<<0) | (0x20<<16);

		return ret;
	}
};

struct GcmTileInfo
{
	u32 location;
	u32 offset;
	u32 size;
	u32 pitch;
	u32 comp;
	u32 base;
	u32 bank;
	bool binded;

	GcmTileInfo()
	{
		memset(this, 0, sizeof(*this));
	}

	CellGcmTileInfo pack() const
	{
		CellGcmTileInfo ret;

		ret.tile = (location + 1) | (bank << 4) | ((offset / 0x10000) << 16) | (location << 31);
		ret.limit = ((offset + size - 1) / 0x10000) << 16 | (location << 31);
		ret.pitch = (pitch / 0x100) << 8;
		ret.format = base | ((base + ((size - 1) / 0x10000)) << 13) | (comp << 26) | (1 << 30);

		return ret;
	}
};

enum
{
	// NV40_CHANNEL_DMA (NV406E)
	NV406E_SET_REFERENCE                     = 0x00000050 >> 2,
	NV406E_SET_CONTEXT_DMA_SEMAPHORE         = 0x00000060 >> 2,
	NV406E_SEMAPHORE_OFFSET                  = 0x00000064 >> 2,
	NV406E_SEMAPHORE_ACQUIRE                 = 0x00000068 >> 2,
	NV406E_SEMAPHORE_RELEASE                 = 0x0000006c >> 2,

	// NV40_CURIE_PRIMITIVE	(NV4097)
	NV4097_SET_OBJECT                        = 0x00000000 >> 2,
	NV4097_NO_OPERATION                      = 0x00000100 >> 2,
	NV4097_NOTIFY                            = 0x00000104 >> 2,
	NV4097_WAIT_FOR_IDLE                     = 0x00000110 >> 2,
	NV4097_PM_TRIGGER                        = 0x00000140 >> 2,
	NV4097_SET_CONTEXT_DMA_NOTIFIES          = 0x00000180 >> 2,
	NV4097_SET_CONTEXT_DMA_A                 = 0x00000184 >> 2,
	NV4097_SET_CONTEXT_DMA_B                 = 0x00000188 >> 2,
	NV4097_SET_CONTEXT_DMA_COLOR_B           = 0x0000018c >> 2,
	NV4097_SET_CONTEXT_DMA_STATE             = 0x00000190 >> 2,
	NV4097_SET_CONTEXT_DMA_COLOR_A           = 0x00000194 >> 2,
	NV4097_SET_CONTEXT_DMA_ZETA              = 0x00000198 >> 2,
	NV4097_SET_CONTEXT_DMA_VERTEX_A          = 0x0000019c >> 2,
	NV4097_SET_CONTEXT_DMA_VERTEX_B          = 0x000001a0 >> 2,
	NV4097_SET_CONTEXT_DMA_SEMAPHORE         = 0x000001a4 >> 2,
	NV4097_SET_CONTEXT_DMA_REPORT            = 0x000001a8 >> 2,
	NV4097_SET_CONTEXT_DMA_CLIP_ID           = 0x000001ac >> 2,
	NV4097_SET_CONTEXT_DMA_CULL_DATA         = 0x000001b0 >> 2,
	NV4097_SET_CONTEXT_DMA_COLOR_C           = 0x000001b4 >> 2,
	NV4097_SET_CONTEXT_DMA_COLOR_D           = 0x000001b8 >> 2,
	NV4097_SET_SURFACE_CLIP_HORIZONTAL       = 0x00000200 >> 2,
	NV4097_SET_SURFACE_CLIP_VERTICAL         = 0x00000204 >> 2,
	NV4097_SET_SURFACE_FORMAT                = 0x00000208 >> 2,
	NV4097_SET_SURFACE_PITCH_A               = 0x0000020c >> 2,
	NV4097_SET_SURFACE_COLOR_AOFFSET         = 0x00000210 >> 2,
	NV4097_SET_SURFACE_ZETA_OFFSET           = 0x00000214 >> 2,
	NV4097_SET_SURFACE_COLOR_BOFFSET         = 0x00000218 >> 2,
	NV4097_SET_SURFACE_PITCH_B               = 0x0000021c >> 2,
	NV4097_SET_SURFACE_COLOR_TARGET          = 0x00000220 >> 2,
	NV4097_SET_SURFACE_PITCH_Z               = 0x0000022c >> 2,
	NV4097_INVALIDATE_ZCULL                  = 0x00000234 >> 2,
	NV4097_SET_CYLINDRICAL_WRAP              = 0x00000238 >> 2,
	NV4097_SET_CYLINDRICAL_WRAP1             = 0x0000023c >> 2,
	NV4097_SET_SURFACE_PITCH_C               = 0x00000280 >> 2,
	NV4097_SET_SURFACE_PITCH_D               = 0x00000284 >> 2,
	NV4097_SET_SURFACE_COLOR_COFFSET         = 0x00000288 >> 2,
	NV4097_SET_SURFACE_COLOR_DOFFSET         = 0x0000028c >> 2,
	NV4097_SET_WINDOW_OFFSET                 = 0x000002b8 >> 2,
	NV4097_SET_WINDOW_CLIP_TYPE              = 0x000002bc >> 2,
	NV4097_SET_WINDOW_CLIP_HORIZONTAL        = 0x000002c0 >> 2,
	NV4097_SET_WINDOW_CLIP_VERTICAL          = 0x000002c4 >> 2,
	NV4097_SET_DITHER_ENABLE                 = 0x00000300 >> 2,
	NV4097_SET_ALPHA_TEST_ENABLE             = 0x00000304 >> 2,
	NV4097_SET_ALPHA_FUNC                    = 0x00000308 >> 2,
	NV4097_SET_ALPHA_REF                     = 0x0000030c >> 2,
	NV4097_SET_BLEND_ENABLE                  = 0x00000310 >> 2,
	NV4097_SET_BLEND_FUNC_SFACTOR            = 0x00000314 >> 2,
	NV4097_SET_BLEND_FUNC_DFACTOR            = 0x00000318 >> 2,
	NV4097_SET_BLEND_COLOR                   = 0x0000031c >> 2,
	NV4097_SET_BLEND_EQUATION                = 0x00000320 >> 2,
	NV4097_SET_COLOR_MASK                    = 0x00000324 >> 2,
	NV4097_SET_STENCIL_TEST_ENABLE           = 0x00000328 >> 2,
	NV4097_SET_STENCIL_MASK                  = 0x0000032c >> 2,
	NV4097_SET_STENCIL_FUNC                  = 0x00000330 >> 2,
	NV4097_SET_STENCIL_FUNC_REF              = 0x00000334 >> 2,
	NV4097_SET_STENCIL_FUNC_MASK             = 0x00000338 >> 2,
	NV4097_SET_STENCIL_OP_FAIL               = 0x0000033c >> 2,
	NV4097_SET_STENCIL_OP_ZFAIL              = 0x00000340 >> 2,
	NV4097_SET_STENCIL_OP_ZPASS              = 0x00000344 >> 2,
	NV4097_SET_TWO_SIDED_STENCIL_TEST_ENABLE = 0x00000348 >> 2,
	NV4097_SET_BACK_STENCIL_MASK             = 0x0000034c >> 2,
	NV4097_SET_BACK_STENCIL_FUNC             = 0x00000350 >> 2,
	NV4097_SET_BACK_STENCIL_FUNC_REF         = 0x00000354 >> 2,
	NV4097_SET_BACK_STENCIL_FUNC_MASK        = 0x00000358 >> 2,
	NV4097_SET_BACK_STENCIL_OP_FAIL          = 0x0000035c >> 2,
	NV4097_SET_BACK_STENCIL_OP_ZFAIL         = 0x00000360 >> 2,
	NV4097_SET_BACK_STENCIL_OP_ZPASS         = 0x00000364 >> 2,
	NV4097_SET_SHADE_MODE                    = 0x00000368 >> 2,
	NV4097_SET_BLEND_ENABLE_MRT              = 0x0000036c >> 2,
	NV4097_SET_COLOR_MASK_MRT                = 0x00000370 >> 2,
	NV4097_SET_LOGIC_OP_ENABLE               = 0x00000374 >> 2,
	NV4097_SET_LOGIC_OP                      = 0x00000378 >> 2,
	NV4097_SET_BLEND_COLOR2                  = 0x0000037c >> 2,
	NV4097_SET_DEPTH_BOUNDS_TEST_ENABLE      = 0x00000380 >> 2,
	NV4097_SET_DEPTH_BOUNDS_MIN              = 0x00000384 >> 2,
	NV4097_SET_DEPTH_BOUNDS_MAX              = 0x00000388 >> 2,
	NV4097_SET_CLIP_MIN                      = 0x00000394 >> 2,
	NV4097_SET_CLIP_MAX                      = 0x00000398 >> 2,
	NV4097_SET_CONTROL0                      = 0x000003b0 >> 2,
	NV4097_SET_LINE_WIDTH                    = 0x000003b8 >> 2,
	NV4097_SET_LINE_SMOOTH_ENABLE            = 0x000003bc >> 2,
	NV4097_SET_ANISO_SPREAD                  = 0x000003c0 >> 2,
	NV4097_SET_SCISSOR_HORIZONTAL            = 0x000008c0 >> 2,
	NV4097_SET_SCISSOR_VERTICAL              = 0x000008c4 >> 2,
	NV4097_SET_FOG_MODE                      = 0x000008cc >> 2,
	NV4097_SET_FOG_PARAMS                    = 0x000008d0 >> 2,
	NV4097_SET_SHADER_PROGRAM                = 0x000008e4 >> 2,
	NV4097_SET_VERTEX_TEXTURE_OFFSET         = 0x00000900 >> 2,
	NV4097_SET_VERTEX_TEXTURE_FORMAT         = 0x00000904 >> 2,
	NV4097_SET_VERTEX_TEXTURE_ADDRESS        = 0x00000908 >> 2,
	NV4097_SET_VERTEX_TEXTURE_CONTROL0       = 0x0000090c >> 2,
	NV4097_SET_VERTEX_TEXTURE_CONTROL3       = 0x00000910 >> 2,
	NV4097_SET_VERTEX_TEXTURE_FILTER         = 0x00000914 >> 2,
	NV4097_SET_VERTEX_TEXTURE_IMAGE_RECT     = 0x00000918 >> 2,
	NV4097_SET_VERTEX_TEXTURE_BORDER_COLOR   = 0x0000091c >> 2,
	NV4097_SET_VIEWPORT_HORIZONTAL           = 0x00000a00 >> 2,
	NV4097_SET_VIEWPORT_VERTICAL             = 0x00000a04 >> 2,
	NV4097_SET_POINT_CENTER_MODE             = 0x00000a0c >> 2,
	NV4097_ZCULL_SYNC                        = 0x00000a1c >> 2,
	NV4097_SET_VIEWPORT_OFFSET               = 0x00000a20 >> 2,
	NV4097_SET_VIEWPORT_SCALE                = 0x00000a30 >> 2,
	NV4097_SET_POLY_OFFSET_POINT_ENABLE      = 0x00000a60 >> 2,
	NV4097_SET_POLY_OFFSET_LINE_ENABLE       = 0x00000a64 >> 2,
	NV4097_SET_POLY_OFFSET_FILL_ENABLE       = 0x00000a68 >> 2,
	NV4097_SET_DEPTH_FUNC                    = 0x00000a6c >> 2,
	NV4097_SET_DEPTH_MASK                    = 0x00000a70 >> 2,
	NV4097_SET_DEPTH_TEST_ENABLE             = 0x00000a74 >> 2,
	NV4097_SET_POLYGON_OFFSET_SCALE_FACTOR   = 0x00000a78 >> 2,
	NV4097_SET_POLYGON_OFFSET_BIAS           = 0x00000a7c >> 2,
	NV4097_SET_VERTEX_DATA_SCALED4S_M        = 0x00000a80 >> 2,
	NV4097_SET_TEXTURE_CONTROL2              = 0x00000b00 >> 2,
	NV4097_SET_TEX_COORD_CONTROL             = 0x00000b40 >> 2,
	NV4097_SET_TRANSFORM_PROGRAM             = 0x00000b80 >> 2,
	NV4097_SET_SPECULAR_ENABLE               = 0x00001428 >> 2,
	NV4097_SET_TWO_SIDE_LIGHT_EN             = 0x0000142c >> 2,
	NV4097_CLEAR_ZCULL_SURFACE               = 0x00001438 >> 2,
	NV4097_SET_PERFORMANCE_PARAMS            = 0x00001450 >> 2,
	NV4097_SET_FLAT_SHADE_OP                 = 0x00001454 >> 2,
	NV4097_SET_EDGE_FLAG                     = 0x0000145c >> 2,
	NV4097_SET_USER_CLIP_PLANE_CONTROL       = 0x00001478 >> 2,
	NV4097_SET_POLYGON_STIPPLE               = 0x0000147c >> 2,
	NV4097_SET_POLYGON_STIPPLE_PATTERN       = 0x00001480 >> 2,
	NV4097_SET_VERTEX_DATA3F_M               = 0x00001500 >> 2,
	NV4097_SET_VERTEX_DATA_ARRAY_OFFSET      = 0x00001680 >> 2,
	NV4097_INVALIDATE_VERTEX_CACHE_FILE      = 0x00001710 >> 2,
	NV4097_INVALIDATE_VERTEX_FILE            = 0x00001714 >> 2,
	NV4097_PIPE_NOP                          = 0x00001718 >> 2,
	NV4097_SET_VERTEX_DATA_BASE_OFFSET       = 0x00001738 >> 2,
	NV4097_SET_VERTEX_DATA_BASE_INDEX        = 0x0000173c >> 2,
	NV4097_SET_VERTEX_DATA_ARRAY_FORMAT      = 0x00001740 >> 2,
	NV4097_CLEAR_REPORT_VALUE                = 0x000017c8 >> 2,
	NV4097_SET_ZPASS_PIXEL_COUNT_ENABLE      = 0x000017cc >> 2,
	NV4097_GET_REPORT                        = 0x00001800 >> 2,
	NV4097_SET_ZCULL_STATS_ENABLE            = 0x00001804 >> 2,
	NV4097_SET_BEGIN_END                     = 0x00001808 >> 2,
	NV4097_ARRAY_ELEMENT16                   = 0x0000180c >> 2,
	NV4097_ARRAY_ELEMENT32                   = 0x00001810 >> 2,
	NV4097_DRAW_ARRAYS                       = 0x00001814 >> 2,
	NV4097_INLINE_ARRAY                      = 0x00001818 >> 2,
	NV4097_SET_INDEX_ARRAY_ADDRESS           = 0x0000181c >> 2,
	NV4097_SET_INDEX_ARRAY_DMA               = 0x00001820 >> 2,
	NV4097_DRAW_INDEX_ARRAY                  = 0x00001824 >> 2,
	NV4097_SET_FRONT_POLYGON_MODE            = 0x00001828 >> 2,
	NV4097_SET_BACK_POLYGON_MODE             = 0x0000182c >> 2,
	NV4097_SET_CULL_FACE                     = 0x00001830 >> 2,
	NV4097_SET_FRONT_FACE                    = 0x00001834 >> 2,
	NV4097_SET_POLY_SMOOTH_ENABLE            = 0x00001838 >> 2,
	NV4097_SET_CULL_FACE_ENABLE              = 0x0000183c >> 2,
	NV4097_SET_TEXTURE_CONTROL3              = 0x00001840 >> 2,
	NV4097_SET_VERTEX_DATA2F_M               = 0x00001880 >> 2,
	NV4097_SET_VERTEX_DATA2S_M               = 0x00001900 >> 2,
	NV4097_SET_VERTEX_DATA4UB_M              = 0x00001940 >> 2,
	NV4097_SET_VERTEX_DATA4S_M               = 0x00001980 >> 2,
	NV4097_SET_TEXTURE_OFFSET                = 0x00001a00 >> 2,
	NV4097_SET_TEXTURE_FORMAT                = 0x00001a04 >> 2,
	NV4097_SET_TEXTURE_ADDRESS               = 0x00001a08 >> 2,
	NV4097_SET_TEXTURE_CONTROL0              = 0x00001a0c >> 2,
	NV4097_SET_TEXTURE_CONTROL1              = 0x00001a10 >> 2,
	NV4097_SET_TEXTURE_FILTER                = 0x00001a14 >> 2,
	NV4097_SET_TEXTURE_IMAGE_RECT            = 0x00001a18 >> 2,
	NV4097_SET_TEXTURE_BORDER_COLOR          = 0x00001a1c >> 2,
	NV4097_SET_VERTEX_DATA4F_M               = 0x00001c00 >> 2,
	NV4097_SET_COLOR_KEY_COLOR               = 0x00001d00 >> 2,
	NV4097_SET_SHADER_CONTROL                = 0x00001d60 >> 2,
	NV4097_SET_INDEXED_CONSTANT_READ_LIMITS  = 0x00001d64 >> 2,
	NV4097_SET_SEMAPHORE_OFFSET              = 0x00001d6c >> 2,
	NV4097_BACK_END_WRITE_SEMAPHORE_RELEASE  = 0x00001d70 >> 2,
	NV4097_TEXTURE_READ_SEMAPHORE_RELEASE    = 0x00001d74 >> 2,
	NV4097_SET_ZMIN_MAX_CONTROL              = 0x00001d78 >> 2,
	NV4097_SET_ANTI_ALIASING_CONTROL         = 0x00001d7c >> 2,
	NV4097_SET_SURFACE_COMPRESSION           = 0x00001d80 >> 2,
	NV4097_SET_ZCULL_EN                      = 0x00001d84 >> 2,
	NV4097_SET_SHADER_WINDOW                 = 0x00001d88 >> 2,
	NV4097_SET_ZSTENCIL_CLEAR_VALUE          = 0x00001d8c >> 2,
	NV4097_SET_COLOR_CLEAR_VALUE             = 0x00001d90 >> 2,
	NV4097_CLEAR_SURFACE                     = 0x00001d94 >> 2,
	NV4097_SET_CLEAR_RECT_HORIZONTAL         = 0x00001d98 >> 2,
	NV4097_SET_CLEAR_RECT_VERTICAL           = 0x00001d9c >> 2,
	NV4097_SET_CLIP_ID_TEST_ENABLE           = 0x00001da4 >> 2,
	NV4097_SET_RESTART_INDEX_ENABLE          = 0x00001dac >> 2,
	NV4097_SET_RESTART_INDEX                 = 0x00001db0 >> 2,
	NV4097_SET_LINE_STIPPLE                  = 0x00001db4 >> 2,
	NV4097_SET_LINE_STIPPLE_PATTERN          = 0x00001db8 >> 2,
	NV4097_SET_VERTEX_DATA1F_M               = 0x00001e40 >> 2,
	NV4097_SET_TRANSFORM_EXECUTION_MODE      = 0x00001e94 >> 2,
	NV4097_SET_RENDER_ENABLE                 = 0x00001e98 >> 2,
	NV4097_SET_TRANSFORM_PROGRAM_LOAD        = 0x00001e9c >> 2,
	NV4097_SET_TRANSFORM_PROGRAM_START       = 0x00001ea0 >> 2,
	NV4097_SET_ZCULL_CONTROL0                = 0x00001ea4 >> 2,
	NV4097_SET_ZCULL_CONTROL1                = 0x00001ea8 >> 2,
	NV4097_SET_SCULL_CONTROL                 = 0x00001eac >> 2,
	NV4097_SET_POINT_SIZE                    = 0x00001ee0 >> 2,
	NV4097_SET_POINT_PARAMS_ENABLE           = 0x00001ee4 >> 2,
	NV4097_SET_POINT_SPRITE_CONTROL          = 0x00001ee8 >> 2,
	NV4097_SET_TRANSFORM_TIMEOUT             = 0x00001ef8 >> 2,
	NV4097_SET_TRANSFORM_CONSTANT_LOAD       = 0x00001efc >> 2,
	NV4097_SET_TRANSFORM_CONSTANT            = 0x00001f00 >> 2,
	NV4097_SET_FREQUENCY_DIVIDER_OPERATION   = 0x00001fc0 >> 2,
	NV4097_SET_ATTRIB_COLOR                  = 0x00001fc4 >> 2,
	NV4097_SET_ATTRIB_TEX_COORD              = 0x00001fc8 >> 2,
	NV4097_SET_ATTRIB_TEX_COORD_EX           = 0x00001fcc >> 2,
	NV4097_SET_ATTRIB_UCLIP0                 = 0x00001fd0 >> 2,
	NV4097_SET_ATTRIB_UCLIP1                 = 0x00001fd4 >> 2,
	NV4097_INVALIDATE_L2                     = 0x00001fd8 >> 2,
	NV4097_SET_REDUCE_DST_COLOR              = 0x00001fe0 >> 2,
	NV4097_SET_NO_PARANOID_TEXTURE_FETCHES   = 0x00001fe8 >> 2,
	NV4097_SET_SHADER_PACKER                 = 0x00001fec >> 2,
	NV4097_SET_VERTEX_ATTRIB_INPUT_MASK      = 0x00001ff0 >> 2,
	NV4097_SET_VERTEX_ATTRIB_OUTPUT_MASK     = 0x00001ff4 >> 2,
	NV4097_SET_TRANSFORM_BRANCH_BITS         = 0x00001ff8 >> 2,

	// NV03_MEMORY_TO_MEMORY_FORMAT	(NV0039)
	NV0039_SET_OBJECT                        = 0x00002000 >> 2,
	NV0039_SET_CONTEXT_DMA_NOTIFIES          = 0x00002180 >> 2,
	NV0039_SET_CONTEXT_DMA_BUFFER_IN         = 0x00002184 >> 2,
	NV0039_SET_CONTEXT_DMA_BUFFER_OUT        = 0x00002188 >> 2,
	NV0039_OFFSET_IN                         = 0x0000230C >> 2,
	NV0039_OFFSET_OUT                        = 0x00002310 >> 2,
	NV0039_PITCH_IN                          = 0x00002314 >> 2,
	NV0039_PITCH_OUT                         = 0x00002318 >> 2,
	NV0039_LINE_LENGTH_IN                    = 0x0000231C >> 2,
	NV0039_LINE_COUNT                        = 0x00002320 >> 2,
	NV0039_FORMAT                            = 0x00002324 >> 2,
	NV0039_BUFFER_NOTIFY                     = 0x00002328 >> 2,

	// NV30_CONTEXT_SURFACES_2D	(NV3062)
	NV3062_SET_OBJECT                        = 0x00006000 >> 2,
	NV3062_SET_CONTEXT_DMA_NOTIFIES          = 0x00006180 >> 2,
	NV3062_SET_CONTEXT_DMA_IMAGE_SOURCE      = 0x00006184 >> 2,
	NV3062_SET_CONTEXT_DMA_IMAGE_DESTIN      = 0x00006188 >> 2,
	NV3062_SET_COLOR_FORMAT                  = 0x00006300 >> 2,
	NV3062_SET_PITCH                         = 0x00006304 >> 2,
	NV3062_SET_OFFSET_SOURCE                 = 0x00006308 >> 2,
	NV3062_SET_OFFSET_DESTIN                 = 0x0000630C >> 2,

	// NV30_CONTEXT_SURFACE_SWIZZLED (NV309E)
	NV309E_SET_OBJECT                        = 0x00008000 >> 2,
	NV309E_SET_CONTEXT_DMA_NOTIFIES          = 0x00008180 >> 2,
	NV309E_SET_CONTEXT_DMA_IMAGE             = 0x00008184 >> 2,
	NV309E_SET_FORMAT                        = 0x00008300 >> 2,
	NV309E_SET_OFFSET                        = 0x00008304 >> 2,

	// NV30_IMAGE_FROM_CPU (NV308A)
	NV308A_SET_OBJECT                        = 0x0000A000 >> 2,
	NV308A_SET_CONTEXT_DMA_NOTIFIES          = 0x0000A180 >> 2,
	NV308A_SET_CONTEXT_COLOR_KEY             = 0x0000A184 >> 2,
	NV308A_SET_CONTEXT_CLIP_RECTANGLE        = 0x0000A188 >> 2,
	NV308A_SET_CONTEXT_PATTERN               = 0x0000A18C >> 2,
	NV308A_SET_CONTEXT_ROP                   = 0x0000A190 >> 2,
	NV308A_SET_CONTEXT_BETA1                 = 0x0000A194 >> 2,
	NV308A_SET_CONTEXT_BETA4                 = 0x0000A198 >> 2,
	NV308A_SET_CONTEXT_SURFACE               = 0x0000A19C >> 2,
	NV308A_SET_COLOR_CONVERSION              = 0x0000A2F8 >> 2,
	NV308A_SET_OPERATION                     = 0x0000A2FC >> 2,
	NV308A_SET_COLOR_FORMAT                  = 0x0000A300 >> 2,
	NV308A_POINT                             = 0x0000A304 >> 2,
	NV308A_SIZE_OUT                          = 0x0000A308 >> 2,
	NV308A_SIZE_IN                           = 0x0000A30C >> 2,
	NV308A_COLOR                             = 0x0000A400 >> 2,

	// NV30_SCALED_IMAGE_FROM_MEMORY (NV3089)
	NV3089_SET_OBJECT                        = 0x0000C000 >> 2,
	NV3089_SET_CONTEXT_DMA_NOTIFIES          = 0x0000C180 >> 2,
	NV3089_SET_CONTEXT_DMA_IMAGE             = 0x0000C184 >> 2,
	NV3089_SET_CONTEXT_PATTERN               = 0x0000C188 >> 2,
	NV3089_SET_CONTEXT_ROP                   = 0x0000C18C >> 2,
	NV3089_SET_CONTEXT_BETA1                 = 0x0000C190 >> 2,
	NV3089_SET_CONTEXT_BETA4                 = 0x0000C194 >> 2,
	NV3089_SET_CONTEXT_SURFACE               = 0x0000C198 >> 2,
	NV3089_SET_COLOR_CONVERSION              = 0x0000C2FC >> 2,
	NV3089_SET_COLOR_FORMAT                  = 0x0000C300 >> 2,
	NV3089_SET_OPERATION                     = 0x0000C304 >> 2,
	NV3089_CLIP_POINT                        = 0x0000C308 >> 2,
	NV3089_CLIP_SIZE                         = 0x0000C30C >> 2,
	NV3089_IMAGE_OUT_POINT                   = 0x0000C310 >> 2,
	NV3089_IMAGE_OUT_SIZE                    = 0x0000C314 >> 2,
	NV3089_DS_DX                             = 0x0000C318 >> 2,
	NV3089_DT_DY                             = 0x0000C31C >> 2,
	NV3089_IMAGE_IN_SIZE                     = 0x0000C400 >> 2,
	NV3089_IMAGE_IN_FORMAT                   = 0x0000C404 >> 2,
	NV3089_IMAGE_IN_OFFSET                   = 0x0000C408 >> 2,
	NV3089_IMAGE_IN                          = 0x0000C40C >> 2,

	GCM_SET_USER_COMMAND                     = 0x0000EB00 >> 2,

	GCM_FLIP_COMMAND                         = 0x0000FEAC >> 2
};


enum Method
{
	CELL_GCM_METHOD_FLAG_NON_INCREMENT = 0x40000000,
	CELL_GCM_METHOD_FLAG_JUMP = 0x20000000,
	CELL_GCM_METHOD_FLAG_CALL = 0x00000002,
	CELL_GCM_METHOD_FLAG_RETURN = 0x00020000,
};

namespace rsx
{
	template<typename AT>
	static inline u32 make_command(vm::_ptr_base<be_t<u32>, AT>& dst, u32 start_register, std::initializer_list<any32> values)
	{
		*dst++ = start_register << 2 | static_cast<u32>(values.size()) << 18;

		for (const any32& cmd : values)
		{
			*dst++ = cmd.as<u32>();
		}

		return SIZE_32(u32) * (static_cast<u32>(values.size()) + 1);
	}

	template<typename AT>
	static inline u32 make_jump(vm::_ptr_base<be_t<u32>, AT>& dst, u32 offset)
	{
		*dst++ = CELL_GCM_METHOD_FLAG_JUMP | offset;

		return SIZE_32(u32);
	}

	std::string get_method_name(const u32 id);

	std::function<std::string(u32)> get_pretty_printing_function(const u32 id);
}