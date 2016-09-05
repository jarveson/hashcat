/**
 * Author......: Jens Steube <jens.steube@gmail.com>
 * License.....: MIT
 */

#define COMBINATOR_MODE_BASE_LEFT  10001
#define COMBINATOR_MODE_BASE_RIGHT 10002

#ifdef SHARED_H
#define _AXCRYPT_
#define _BCRYPT_
#define _CLOUDKEY_
#define _KECCAK_
#define _MD4_
#define _MD5_
#define _MD5H_
#define _OFFICE2013_
#define _PBKDF2_MD5_
#define _PBKDF2_SHA1_
#define _PBKDF2_SHA256_
#define _PBKDF2_SHA512_
#define _PDF17L8_
#define _RAR3_
#define _RIPEMD160_
#define _SAPB_
#define _SHA1_
#define _SHA256_
#define _SHA384_
#define _SHA512_
#define _SIPHASH_
#define _WHIRLPOOL_
#define _ZIP2_
#endif

#ifdef _SIPHASH_

typedef enum siphash_constants
{
  SIPHASHM_0=0x736f6d6570736575,
  SIPHASHM_1=0x646f72616e646f6d,
  SIPHASHM_2=0x6c7967656e657261,
  SIPHASHM_3=0x7465646279746573

} siphash_constants_t;

#endif

#if defined _BCRYPT_ || defined _PSAFE2_

typedef enum bcrypt_constants
{
  BCRYPTM_0=0x4F727068u,
  BCRYPTM_1=0x65616E42u,
  BCRYPTM_2=0x65686F6Cu,
  BCRYPTM_3=0x64657253u,
  BCRYPTM_4=0x63727944u,
  BCRYPTM_5=0x6F756274u

} bcrypt_constants_t;

#endif

#if defined _SHA1_ || defined _SAPG_ || defined _OFFICE2007_ || defined _OFFICE2010_ || defined _OLDOFFICE34_ || defined _ANDROIDFDE_ || defined _DCC2_ || defined _WPA_ || defined _MD5_SHA1_ || defined _SHA1_MD5_ || defined _PSAFE2_ || defined _LOTUS8_ || defined _PBKDF2_SHA1_ || defined _RAR3_ || defined _SHA256_SHA1_ || defined _ZIP2_ || defined _AXCRYPT_

typedef enum sha1_constants
{
  SHA1M_A=0x67452301u,
  SHA1M_B=0xefcdab89u,
  SHA1M_C=0x98badcfeu,
  SHA1M_D=0x10325476u,
  SHA1M_E=0xc3d2e1f0u,

  SHA1C00=0x5a827999u,
  SHA1C01=0x6ed9eba1u,
  SHA1C02=0x8f1bbcdcu,
  SHA1C03=0xca62c1d6u

} sha1_constants_t;

#endif

#if defined _SHA256_ || defined _PDF17L8_ || defined _SEVEN_ZIP_ || defined _ANDROIDFDE_ || defined _CLOUDKEY_ || defined _SCRYPT_ || defined _PBKDF2_SHA256_ || defined _SHA256_SHA1_ || defined _MS_DRSR_ || defined _ANDROIDFDE_SAMSUNG_ || defined _RAR5_ || defined _KEEPASS_

typedef enum sha256_constants
{
  SHA256M_A=0x6a09e667u,
  SHA256M_B=0xbb67ae85u,
  SHA256M_C=0x3c6ef372u,
  SHA256M_D=0xa54ff53au,
  SHA256M_E=0x510e527fu,
  SHA256M_F=0x9b05688cu,
  SHA256M_G=0x1f83d9abu,
  SHA256M_H=0x5be0cd19u,

  SHA256C00=0x428a2f98u,
  SHA256C01=0x71374491u,
  SHA256C02=0xb5c0fbcfu,
  SHA256C03=0xe9b5dba5u,
  SHA256C04=0x3956c25bu,
  SHA256C05=0x59f111f1u,
  SHA256C06=0x923f82a4u,
  SHA256C07=0xab1c5ed5u,
  SHA256C08=0xd807aa98u,
  SHA256C09=0x12835b01u,
  SHA256C0a=0x243185beu,
  SHA256C0b=0x550c7dc3u,
  SHA256C0c=0x72be5d74u,
  SHA256C0d=0x80deb1feu,
  SHA256C0e=0x9bdc06a7u,
  SHA256C0f=0xc19bf174u,
  SHA256C10=0xe49b69c1u,
  SHA256C11=0xefbe4786u,
  SHA256C12=0x0fc19dc6u,
  SHA256C13=0x240ca1ccu,
  SHA256C14=0x2de92c6fu,
  SHA256C15=0x4a7484aau,
  SHA256C16=0x5cb0a9dcu,
  SHA256C17=0x76f988dau,
  SHA256C18=0x983e5152u,
  SHA256C19=0xa831c66du,
  SHA256C1a=0xb00327c8u,
  SHA256C1b=0xbf597fc7u,
  SHA256C1c=0xc6e00bf3u,
  SHA256C1d=0xd5a79147u,
  SHA256C1e=0x06ca6351u,
  SHA256C1f=0x14292967u,
  SHA256C20=0x27b70a85u,
  SHA256C21=0x2e1b2138u,
  SHA256C22=0x4d2c6dfcu,
  SHA256C23=0x53380d13u,
  SHA256C24=0x650a7354u,
  SHA256C25=0x766a0abbu,
  SHA256C26=0x81c2c92eu,
  SHA256C27=0x92722c85u,
  SHA256C28=0xa2bfe8a1u,
  SHA256C29=0xa81a664bu,
  SHA256C2a=0xc24b8b70u,
  SHA256C2b=0xc76c51a3u,
  SHA256C2c=0xd192e819u,
  SHA256C2d=0xd6990624u,
  SHA256C2e=0xf40e3585u,
  SHA256C2f=0x106aa070u,
  SHA256C30=0x19a4c116u,
  SHA256C31=0x1e376c08u,
  SHA256C32=0x2748774cu,
  SHA256C33=0x34b0bcb5u,
  SHA256C34=0x391c0cb3u,
  SHA256C35=0x4ed8aa4au,
  SHA256C36=0x5b9cca4fu,
  SHA256C37=0x682e6ff3u,
  SHA256C38=0x748f82eeu,
  SHA256C39=0x78a5636fu,
  SHA256C3a=0x84c87814u,
  SHA256C3b=0x8cc70208u,
  SHA256C3c=0x90befffau,
  SHA256C3d=0xa4506cebu,
  SHA256C3e=0xbef9a3f7u,
  SHA256C3f=0xc67178f2u

} sha256_constants_t;

#endif

#if defined _MD4_ || defined _DCC2_ || defined _NETNTLMV2_ || defined _KRB5PA_ || defined _MS_DRSR_ || defined _KRB5TGS_

typedef enum md4_constants
{
  MD4M_A=0x67452301u,
  MD4M_B=0xefcdab89u,
  MD4M_C=0x98badcfeu,
  MD4M_D=0x10325476u,

  MD4S00=3u,
  MD4S01=7u,
  MD4S02=11u,
  MD4S03=19u,
  MD4S10=3u,
  MD4S11=5u,
  MD4S12=9u,
  MD4S13=13u,
  MD4S20=3u,
  MD4S21=9u,
  MD4S22=11u,
  MD4S23=15u,

  MD4C00=0x00000000u,
  MD4C01=0x5a827999u,
  MD4C02=0x6ed9eba1u

} md4_constants_t;

#endif

#if defined _MD5_ || defined _MD5H_ || defined _SAPB_ || defined _OLDOFFICE01_ || defined _WPA_ || defined _MD5_SHA1_ || defined _SHA1_MD5_ || defined _NETNTLMV2_ || defined _KRB5PA_  || defined _PBKDF2_MD5_ || defined _KRB5TGS_

typedef enum md5_constants
{
  MD5M_A=0x67452301u,
  MD5M_B=0xefcdab89u,
  MD5M_C=0x98badcfeu,
  MD5M_D=0x10325476u,

  MD5S00=7u,
  MD5S01=12u,
  MD5S02=17u,
  MD5S03=22u,
  MD5S10=5u,
  MD5S11=9u,
  MD5S12=14u,
  MD5S13=20u,
  MD5S20=4u,
  MD5S21=11u,
  MD5S22=16u,
  MD5S23=23u,
  MD5S30=6u,
  MD5S31=10u,
  MD5S32=15u,
  MD5S33=21u,

  MD5C00=0xd76aa478u,
  MD5C01=0xe8c7b756u,
  MD5C02=0x242070dbu,
  MD5C03=0xc1bdceeeu,
  MD5C04=0xf57c0fafu,
  MD5C05=0x4787c62au,
  MD5C06=0xa8304613u,
  MD5C07=0xfd469501u,
  MD5C08=0x698098d8u,
  MD5C09=0x8b44f7afu,
  MD5C0a=0xffff5bb1u,
  MD5C0b=0x895cd7beu,
  MD5C0c=0x6b901122u,
  MD5C0d=0xfd987193u,
  MD5C0e=0xa679438eu,
  MD5C0f=0x49b40821u,
  MD5C10=0xf61e2562u,
  MD5C11=0xc040b340u,
  MD5C12=0x265e5a51u,
  MD5C13=0xe9b6c7aau,
  MD5C14=0xd62f105du,
  MD5C15=0x02441453u,
  MD5C16=0xd8a1e681u,
  MD5C17=0xe7d3fbc8u,
  MD5C18=0x21e1cde6u,
  MD5C19=0xc33707d6u,
  MD5C1a=0xf4d50d87u,
  MD5C1b=0x455a14edu,
  MD5C1c=0xa9e3e905u,
  MD5C1d=0xfcefa3f8u,
  MD5C1e=0x676f02d9u,
  MD5C1f=0x8d2a4c8au,
  MD5C20=0xfffa3942u,
  MD5C21=0x8771f681u,
  MD5C22=0x6d9d6122u,
  MD5C23=0xfde5380cu,
  MD5C24=0xa4beea44u,
  MD5C25=0x4bdecfa9u,
  MD5C26=0xf6bb4b60u,
  MD5C27=0xbebfbc70u,
  MD5C28=0x289b7ec6u,
  MD5C29=0xeaa127fau,
  MD5C2a=0xd4ef3085u,
  MD5C2b=0x04881d05u,
  MD5C2c=0xd9d4d039u,
  MD5C2d=0xe6db99e5u,
  MD5C2e=0x1fa27cf8u,
  MD5C2f=0xc4ac5665u,
  MD5C30=0xf4292244u,
  MD5C31=0x432aff97u,
  MD5C32=0xab9423a7u,
  MD5C33=0xfc93a039u,
  MD5C34=0x655b59c3u,
  MD5C35=0x8f0ccc92u,
  MD5C36=0xffeff47du,
  MD5C37=0x85845dd1u,
  MD5C38=0x6fa87e4fu,
  MD5C39=0xfe2ce6e0u,
  MD5C3a=0xa3014314u,
  MD5C3b=0x4e0811a1u,
  MD5C3c=0xf7537e82u,
  MD5C3d=0xbd3af235u,
  MD5C3e=0x2ad7d2bbu,
  MD5C3f=0xeb86d391u

} md5_constants_t;

#endif

#if defined _SHA384_ || defined _PDF17L8_

typedef enum sha384_constants
{
  SHA384M_A=0xcbbb9d5dc1059ed8,
  SHA384M_B=0x629a292a367cd507,
  SHA384M_C=0x9159015a3070dd17,
  SHA384M_D=0x152fecd8f70e5939,
  SHA384M_E=0x67332667ffc00b31,
  SHA384M_F=0x8eb44a8768581511,
  SHA384M_G=0xdb0c2e0d64f98fa7,
  SHA384M_H=0x47b5481dbefa4fa4,

  SHA384C00=0x428a2f98d728ae22,
  SHA384C01=0x7137449123ef65cd,
  SHA384C02=0xb5c0fbcfec4d3b2f,
  SHA384C03=0xe9b5dba58189dbbc,
  SHA384C04=0x3956c25bf348b538,
  SHA384C05=0x59f111f1b605d019,
  SHA384C06=0x923f82a4af194f9b,
  SHA384C07=0xab1c5ed5da6d8118,
  SHA384C08=0xd807aa98a3030242,
  SHA384C09=0x12835b0145706fbe,
  SHA384C0a=0x243185be4ee4b28c,
  SHA384C0b=0x550c7dc3d5ffb4e2,
  SHA384C0c=0x72be5d74f27b896f,
  SHA384C0d=0x80deb1fe3b1696b1,
  SHA384C0e=0x9bdc06a725c71235,
  SHA384C0f=0xc19bf174cf692694,
  SHA384C10=0xe49b69c19ef14ad2,
  SHA384C11=0xefbe4786384f25e3,
  SHA384C12=0x0fc19dc68b8cd5b5,
  SHA384C13=0x240ca1cc77ac9c65,
  SHA384C14=0x2de92c6f592b0275,
  SHA384C15=0x4a7484aa6ea6e483,
  SHA384C16=0x5cb0a9dcbd41fbd4,
  SHA384C17=0x76f988da831153b5,
  SHA384C18=0x983e5152ee66dfab,
  SHA384C19=0xa831c66d2db43210,
  SHA384C1a=0xb00327c898fb213f,
  SHA384C1b=0xbf597fc7beef0ee4,
  SHA384C1c=0xc6e00bf33da88fc2,
  SHA384C1d=0xd5a79147930aa725,
  SHA384C1e=0x06ca6351e003826f,
  SHA384C1f=0x142929670a0e6e70,
  SHA384C20=0x27b70a8546d22ffc,
  SHA384C21=0x2e1b21385c26c926,
  SHA384C22=0x4d2c6dfc5ac42aed,
  SHA384C23=0x53380d139d95b3df,
  SHA384C24=0x650a73548baf63de,
  SHA384C25=0x766a0abb3c77b2a8,
  SHA384C26=0x81c2c92e47edaee6,
  SHA384C27=0x92722c851482353b,
  SHA384C28=0xa2bfe8a14cf10364,
  SHA384C29=0xa81a664bbc423001,
  SHA384C2a=0xc24b8b70d0f89791,
  SHA384C2b=0xc76c51a30654be30,
  SHA384C2c=0xd192e819d6ef5218,
  SHA384C2d=0xd69906245565a910,
  SHA384C2e=0xf40e35855771202a,
  SHA384C2f=0x106aa07032bbd1b8,
  SHA384C30=0x19a4c116b8d2d0c8,
  SHA384C31=0x1e376c085141ab53,
  SHA384C32=0x2748774cdf8eeb99,
  SHA384C33=0x34b0bcb5e19b48a8,
  SHA384C34=0x391c0cb3c5c95a63,
  SHA384C35=0x4ed8aa4ae3418acb,
  SHA384C36=0x5b9cca4f7763e373,
  SHA384C37=0x682e6ff3d6b2b8a3,
  SHA384C38=0x748f82ee5defb2fc,
  SHA384C39=0x78a5636f43172f60,
  SHA384C3a=0x84c87814a1f0ab72,
  SHA384C3b=0x8cc702081a6439ec,
  SHA384C3c=0x90befffa23631e28,
  SHA384C3d=0xa4506cebde82bde9,
  SHA384C3e=0xbef9a3f7b2c67915,
  SHA384C3f=0xc67178f2e372532b,
  SHA384C40=0xca273eceea26619c,
  SHA384C41=0xd186b8c721c0c207,
  SHA384C42=0xeada7dd6cde0eb1e,
  SHA384C43=0xf57d4f7fee6ed178,
  SHA384C44=0x06f067aa72176fba,
  SHA384C45=0x0a637dc5a2c898a6,
  SHA384C46=0x113f9804bef90dae,
  SHA384C47=0x1b710b35131c471b,
  SHA384C48=0x28db77f523047d84,
  SHA384C49=0x32caab7b40c72493,
  SHA384C4a=0x3c9ebe0a15c9bebc,
  SHA384C4b=0x431d67c49c100d4c,
  SHA384C4c=0x4cc5d4becb3e42b6,
  SHA384C4d=0x597f299cfc657e2a,
  SHA384C4e=0x5fcb6fab3ad6faec,
  SHA384C4f=0x6c44198c4a475817

} sha384_constants_t;

#endif

#if defined _SHA512_ || defined _CLOUDKEY_ || defined _OFFICE2013_ || defined _PDF17L8_ || defined _PBKDF2_SHA512_

typedef enum sha512_constants
{
  SHA512M_A=0x6a09e667f3bcc908,
  SHA512M_B=0xbb67ae8584caa73b,
  SHA512M_C=0x3c6ef372fe94f82b,
  SHA512M_D=0xa54ff53a5f1d36f1,
  SHA512M_E=0x510e527fade682d1,
  SHA512M_F=0x9b05688c2b3e6c1f,
  SHA512M_G=0x1f83d9abfb41bd6b,
  SHA512M_H=0x5be0cd19137e2179,

  SHA512C00=0x428a2f98d728ae22,
  SHA512C01=0x7137449123ef65cd,
  SHA512C02=0xb5c0fbcfec4d3b2f,
  SHA512C03=0xe9b5dba58189dbbc,
  SHA512C04=0x3956c25bf348b538,
  SHA512C05=0x59f111f1b605d019,
  SHA512C06=0x923f82a4af194f9b,
  SHA512C07=0xab1c5ed5da6d8118,
  SHA512C08=0xd807aa98a3030242,
  SHA512C09=0x12835b0145706fbe,
  SHA512C0a=0x243185be4ee4b28c,
  SHA512C0b=0x550c7dc3d5ffb4e2,
  SHA512C0c=0x72be5d74f27b896f,
  SHA512C0d=0x80deb1fe3b1696b1,
  SHA512C0e=0x9bdc06a725c71235,
  SHA512C0f=0xc19bf174cf692694,
  SHA512C10=0xe49b69c19ef14ad2,
  SHA512C11=0xefbe4786384f25e3,
  SHA512C12=0x0fc19dc68b8cd5b5,
  SHA512C13=0x240ca1cc77ac9c65,
  SHA512C14=0x2de92c6f592b0275,
  SHA512C15=0x4a7484aa6ea6e483,
  SHA512C16=0x5cb0a9dcbd41fbd4,
  SHA512C17=0x76f988da831153b5,
  SHA512C18=0x983e5152ee66dfab,
  SHA512C19=0xa831c66d2db43210,
  SHA512C1a=0xb00327c898fb213f,
  SHA512C1b=0xbf597fc7beef0ee4,
  SHA512C1c=0xc6e00bf33da88fc2,
  SHA512C1d=0xd5a79147930aa725,
  SHA512C1e=0x06ca6351e003826f,
  SHA512C1f=0x142929670a0e6e70,
  SHA512C20=0x27b70a8546d22ffc,
  SHA512C21=0x2e1b21385c26c926,
  SHA512C22=0x4d2c6dfc5ac42aed,
  SHA512C23=0x53380d139d95b3df,
  SHA512C24=0x650a73548baf63de,
  SHA512C25=0x766a0abb3c77b2a8,
  SHA512C26=0x81c2c92e47edaee6,
  SHA512C27=0x92722c851482353b,
  SHA512C28=0xa2bfe8a14cf10364,
  SHA512C29=0xa81a664bbc423001,
  SHA512C2a=0xc24b8b70d0f89791,
  SHA512C2b=0xc76c51a30654be30,
  SHA512C2c=0xd192e819d6ef5218,
  SHA512C2d=0xd69906245565a910,
  SHA512C2e=0xf40e35855771202a,
  SHA512C2f=0x106aa07032bbd1b8,
  SHA512C30=0x19a4c116b8d2d0c8,
  SHA512C31=0x1e376c085141ab53,
  SHA512C32=0x2748774cdf8eeb99,
  SHA512C33=0x34b0bcb5e19b48a8,
  SHA512C34=0x391c0cb3c5c95a63,
  SHA512C35=0x4ed8aa4ae3418acb,
  SHA512C36=0x5b9cca4f7763e373,
  SHA512C37=0x682e6ff3d6b2b8a3,
  SHA512C38=0x748f82ee5defb2fc,
  SHA512C39=0x78a5636f43172f60,
  SHA512C3a=0x84c87814a1f0ab72,
  SHA512C3b=0x8cc702081a6439ec,
  SHA512C3c=0x90befffa23631e28,
  SHA512C3d=0xa4506cebde82bde9,
  SHA512C3e=0xbef9a3f7b2c67915,
  SHA512C3f=0xc67178f2e372532b,
  SHA512C40=0xca273eceea26619c,
  SHA512C41=0xd186b8c721c0c207,
  SHA512C42=0xeada7dd6cde0eb1e,
  SHA512C43=0xf57d4f7fee6ed178,
  SHA512C44=0x06f067aa72176fba,
  SHA512C45=0x0a637dc5a2c898a6,
  SHA512C46=0x113f9804bef90dae,
  SHA512C47=0x1b710b35131c471b,
  SHA512C48=0x28db77f523047d84,
  SHA512C49=0x32caab7b40c72493,
  SHA512C4a=0x3c9ebe0a15c9bebc,
  SHA512C4b=0x431d67c49c100d4c,
  SHA512C4c=0x4cc5d4becb3e42b6,
  SHA512C4d=0x597f299cfc657e2a,
  SHA512C4e=0x5fcb6fab3ad6faec,
  SHA512C4f=0x6c44198c4a475817

} sha512_constants_t;

#endif

#ifdef _RIPEMD160_

typedef enum ripemd160_constants
{
  RIPEMD160M_A=0x67452301u,
  RIPEMD160M_B=0xefcdab89u,
  RIPEMD160M_C=0x98badcfeu,
  RIPEMD160M_D=0x10325476u,
  RIPEMD160M_E=0xc3d2e1f0u,

  RIPEMD160C00=0x00000000u,
  RIPEMD160C10=0x5a827999u,
  RIPEMD160C20=0x6ed9eba1u,
  RIPEMD160C30=0x8f1bbcdcu,
  RIPEMD160C40=0xa953fd4eu,
  RIPEMD160C50=0x50a28be6u,
  RIPEMD160C60=0x5c4dd124u,
  RIPEMD160C70=0x6d703ef3u,
  RIPEMD160C80=0x7a6d76e9u,
  RIPEMD160C90=0x00000000u,

  RIPEMD160S00=11u,
  RIPEMD160S01=14u,
  RIPEMD160S02=15u,
  RIPEMD160S03=12u,
  RIPEMD160S04=5u,
  RIPEMD160S05=8u,
  RIPEMD160S06=7u,
  RIPEMD160S07=9u,
  RIPEMD160S08=11u,
  RIPEMD160S09=13u,
  RIPEMD160S0A=14u,
  RIPEMD160S0B=15u,
  RIPEMD160S0C=6u,
  RIPEMD160S0D=7u,
  RIPEMD160S0E=9u,
  RIPEMD160S0F=8u,

  RIPEMD160S10=7u,
  RIPEMD160S11=6u,
  RIPEMD160S12=8u,
  RIPEMD160S13=13u,
  RIPEMD160S14=11u,
  RIPEMD160S15=9u,
  RIPEMD160S16=7u,
  RIPEMD160S17=15u,
  RIPEMD160S18=7u,
  RIPEMD160S19=12u,
  RIPEMD160S1A=15u,
  RIPEMD160S1B=9u,
  RIPEMD160S1C=11u,
  RIPEMD160S1D=7u,
  RIPEMD160S1E=13u,
  RIPEMD160S1F=12u,

  RIPEMD160S20=11u,
  RIPEMD160S21=13u,
  RIPEMD160S22=6u,
  RIPEMD160S23=7u,
  RIPEMD160S24=14u,
  RIPEMD160S25=9u,
  RIPEMD160S26=13u,
  RIPEMD160S27=15u,
  RIPEMD160S28=14u,
  RIPEMD160S29=8u,
  RIPEMD160S2A=13u,
  RIPEMD160S2B=6u,
  RIPEMD160S2C=5u,
  RIPEMD160S2D=12u,
  RIPEMD160S2E=7u,
  RIPEMD160S2F=5u,

  RIPEMD160S30=11u,
  RIPEMD160S31=12u,
  RIPEMD160S32=14u,
  RIPEMD160S33=15u,
  RIPEMD160S34=14u,
  RIPEMD160S35=15u,
  RIPEMD160S36=9u,
  RIPEMD160S37=8u,
  RIPEMD160S38=9u,
  RIPEMD160S39=14u,
  RIPEMD160S3A=5u,
  RIPEMD160S3B=6u,
  RIPEMD160S3C=8u,
  RIPEMD160S3D=6u,
  RIPEMD160S3E=5u,
  RIPEMD160S3F=12u,

  RIPEMD160S40=9u,
  RIPEMD160S41=15u,
  RIPEMD160S42=5u,
  RIPEMD160S43=11u,
  RIPEMD160S44=6u,
  RIPEMD160S45=8u,
  RIPEMD160S46=13u,
  RIPEMD160S47=12u,
  RIPEMD160S48=5u,
  RIPEMD160S49=12u,
  RIPEMD160S4A=13u,
  RIPEMD160S4B=14u,
  RIPEMD160S4C=11u,
  RIPEMD160S4D=8u,
  RIPEMD160S4E=5u,
  RIPEMD160S4F=6u,

  RIPEMD160S50=8u,
  RIPEMD160S51=9u,
  RIPEMD160S52=9u,
  RIPEMD160S53=11u,
  RIPEMD160S54=13u,
  RIPEMD160S55=15u,
  RIPEMD160S56=15u,
  RIPEMD160S57=5u,
  RIPEMD160S58=7u,
  RIPEMD160S59=7u,
  RIPEMD160S5A=8u,
  RIPEMD160S5B=11u,
  RIPEMD160S5C=14u,
  RIPEMD160S5D=14u,
  RIPEMD160S5E=12u,
  RIPEMD160S5F=6u,

  RIPEMD160S60=9u,
  RIPEMD160S61=13u,
  RIPEMD160S62=15u,
  RIPEMD160S63=7u,
  RIPEMD160S64=12u,
  RIPEMD160S65=8u,
  RIPEMD160S66=9u,
  RIPEMD160S67=11u,
  RIPEMD160S68=7u,
  RIPEMD160S69=7u,
  RIPEMD160S6A=12u,
  RIPEMD160S6B=7u,
  RIPEMD160S6C=6u,
  RIPEMD160S6D=15u,
  RIPEMD160S6E=13u,
  RIPEMD160S6F=11u,

  RIPEMD160S70=9u,
  RIPEMD160S71=7u,
  RIPEMD160S72=15u,
  RIPEMD160S73=11u,
  RIPEMD160S74=8u,
  RIPEMD160S75=6u,
  RIPEMD160S76=6u,
  RIPEMD160S77=14u,
  RIPEMD160S78=12u,
  RIPEMD160S79=13u,
  RIPEMD160S7A=5u,
  RIPEMD160S7B=14u,
  RIPEMD160S7C=13u,
  RIPEMD160S7D=13u,
  RIPEMD160S7E=7u,
  RIPEMD160S7F=5u,

  RIPEMD160S80=15u,
  RIPEMD160S81=5u,
  RIPEMD160S82=8u,
  RIPEMD160S83=11u,
  RIPEMD160S84=14u,
  RIPEMD160S85=14u,
  RIPEMD160S86=6u,
  RIPEMD160S87=14u,
  RIPEMD160S88=6u,
  RIPEMD160S89=9u,
  RIPEMD160S8A=12u,
  RIPEMD160S8B=9u,
  RIPEMD160S8C=12u,
  RIPEMD160S8D=5u,
  RIPEMD160S8E=15u,
  RIPEMD160S8F=8u,

  RIPEMD160S90=8u,
  RIPEMD160S91=5u,
  RIPEMD160S92=12u,
  RIPEMD160S93=9u,
  RIPEMD160S94=12u,
  RIPEMD160S95=5u,
  RIPEMD160S96=14u,
  RIPEMD160S97=6u,
  RIPEMD160S98=8u,
  RIPEMD160S99=13u,
  RIPEMD160S9A=6u,
  RIPEMD160S9B=5u,
  RIPEMD160S9C=15u,
  RIPEMD160S9D=13u,
  RIPEMD160S9E=11u,
  RIPEMD160S9F=11u

} ripemd160_constants_t;

#endif

#ifdef _KECCAK_

typedef enum keccak_constants
{
  KECCAK_RNDC_00=0x0000000000000001,
  KECCAK_RNDC_01=0x0000000000008082,
  KECCAK_RNDC_02=0x000000000000808a,
  KECCAK_RNDC_03=0x0000000080008000,
  KECCAK_RNDC_04=0x000000000000808b,
  KECCAK_RNDC_05=0x0000000080000001,
  KECCAK_RNDC_06=0x0000000080008081,
  KECCAK_RNDC_07=0x0000000000008009,
  KECCAK_RNDC_08=0x000000000000008a,
  KECCAK_RNDC_09=0x0000000000000088,
  KECCAK_RNDC_10=0x0000000080008009,
  KECCAK_RNDC_11=0x000000008000000a,
  KECCAK_RNDC_12=0x000000008000808b,
  KECCAK_RNDC_13=0x000000000000008b,
  KECCAK_RNDC_14=0x0000000000008089,
  KECCAK_RNDC_15=0x0000000000008003,
  KECCAK_RNDC_16=0x0000000000008002,
  KECCAK_RNDC_17=0x0000000000000080,
  KECCAK_RNDC_18=0x000000000000800a,
  KECCAK_RNDC_19=0x000000008000000a,
  KECCAK_RNDC_20=0x0000000080008081,
  KECCAK_RNDC_21=0x0000000000008080,
  KECCAK_RNDC_22=0x0000000080000001,
  KECCAK_RNDC_23=0x0000000080008008,

  KECCAK_PILN_00=10u,
  KECCAK_PILN_01=7u,
  KECCAK_PILN_02=11u,
  KECCAK_PILN_03=17u,
  KECCAK_PILN_04=18u,
  KECCAK_PILN_05=3u,
  KECCAK_PILN_06=5u,
  KECCAK_PILN_07=16u,
  KECCAK_PILN_08=8u,
  KECCAK_PILN_09=21u,
  KECCAK_PILN_10=24u,
  KECCAK_PILN_11=4u,
  KECCAK_PILN_12=15u,
  KECCAK_PILN_13=23u,
  KECCAK_PILN_14=19u,
  KECCAK_PILN_15=13u,
  KECCAK_PILN_16=12u,
  KECCAK_PILN_17=2u,
  KECCAK_PILN_18=20u,
  KECCAK_PILN_19=14u,
  KECCAK_PILN_20=22u,
  KECCAK_PILN_21=9u,
  KECCAK_PILN_22=6u,
  KECCAK_PILN_23=1u,

  KECCAK_ROTC_00=1u,
  KECCAK_ROTC_01=3u,
  KECCAK_ROTC_02=6u,
  KECCAK_ROTC_03=10u,
  KECCAK_ROTC_04=15u,
  KECCAK_ROTC_05=21u,
  KECCAK_ROTC_06=28u,
  KECCAK_ROTC_07=36u,
  KECCAK_ROTC_08=45u,
  KECCAK_ROTC_09=55u,
  KECCAK_ROTC_10=2u,
  KECCAK_ROTC_11=14u,
  KECCAK_ROTC_12=27u,
  KECCAK_ROTC_13=41u,
  KECCAK_ROTC_14=56u,
  KECCAK_ROTC_15=8u,
  KECCAK_ROTC_16=25u,
  KECCAK_ROTC_17=43u,
  KECCAK_ROTC_18=62u,
  KECCAK_ROTC_19=18u,
  KECCAK_ROTC_20=39u,
  KECCAK_ROTC_21=61u,
  KECCAK_ROTC_22=20u,
  KECCAK_ROTC_23=44u,

} keccak_constants_t;

#endif

#ifdef _MYSQL323_

typedef enum mysql323_constants
{
  MYSQL323_A=0x50305735u,
  MYSQL323_B=0x12345671u

} mysql323_constants_t;

#endif