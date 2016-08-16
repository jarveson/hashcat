#pragma once

#include "SPUOpcodes.h"

class SPUThread;

using spu_inter_func_t = void(*)(SPUThread& spu, spu_opcode_t op);

struct spu_interpreter
{
	static void UNK(SPUThread&, spu_opcode_t);
	static void set_interrupt_status(SPUThread&, spu_opcode_t);

	static void STOP(SPUThread&, spu_opcode_t);
	static void LNOP(SPUThread&, spu_opcode_t);
	static void SYNC(SPUThread&, spu_opcode_t);
	static void DSYNC(SPUThread&, spu_opcode_t);
	static void MFSPR(SPUThread&, spu_opcode_t);
	static void RDCH(SPUThread&, spu_opcode_t);
	static void RCHCNT(SPUThread&, spu_opcode_t);
	static void SF(SPUThread&, spu_opcode_t);
	static void OR(SPUThread&, spu_opcode_t);
	static void BG(SPUThread&, spu_opcode_t);
	static void SFH(SPUThread&, spu_opcode_t);
	static void NOR(SPUThread&, spu_opcode_t);
	static void ABSDB(SPUThread&, spu_opcode_t);
	static void ROT(SPUThread&, spu_opcode_t);
	static void ROTM(SPUThread&, spu_opcode_t);
	static void ROTMA(SPUThread&, spu_opcode_t);
	static void SHL(SPUThread&, spu_opcode_t);
	static void ROTH(SPUThread&, spu_opcode_t);
	static void ROTHM(SPUThread&, spu_opcode_t);
	static void ROTMAH(SPUThread&, spu_opcode_t);
	static void SHLH(SPUThread&, spu_opcode_t);
	static void ROTI(SPUThread&, spu_opcode_t);
	static void ROTMI(SPUThread&, spu_opcode_t);
	static void ROTMAI(SPUThread&, spu_opcode_t);
	static void SHLI(SPUThread&, spu_opcode_t);
	static void ROTHI(SPUThread&, spu_opcode_t);
	static void ROTHMI(SPUThread&, spu_opcode_t);
	static void ROTMAHI(SPUThread&, spu_opcode_t);
	static void SHLHI(SPUThread&, spu_opcode_t);
	static void A(SPUThread&, spu_opcode_t);
	static void AND(SPUThread&, spu_opcode_t);
	static void CG(SPUThread&, spu_opcode_t);
	static void AH(SPUThread&, spu_opcode_t);
	static void NAND(SPUThread&, spu_opcode_t);
	static void AVGB(SPUThread&, spu_opcode_t);
	static void MTSPR(SPUThread&, spu_opcode_t);
	static void WRCH(SPUThread&, spu_opcode_t);
	static void BIZ(SPUThread&, spu_opcode_t);
	static void BINZ(SPUThread&, spu_opcode_t);
	static void BIHZ(SPUThread&, spu_opcode_t);
	static void BIHNZ(SPUThread&, spu_opcode_t);
	static void STOPD(SPUThread&, spu_opcode_t);
	static void STQX(SPUThread&, spu_opcode_t);
	static void BI(SPUThread&, spu_opcode_t);
	static void BISL(SPUThread&, spu_opcode_t);
	static void IRET(SPUThread&, spu_opcode_t);
	static void BISLED(SPUThread&, spu_opcode_t);
	static void HBR(SPUThread&, spu_opcode_t);
	static void GB(SPUThread&, spu_opcode_t);
	static void GBH(SPUThread&, spu_opcode_t);
	static void GBB(SPUThread&, spu_opcode_t);
	static void FSM(SPUThread&, spu_opcode_t);
	static void FSMH(SPUThread&, spu_opcode_t);
	static void FSMB(SPUThread&, spu_opcode_t);
	static void LQX(SPUThread&, spu_opcode_t);
	static void ROTQBYBI(SPUThread&, spu_opcode_t);
	static void ROTQMBYBI(SPUThread&, spu_opcode_t);
	static void SHLQBYBI(SPUThread&, spu_opcode_t);
	static void CBX(SPUThread&, spu_opcode_t);
	static void CHX(SPUThread&, spu_opcode_t);
	static void CWX(SPUThread&, spu_opcode_t);
	static void CDX(SPUThread&, spu_opcode_t);
	static void ROTQBI(SPUThread&, spu_opcode_t);
	static void ROTQMBI(SPUThread&, spu_opcode_t);
	static void SHLQBI(SPUThread&, spu_opcode_t);
	static void ROTQBY(SPUThread&, spu_opcode_t);
	static void ROTQMBY(SPUThread&, spu_opcode_t);
	static void SHLQBY(SPUThread&, spu_opcode_t);
	static void ORX(SPUThread&, spu_opcode_t);
	static void CBD(SPUThread&, spu_opcode_t);
	static void CHD(SPUThread&, spu_opcode_t);
	static void CWD(SPUThread&, spu_opcode_t);
	static void CDD(SPUThread&, spu_opcode_t);
	static void ROTQBII(SPUThread&, spu_opcode_t);
	static void ROTQMBII(SPUThread&, spu_opcode_t);
	static void SHLQBII(SPUThread&, spu_opcode_t);
	static void ROTQBYI(SPUThread&, spu_opcode_t);
	static void ROTQMBYI(SPUThread&, spu_opcode_t);
	static void SHLQBYI(SPUThread&, spu_opcode_t);
	static void NOP(SPUThread&, spu_opcode_t);
	static void CGT(SPUThread&, spu_opcode_t);
	static void XOR(SPUThread&, spu_opcode_t);
	static void CGTH(SPUThread&, spu_opcode_t);
	static void EQV(SPUThread&, spu_opcode_t);
	static void CGTB(SPUThread&, spu_opcode_t);
	static void SUMB(SPUThread&, spu_opcode_t);
	static void HGT(SPUThread&, spu_opcode_t);
	static void CLZ(SPUThread&, spu_opcode_t);
	static void XSWD(SPUThread&, spu_opcode_t);
	static void XSHW(SPUThread&, spu_opcode_t);
	static void CNTB(SPUThread&, spu_opcode_t);
	static void XSBH(SPUThread&, spu_opcode_t);
	static void CLGT(SPUThread&, spu_opcode_t);
	static void ANDC(SPUThread&, spu_opcode_t);
	static void CLGTH(SPUThread&, spu_opcode_t);
	static void ORC(SPUThread&, spu_opcode_t);
	static void CLGTB(SPUThread&, spu_opcode_t);
	static void HLGT(SPUThread&, spu_opcode_t);
	static void CEQ(SPUThread&, spu_opcode_t);
	static void MPYHHU(SPUThread&, spu_opcode_t);
	static void ADDX(SPUThread&, spu_opcode_t);
	static void SFX(SPUThread&, spu_opcode_t);
	static void CGX(SPUThread&, spu_opcode_t);
	static void BGX(SPUThread&, spu_opcode_t);
	static void MPYHHA(SPUThread&, spu_opcode_t);
	static void MPYHHAU(SPUThread&, spu_opcode_t);
	static void MPY(SPUThread&, spu_opcode_t);
	static void MPYH(SPUThread&, spu_opcode_t);
	static void MPYHH(SPUThread&, spu_opcode_t);
	static void MPYS(SPUThread&, spu_opcode_t);
	static void CEQH(SPUThread&, spu_opcode_t);
	static void MPYU(SPUThread&, spu_opcode_t);
	static void CEQB(SPUThread&, spu_opcode_t);
	static void HEQ(SPUThread&, spu_opcode_t);
	static void BRZ(SPUThread&, spu_opcode_t);
	static void STQA(SPUThread&, spu_opcode_t);
	static void BRNZ(SPUThread&, spu_opcode_t);
	static void BRHZ(SPUThread&, spu_opcode_t);
	static void BRHNZ(SPUThread&, spu_opcode_t);
	static void STQR(SPUThread&, spu_opcode_t);
	static void BRA(SPUThread&, spu_opcode_t);
	static void LQA(SPUThread&, spu_opcode_t);
	static void BRASL(SPUThread&, spu_opcode_t);
	static void BR(SPUThread&, spu_opcode_t);
	static void FSMBI(SPUThread&, spu_opcode_t);
	static void BRSL(SPUThread&, spu_opcode_t);
	static void LQR(SPUThread&, spu_opcode_t);
	static void IL(SPUThread&, spu_opcode_t);
	static void ILHU(SPUThread&, spu_opcode_t);
	static void ILH(SPUThread&, spu_opcode_t);
	static void IOHL(SPUThread&, spu_opcode_t);
	static void ORI(SPUThread&, spu_opcode_t);
	static void ORHI(SPUThread&, spu_opcode_t);
	static void ORBI(SPUThread&, spu_opcode_t);
	static void SFI(SPUThread&, spu_opcode_t);
	static void SFHI(SPUThread&, spu_opcode_t);
	static void ANDI(SPUThread&, spu_opcode_t);
	static void ANDHI(SPUThread&, spu_opcode_t);
	static void ANDBI(SPUThread&, spu_opcode_t);
	static void AI(SPUThread&, spu_opcode_t);
	static void AHI(SPUThread&, spu_opcode_t);
	static void STQD(SPUThread&, spu_opcode_t);
	static void LQD(SPUThread&, spu_opcode_t);
	static void XORI(SPUThread&, spu_opcode_t);
	static void XORHI(SPUThread&, spu_opcode_t);
	static void XORBI(SPUThread&, spu_opcode_t);
	static void CGTI(SPUThread&, spu_opcode_t);
	static void CGTHI(SPUThread&, spu_opcode_t);
	static void CGTBI(SPUThread&, spu_opcode_t);
	static void HGTI(SPUThread&, spu_opcode_t);
	static void CLGTI(SPUThread&, spu_opcode_t);
	static void CLGTHI(SPUThread&, spu_opcode_t);
	static void CLGTBI(SPUThread&, spu_opcode_t);
	static void HLGTI(SPUThread&, spu_opcode_t);
	static void MPYI(SPUThread&, spu_opcode_t);
	static void MPYUI(SPUThread&, spu_opcode_t);
	static void CEQI(SPUThread&, spu_opcode_t);
	static void CEQHI(SPUThread&, spu_opcode_t);
	static void CEQBI(SPUThread&, spu_opcode_t);
	static void HEQI(SPUThread&, spu_opcode_t);
	static void HBRA(SPUThread&, spu_opcode_t);
	static void HBRR(SPUThread&, spu_opcode_t);
	static void ILA(SPUThread&, spu_opcode_t);
	static void SELB(SPUThread&, spu_opcode_t);
	static void SHUFB(SPUThread&, spu_opcode_t);
	static void MPYA(SPUThread&, spu_opcode_t);
	static void DFCGT(SPUThread&, spu_opcode_t);
	static void DFCMGT(SPUThread&, spu_opcode_t);
	static void DFTSV(SPUThread&, spu_opcode_t);
	static void DFCEQ(SPUThread&, spu_opcode_t);
	static void DFCMEQ(SPUThread&, spu_opcode_t);
};

struct spu_interpreter_fast final : spu_interpreter
{
	static void FREST(SPUThread&, spu_opcode_t);
	static void FRSQEST(SPUThread&, spu_opcode_t);
	static void FCGT(SPUThread&, spu_opcode_t);
	static void FA(SPUThread&, spu_opcode_t);
	static void FS(SPUThread&, spu_opcode_t);
	static void FM(SPUThread&, spu_opcode_t);
	static void FCMGT(SPUThread&, spu_opcode_t);
	static void DFA(SPUThread&, spu_opcode_t);
	static void DFS(SPUThread&, spu_opcode_t);
	static void DFM(SPUThread&, spu_opcode_t);
	static void DFMA(SPUThread&, spu_opcode_t);
	static void DFMS(SPUThread&, spu_opcode_t);
	static void DFNMS(SPUThread&, spu_opcode_t);
	static void DFNMA(SPUThread&, spu_opcode_t);
	static void FSCRRD(SPUThread&, spu_opcode_t);
	static void FESD(SPUThread&, spu_opcode_t);
	static void FRDS(SPUThread&, spu_opcode_t);
	static void FSCRWR(SPUThread&, spu_opcode_t);
	static void FCEQ(SPUThread&, spu_opcode_t);
	static void FCMEQ(SPUThread&, spu_opcode_t);
	static void FI(SPUThread&, spu_opcode_t);
	static void CFLTS(SPUThread&, spu_opcode_t);
	static void CFLTU(SPUThread&, spu_opcode_t);
	static void CSFLT(SPUThread&, spu_opcode_t);
	static void CUFLT(SPUThread&, spu_opcode_t);
	static void FNMS(SPUThread&, spu_opcode_t);
	static void FMA(SPUThread&, spu_opcode_t);
	static void FMS(SPUThread&, spu_opcode_t);
};

struct spu_interpreter_precise final : spu_interpreter
{
	static void FREST(SPUThread&, spu_opcode_t);
	static void FRSQEST(SPUThread&, spu_opcode_t);
	static void FCGT(SPUThread&, spu_opcode_t);
	static void FA(SPUThread&, spu_opcode_t);
	static void FS(SPUThread&, spu_opcode_t);
	static void FM(SPUThread&, spu_opcode_t);
	static void FCMGT(SPUThread&, spu_opcode_t);
	static void DFA(SPUThread&, spu_opcode_t);
	static void DFS(SPUThread&, spu_opcode_t);
	static void DFM(SPUThread&, spu_opcode_t);
	static void DFMA(SPUThread&, spu_opcode_t);
	static void DFMS(SPUThread&, spu_opcode_t);
	static void DFNMS(SPUThread&, spu_opcode_t);
	static void DFNMA(SPUThread&, spu_opcode_t);
	static void FSCRRD(SPUThread&, spu_opcode_t);
	static void FESD(SPUThread&, spu_opcode_t);
	static void FRDS(SPUThread&, spu_opcode_t);
	static void FSCRWR(SPUThread&, spu_opcode_t);
	static void FCEQ(SPUThread&, spu_opcode_t);
	static void FCMEQ(SPUThread&, spu_opcode_t);
	static void FI(SPUThread&, spu_opcode_t);
	static void CFLTS(SPUThread&, spu_opcode_t);
	static void CFLTU(SPUThread&, spu_opcode_t);
	static void CSFLT(SPUThread&, spu_opcode_t);
	static void CUFLT(SPUThread&, spu_opcode_t);
	static void FNMS(SPUThread&, spu_opcode_t);
	static void FMA(SPUThread&, spu_opcode_t);
	static void FMS(SPUThread&, spu_opcode_t);
};