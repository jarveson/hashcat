REM put PARAM.SFO and PARAM.PFD in same directory as hashcat executable
.\hashcat64.exe --potfile-disable --hex-wordlist -w3  -m150 -a0 param.sfo paramkey.dict
pause