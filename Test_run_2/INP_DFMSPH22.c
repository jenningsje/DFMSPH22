== this line indicates the length of the title ===========================================
*2020.08.29 ==  12C+120Te ================================================================
63.0   1.2     80
 ECM   r00     iC2b
  6   12    52   120
  ZP   AP   ZT   AT       
   30      1         1        3        0        0.
iter_up   key_exc   key_vNN   key_DD   key_UC   AdelCorr
 5.0    4.0    0.00001    15.1     3.1       0.01
k_up   Crup   eps_iter   RCCstart  RCCfin  RCCstep //RCCstart>RCCfin
0.7      1.2    0.01   -2000.  -1000     20    0.4    1.0    0.01   1.
rWSmin rWSmax rWSstep  VWSmin VWSmax VWSstep aWSmin aWSmax aWSstep chi2WS
1.2       0.85      4
fRBstart fRBfin  deliR //fRBstart>fRBfin
1.00     1.5   0.02   0.4      0.5    0.02  
rGKmin rGKmax rGKstep aGKmin aGKmax aGKstep 
2.    40.   0.5      4    11    0.5      5     15     0.5       1.
A0min A0max A0step  A1min A1max A1step A2min A2max A2step chi2GK

******************************************************

key_vNN - NN forces  
0- M3Y Reid; 
1- M3Y Paris; 
2-MIG; 
11 - NL1 forces P.-G. Reinhard et al. ZPA 323 (1986) 13-25; Hirata PRC 44 (1991) 1467 
12 - NL2 forces P.-G. Reinhard et al. ZPA 323 (1986) 13-25; Hirata PRC 44 (1991) 1467 
13 - NL3 forces G. A. Lalazissis, J. K¨onig and P. Ring, Phys. Rev. C 55 (1997) 540
14 - FSUGold forces B. G. Todd-Rutel and J. Piekarewicz PRL 95 (2005) 122501
21 - TM1 forces Y. Sugahara and H. Toki, Nucl. Phys. A 579 (1994) 557
23 - HS forces C.J. Horowitz, B.D. Serot, Nucl. Phys. A 368 (1981) 503–528

********************************************************************
iter_up - upper limit of the number of iterations
key_exc(0|1) - exchange forces  0 - delta excfor; 1 - finite range excfor
key_vNN(0|1|2|3) - NN forces  0-M3Y Reid; 1- M3Y Paris; 2-MIG; 3 - TM1;
key_DD(0|1|2|3|4|5|6|7|8) - Density Dependence 0-no;
		1-yes,set 1;
		2-yes,set 2;...
		F=CDD*[1+alDD*exp(-beDD*rho)-gDD*rho]
key_UC(0|1) - exchange part of Coulom interection: 0 - no; 1 - attends
k_up -  the upper limit for k-integration
Crup -  the upper limit for r-integration is R0*Crup
chi2WS - if >0. appr the DFP by WSP
chi2GK - if >0. appr the DFP by GK-profile