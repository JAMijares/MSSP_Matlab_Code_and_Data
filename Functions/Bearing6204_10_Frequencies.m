alpha=0;
Bd=0.3125;
Nb=8;
Pd=1.3189;


gamma=(Bd/Pd)*cos(alpha) ;
FTF=0.5*(1-gamma);
BPFI=0.5*Nb*(1+gamma);
BPFO=0.5*Nb*(1-gamma);
BSF=0.5*(Pd/Bd)*(1-gamma^2);
BDF=2*BSF;

% BPFO=3.08;
% BPFI=4.92;
% FTF=0.4%0.39;
% BSF=2.06;