function F=GFPDF_Twomey_matrix_double_integrand(normZp_full_range,DMAomega,FIMSomega,Zpbinc,ZpDMA,chg,g,D1_nm,T_FIMS,T_DMA,P_FIMS,P_DMA,i)
%**************************************************************************
%       Load transfer function data and counts in different bins
%**************************************************************************
% global FIMS
%**************************************************************************
%                  Integration of the inner integrand
%**************************************************************************
% T=average(FIMS.Ta);
% pressure=average(FIMS.Pabs);
% chg=1;

DMA_normZp=normZp_full_range;
DMA_omega=DMAomega;
FIMS_omega=FIMSomega(i,:);
FIMS_normZp=normZp_full_range;
% DMA_normZp=FIMS.normZp_full_range;
% DMA_omega=FIMS.GFPDF.DMAomega;
% FIMS_omega=FIMS.GFPDF.omega_norm(i,:);
% FIMS_normZp=FIMS.normZp_full_range;

Zpc=Zpbinc(i);
% D2=D2_nm/1e9;
D2=g.*D1_nm/1e9;
D1=D1_nm/1e9;
Zp2=Dp2Zp(D2,T_FIMS+273.15,P_FIMS,chg);
normZp2=Zp2/Zpc;
Zp1=Dp2Zp(D1,T_DMA+273.15,P_DMA,chg);
normZp1=Zp1/ZpDMA;
% normZp1=Zp1/nanmean(FIMS.ZpDMA);

omega_FIMS=interp1(FIMS_normZp,FIMS_omega,normZp2);
tmpidx=isnan(omega_FIMS);
omega_FIMS(tmpidx)=0;
tmpidx2=omega_FIMS<0;
omega_FIMS(tmpidx2)=0;

omega_DMA=interp1(DMA_normZp,DMA_omega,normZp1);
tmpidx3=isnan(omega_DMA);
omega_DMA(tmpidx3)=0;
tmpidx4=omega_DMA<0;
omega_DMA(tmpidx4)=0;

F=omega_FIMS.*omega_DMA;
