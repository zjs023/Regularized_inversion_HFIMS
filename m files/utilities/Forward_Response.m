function [cPDF_sim,R_sim]=Forward_Response(case_num)
%%% this function uses input of x to create a non linear function set so
%%% that least sequare method can be used
global FIMS paramt_FIMS
%%% set boundaries for x
%%% remember to change this accordingly! f,G,sigma
action.plot_simulated_GFPDF = 0;
switch case_num
    case 4 
        fcond=  [1];
        G=      [1.0];
        sigma=  [1.05];
    case 3
        fcond=  [0.39,0.32,0.29];
        G=      [1.05,1.40,1.70];
        sigma=  [1.10,1.05,1.10];
    case 2
        fcond=  [0.45 0.55];
        G=      [1.10 1.30];
        sigma=  [1.05 1.05];
    case 1
        fcond=  [1];
        G=      [1.40];
        sigma=  [1.15];
end

fcond=fcond/sum(fcond);
D2norm=[0.8:0.01:2.0];
cPDF_sim=GFPDF_log(D2norm,1,fcond,G,sigma)./D2norm;
if action.plot_simulated_GFPDF
    figure; plot(D2norm,cPDF_sim,'LineWidth',3);
    xlabel('g'); ylabel('GF-PDF');
end

T_FIMS=FIMS.GF.Tsep;
T_DMA=FIMS.GF.DMA_T;
P_FIMS=FIMS.GF.Pamb;
P_DMA=FIMS.GF.DMA_P;
chg=1;
num=length(FIMS.GF.R);

R_F=zeros(1,num);
for j1=1:1:num    
    Zpc=FIMS.Zpbinc(j1);
    Dplobnd_FIMS=Zp2Dp(Zpc*1.2,T_FIMS+273.15,P_FIMS,chg)*1e9;
    Dpupbnd_FIMS=Zp2Dp(Zpc*0.8,T_FIMS+273.15,P_FIMS,chg)*1e9;
    R_F(j1)=FIMS.GF.factor*FIMS.GF.a/FIMS.GF.DMA_beta*quad2d(@(D2,D1) GFPDF_log_double_integrand(FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
            FIMS.GF.Zpbinc,FIMS.GF.ZpDMA,FIMS.GF.chg,D2,D1,fcond,G,sigma,T_FIMS,T_DMA,P_FIMS,P_DMA,j1),Dplobnd_FIMS,Dpupbnd_FIMS,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
end
        
R_sim=R_F/max(R_F);