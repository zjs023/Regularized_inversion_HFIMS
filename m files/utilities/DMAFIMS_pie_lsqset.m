function F=DMAFIMS_pie_lsqset(x,fit_to,Zpbinc,ZpDMA,factor,T_FIMS,T_DMA,P_FIMS,P_DMA,range)
%%% this function uses input of x to create a non linear function set so
%%% that least sequare method can be used
%%% x(i) are values of the GF-PDF in the non-zero location


%%% read FIMS file
    global FIMS 

    meas_norm=fit_to/max(fit_to)+0.001;
    R_F=zeros(1,FIMS.GF.num);


    for j1=1:1:FIMS.GF.num 
        Zpc=Zpbinc(j1);
        Dplobnd_FIMS=Zp2Dp(Zpc*1.2,T_FIMS+273.15,P_FIMS,FIMS.GF.chg)*1e9;
        Dpupbnd_FIMS=Zp2Dp(Zpc*0.8,T_FIMS+273.15,P_FIMS,FIMS.GF.chg)*1e9;
        R_F(j1)=factor*FIMS.GF.a/FIMS.GF.DMA_beta*quad2d(@(D2,D1) GFPDF_pie_double_integrand(FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                        Zpbinc,ZpDMA,FIMS.GF.chg,D2,D1,x,range,T_FIMS,T_DMA,P_FIMS,P_DMA,j1),Dplobnd_FIMS,Dpupbnd_FIMS,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
    end
    theo_norm=R_F/max(R_F)+0.001;

%     relative=((theo_norm-meas_norm)./(meas_norm))*1e5;
    relative=(theo_norm-meas_norm)./sqrt(meas_norm);
    % absolute=(theo_norm-meas_norm)*1e5;
    % F=max(relative,absolute);
    F=relative;
end
% function F=DMAFIMS_pie_lsqset(x)
% %%% this function uses input of x to create a non linear function set so
% %%% that least sequare method can be used
% %%% x(i) are values of the GF-PDF in the non-zero location
% 
% 
% %%% read FIMS file
% global FIMS paramt_FIMS
% 
% fit_to=FIMS.GFPDF.R_avg;
% meas_norm=fit_to/max(fit_to)+0.001;
% R_F=zeros(1,FIMS.GFPDF.num);
% 
% 
% for j1=1:1:FIMS.GFPDF.num    
%     Zpc=FIMS.GFPDF.Zpbinc(j1);
%     Dplobnd_FIMS=Zp2Dp(Zpc*1.2,FIMS.GFPDF.T+273.15,FIMS.GFPDF.pressure,FIMS.GFPDF.chg)*1e9;
%     Dpupbnd_FIMS=Zp2Dp(Zpc*0.8,FIMS.GFPDF.T+273.15,FIMS.GFPDF.pressure,FIMS.GFPDF.chg)*1e9;
%     R_F(j1)=FIMS.GFPDF.factor*FIMS.GFPDF.a/FIMS.GFPDF.DMA_beta*quad2d(@(D2,D1) GFPDF_pie_double_integrand(D2,D1,x,j1),Dplobnd_FIMS,Dpupbnd_FIMS,FIMS.GFPDF.Dplobnd_DMA,FIMS.GFPDF.Dpupbnd_DMA);
% end
% 
% theo_norm=R_F/max(R_F)+0.001;
% 
% relative=((theo_norm-meas_norm)./(meas_norm))*1e5;
% % absolute=(theo_norm-meas_norm)*1e5;
% % F=max(relative,absolute);
% F=relative;