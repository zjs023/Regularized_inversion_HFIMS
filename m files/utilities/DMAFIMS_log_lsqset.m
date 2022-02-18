function F=DMAFIMS_log_lsqset(x,fit_to,Zpbinc,ZpDMA,factor,T_FIMS,T_DMA,P_FIMS,P_DMA)
%%% this function uses input of x to create a non linear function set so
%%% that least sequare method can be used
%%% x(1),x(2),x(3),x(4),x(5),x(6), ... are f1,g1,sigma1,f2,g2,sigma2

    n=length(x)/3;
    fcond=zeros(1,n);
    G=zeros(1,n);
    sigma=zeros(1,n);
    for i=1:1:n
        fcond(i)=x(1+(i-1)*3);
        G(i)=x(2+(i-1)*3);
        sigma(i)=x(3+(i-1)*3);
    end

    %%% read FIMS file
    global FIMS

    % fit_to=FIMS.GF.R_avg;
    meas_norm=fit_to/max(fit_to)+0.001;
    R_F=zeros(1,FIMS.GF.num);

    for j1=1:1:FIMS.GF.num    
        Zpc=Zpbinc(j1);
        Dplobnd_FIMS=Zp2Dp(Zpc*1.2,T_FIMS+273.15,P_FIMS,FIMS.GF.chg)*1e9;
        Dpupbnd_FIMS=Zp2Dp(Zpc*0.8,T_FIMS+273.15,P_FIMS,FIMS.GF.chg)*1e9;
        R_F(j1)=factor*FIMS.GF.a/FIMS.GF.DMA_beta*quad2d(@(D2,D1) GFPDF_log_double_integrand(FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                        Zpbinc,ZpDMA,FIMS.GF.chg,D2,D1,fcond,G,sigma,T_FIMS,T_DMA,P_FIMS,P_DMA,j1),Dplobnd_FIMS,Dpupbnd_FIMS,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
    %     Zpc=FIMS.GF.Zpbinc(j1);
    %     Dplobnd_FIMS=Zp2Dp(Zpc*1.2,FIMS.GF.T+273.15,FIMS.GF.pressure,FIMS.GF.chg)*1e9;
    %     Dpupbnd_FIMS=Zp2Dp(Zpc*0.8,FIMS.GF.T+273.15,FIMS.GF.pressure,FIMS.GF.chg)*1e9;
    %     R_F(j1)=FIMS.GF.factor*FIMS.GF.a/FIMS.GF.DMA_beta*quad2d(@(D2,D1) GFPDF_log_double_integrand(D2,D1,fcond,G,sigma,j1),Dplobnd_FIMS,Dpupbnd_FIMS,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
    end

    theo_norm=R_F/max(R_F)+0.001;
%     relative=((theo_norm-meas_norm).^2./(meas_norm))*1e5;
    relative=((theo_norm-meas_norm)./sqrt(meas_norm));
    % absolute=(theo_norm-meas_norm)*1e5;
    % F=max(relative,absolute);
    F=relative;
end

% function F=DMAFIMS_log_lsqset(x)
% %%% this function uses input of x to create a non linear function set so
% %%% that least sequare method can be used
% %%% x(1),x(2),x(3),x(4),x(5),x(6), ... are f1,g1,sigma1,f2,g2,sigma2
% 
% n=length(x)/3;
% fcond=zeros(1,n);
% G=zeros(1,n);
% sigma=zeros(1,n);
% for i=1:1:n
%     fcond(i)=x(1+(i-1)*3);
%     G(i)=x(2+(i-1)*3);
%     sigma(i)=x(3+(i-1)*3);
% end
% 
% %%% read FIMS file
% global FIMS paramt_FIMS
% 
% fit_to=FIMS.GF.R_avg;
% meas_norm=fit_to/max(fit_to)+0.001;
% R_F=zeros(1,FIMS.GF.num);
% 
% for j1=1:1:FIMS.GF.num    
%     Zpc=FIMS.GF.Zpbinc(j1);
%     Dplobnd_FIMS=Zp2Dp(Zpc*1.2,FIMS.GF.T+273.15,FIMS.GF.pressure,FIMS.GF.chg)*1e9;
%     Dpupbnd_FIMS=Zp2Dp(Zpc*0.8,FIMS.GF.T+273.15,FIMS.GF.pressure,FIMS.GF.chg)*1e9;
%     R_F(j1)=FIMS.GF.factor*FIMS.GF.a/FIMS.GF.DMA_beta*quad2d(@(D2,D1) GFPDF_log_double_integrand(D2,D1,fcond,G,sigma,j1),Dplobnd_FIMS,Dpupbnd_FIMS,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
% end
% 
% theo_norm=R_F/max(R_F)+0.001;
% relative=((theo_norm-meas_norm).^2./(meas_norm))*1e5;
% % absolute=(theo_norm-meas_norm)*1e5;
% % F=max(relative,absolute);
% F=relative;
