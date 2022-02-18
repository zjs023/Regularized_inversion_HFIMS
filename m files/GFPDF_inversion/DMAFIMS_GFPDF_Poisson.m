function [g,cPDF,R_inv,t_run,varargout]=DMAFIMS_GFPDF_Poisson(FIMS,option,varargin)

if nargin>1 %change parameters if specified
    fns=fieldnames(option);
    for j1=1:length(fns), fn=fns{j1};	eval([fn,'=option.',fn,';']); end
end

action.show_fitting=                        0;
action.show_xy_histogram=                   0;

T_FIMS=FIMS.GF.Tsep;
T_DMA=FIMS.GF.DMA_T;
P_FIMS=FIMS.GF.Pamb;
P_DMA=FIMS.GF.DMA_P;
ZpDMA=FIMS.GF.ZpDMA;
factor=FIMS.GF.factor;

%%% lsq inversion
if option.lsq_inversion
    GF_channel_bin = varargin{1};
    GFbinc=linspace(0.8,2.0,GF_channel_bin);
    GFbin_bnds=Dp2Dp_set(GFbinc);
    tic
    for i1=1:GF_channel_bin
        GFlobnd=GFbin_bnds(1,i1);
        GFupbnd=GFbin_bnds(3,i1);
        for j1=1:1:FIMS.GF.num 
            Gamma(j1,i1)=quad2d(@(g,D1) GFPDF_Twomey_matrix_double_integrand ...
                (FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                FIMS.Zpbinc,ZpDMA,FIMS.GF.chg,g,D1,T_FIMS,T_DMA,P_FIMS,P_DMA,j1), ...
                GFlobnd,GFupbnd,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
        end
    end
    R = FIMS.GF.R(:)/max(FIMS.GF.R(:));
    c_real = interp1(FIMS.GF.gNorm, FIMS.GF.c_sim, GFbinc);
    FIMS.GF.cPDF=lsqnonneg(Gamma,R(:));
	FIMS.GF.cPDF=FIMS.GF.cPDF(:)./(sum(FIMS.GF.cPDF(:))*(GFbinc(2)-GFbinc(1)));
    g=GFbinc;
    FIMS.GF.R_inv_lsq=Gamma*FIMS.GF.cPDF;
    R_inv = FIMS.GF.R_inv_lsq(:)/max(FIMS.GF.R_inv_lsq(:));
    cPDF = FIMS.GF.cPDF(:);
    FIMS.GF.residue_lsq=sum(power((FIMS.GF.R(:)./max(FIMS.GF.R(:))-FIMS.GF.R_inv_lsq./max(FIMS.GF.R_inv_lsq)),2));
    t_run = toc;
    if action.show_fitting
        figure; 
        subplot(2,1,1); plot(FIMS.GF.Dpbinc,FIMS.GF.R/max(FIMS.GF.R),'o',FIMS.GF.Dpbinc,FIMS.GF.R_inv_lsq/max(FIMS.GF.R_inv_lsq),'.-');
        title(['LSQ: ','RH = ',num2str(FIMS.GF.RHmix),'%, Dry size = ', num2str(FIMS.GF.DpDMA*1e9), 'nm']);
        xlabel('Dp (nm)');ylabel('Normalized R')
        subplot(2,1,2); plot(g,cPDF);hold on; plot(g,c_real);legend('{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
        xlabel('g'); ylabel('c(g,Dp1)');
    end
end
%%% tikhonov 0th inversion
if option.tikhonov0_inversion
    GF_channel_bin = varargin{1};
    GFbinc=linspace(0.8,2.0,GF_channel_bin);
    GFbin_bnds=Dp2Dp_set(GFbinc);
    tic
    for i1=1:GF_channel_bin
        GFlobnd=GFbin_bnds(1,i1);
        GFupbnd=GFbin_bnds(3,i1);
        for j1=1:1:FIMS.GF.num 
            Gamma(j1,i1)=quad2d(@(g,D1) GFPDF_Twomey_matrix_double_integrand ...
                (FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                FIMS.Zpbinc,ZpDMA,FIMS.GF.chg,g,D1,T_FIMS,T_DMA,P_FIMS,P_DMA,j1), ...
                GFlobnd,GFupbnd,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
        end
    end
    R = FIMS.GF.R(:)/max(FIMS.GF.R(:));
    order = 0;
    c_real = interp1(FIMS.GF.gNorm, FIMS.GF.c_sim, GFbinc);
    L = full(get_l(length(Gamma(1,:)), order));
    [U, sm, XX, V] = cgsvd(Gamma, L); % Generalized singular value decomposition in compact form
    lc0 = l_curve(U,sm,R,'Tikh');
%    lambda = my_logspace(1e-4,1,100);
%     lc0 = hankraus(Gamma, R, lambda, L, GFbinc);
%     lc0 = optimal_lambda(Gamma, R, lambda, L, GFbinc, c_real);
%     fprintf('Regularization parameter from L-Curve = %f\n',lc0);
    FIMS.GF.cPDF = ncsolve(Gamma, R, lc0, L);
    FIMS.GF.cPDF=FIMS.GF.cPDF(:)./(sum(FIMS.GF.cPDF(:))*(GFbinc(2)-GFbinc(1)));
    g=GFbinc;
    FIMS.GF.R_inv_tk=Gamma*FIMS.GF.cPDF;
    R_inv = FIMS.GF.R_inv_tk(:)/max(FIMS.GF.R_inv_tk(:));
    cPDF = FIMS.GF.cPDF(:);
    FIMS.GF.residue_tk=sum(power((FIMS.GF.R(:)/max(FIMS.GF.R(:))-FIMS.GF.R_inv_tk/max(FIMS.GF.R_inv_tk)),2));
    t_run = toc;
    if action.show_fitting
        figure; 
        subplot(2,1,1); plot(FIMS.GF.Dpbinc,FIMS.GF.R/max(FIMS.GF.R),'o',FIMS.GF.Dpbinc,FIMS.GF.R_inv_tk/max(FIMS.GF.R_inv_tk),'.-');
        title(['Tikhonov: ', 'RH = ',num2str(FIMS.GF.RHmix),'%, Dry size = ', num2str(FIMS.GF.DpDMA*1e9), 'nm']);
        xlabel('Dp (nm)');ylabel('Normalized R')
        subplot(2,1,2); plot(g,cPDF);hold on; plot(g,c_real);legend('{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
        xlabel('g'); ylabel('c(g,Dp1)');
        hold on;
    end
end
%%% tikhonov 1th inversion
if option.tikhonov1_inversion
    GF_channel_bin = varargin{1};
    GFbinc=linspace(0.8,2.0,GF_channel_bin);
    GFbin_bnds=Dp2Dp_set(GFbinc);
    tic
    for i1=1:GF_channel_bin
        GFlobnd=GFbin_bnds(1,i1);
        GFupbnd=GFbin_bnds(3,i1);
        for j1=1:1:FIMS.GF.num 
            Gamma(j1,i1)=quad2d(@(g,D1) GFPDF_Twomey_matrix_double_integrand ...
                (FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                FIMS.Zpbinc,ZpDMA,FIMS.GF.chg,g,D1,T_FIMS,T_DMA,P_FIMS,P_DMA,j1), ...
                GFlobnd,GFupbnd,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
        end
    end
    R = FIMS.GF.R(:)/max(FIMS.GF.R(:));
    c_real = interp1(FIMS.GF.gNorm, FIMS.GF.c_sim, GFbinc);
    order = 1;
    L = full(get_l(length(Gamma(1,:)), order));
   [U, sm, XX, V] = cgsvd(Gamma, L); % Generalized singular value decomposition in compact form
   lc0 = l_curve(U,sm,R,'Tikh');
%     lambda = my_logspace(1e-4,1,100);
%     lc0 = hankraus(Gamma, R, lambda, L, GFbinc);
%     lc0 = optimal_lambda(Gamma, R, lambda, L, GFbinc, c_real);
%     fprintf('Regularization parameter from H-R rule = %f\n',lc0);
    FIMS.GF.cPDF = ncsolve(Gamma, R, lc0, L);
    FIMS.GF.cPDF=FIMS.GF.cPDF(:)./(sum(FIMS.GF.cPDF(:))*(GFbinc(2)-GFbinc(1)));
    g=GFbinc;
    FIMS.GF.R_inv_tk=Gamma*FIMS.GF.cPDF;
    R_inv = FIMS.GF.R_inv_tk(:)/max(FIMS.GF.R_inv_tk(:));
    cPDF = FIMS.GF.cPDF(:);
    FIMS.GF.residue_tk=sum(power((FIMS.GF.R(:)/max(FIMS.GF.R(:))-FIMS.GF.R_inv_tk/max(FIMS.GF.R_inv_tk)),2));
    t_run = toc;
    if action.show_fitting
        figure; 
        subplot(2,1,1); plot(FIMS.GF.Dpbinc,FIMS.GF.R/max(FIMS.GF.R),'o',FIMS.GF.Dpbinc,FIMS.GF.R_inv_tk/max(FIMS.GF.R_inv_tk),'.-');
        title(['Tikhonov: ', 'RH = ',num2str(FIMS.GF.RHmix),'%, Dry size = ', num2str(FIMS.GF.DpDMA*1e9), 'nm']);
        xlabel('Dp (nm)');ylabel('Normalized R')
        subplot(2,1,2); plot(g,cPDF); hold on; plot(g,c_real); legend('{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
        xlabel('g'); ylabel('c(g,Dp1)');
        hold on;
    end
end
%%% tikhonov 2th inversion
if option.tikhonov2_inversion
    GF_channel_bin = varargin{1};
    GFbinc=linspace(0.8,2.0,GF_channel_bin);
    GFbin_bnds=Dp2Dp_set(GFbinc);
    tic
    for i1=1:GF_channel_bin
        GFlobnd=GFbin_bnds(1,i1);
        GFupbnd=GFbin_bnds(3,i1);
        for j1=1:1:FIMS.GF.num 
            Gamma(j1,i1)=quad2d(@(g,D1) GFPDF_Twomey_matrix_double_integrand ...
                (FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                FIMS.Zpbinc,ZpDMA,FIMS.GF.chg,g,D1,T_FIMS,T_DMA,P_FIMS,P_DMA,j1), ...
                GFlobnd,GFupbnd,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
        end
    end
    R = FIMS.GF.R(:)/max(FIMS.GF.R(:));
    order = 2;
    c_real = interp1(FIMS.GF.gNorm, FIMS.GF.c_sim, GFbinc);
    L = full(get_l(length(Gamma(1,:)), order));
    [U, sm, XX, V] = cgsvd(Gamma, L); % Generalized singular value decomposition in compact form
    lc0 = l_curve(U,sm,R,'Tikh');
%     lambda = my_logspace(1e-4,1,100);
%     lc0 = hankraus(Gamma, R, lambda, L, GFbinc);
%     lc0 = optimal_lambda(Gamma, R, lambda, L, GFbinc, c_real);
%     fprintf('Regularization parameter from comparing with true GF-PDF = %f\n',lc0);
    FIMS.GF.cPDF = ncsolve(Gamma, R, lc0, L);
    FIMS.GF.cPDF=FIMS.GF.cPDF(:)./(sum(FIMS.GF.cPDF(:))*(GFbinc(2)-GFbinc(1)));
    g=GFbinc;
    FIMS.GF.R_inv_tk=Gamma*FIMS.GF.cPDF;
    R_inv = FIMS.GF.R_inv_tk(:)/max(FIMS.GF.R_inv_tk(:));
    cPDF = FIMS.GF.cPDF(:);
    FIMS.GF.residue_tk=sum(power((FIMS.GF.R(:)/max(FIMS.GF.R(:))-FIMS.GF.R_inv_tk/max(FIMS.GF.R_inv_tk)),2));
    t_run = toc;
    if action.show_fitting
        figure; 
        subplot(2,1,1); plot(FIMS.GF.Dpbinc,FIMS.GF.R/max(FIMS.GF.R),'o',FIMS.GF.Dpbinc,FIMS.GF.R_inv_tk/max(FIMS.GF.R_inv_tk),'.-');
        title(['Tikhonov: ', 'RH = ',num2str(FIMS.GF.RHmix),'%, Dry size = ', num2str(FIMS.GF.DpDMA*1e9), 'nm']);
        xlabel('Dp (nm)');ylabel('Normalized R')
        subplot(2,1,2); plot(g,cPDF);hold on; plot(g,c_real);legend('{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
        xlabel('g'); ylabel('c(g,Dp1)');
        hold on;
    end
end
%%% twomey inversion
if option.twomey_inversion
    GF_channel_bin = varargin{1};
    GFbinc=linspace(0.8,2.0,GF_channel_bin);
    GFbin_bnds=Dp2Dp_set(GFbinc);
    tic
    for i1=1:GF_channel_bin
        GFlobnd=GFbin_bnds(1,i1);
        GFupbnd=GFbin_bnds(3,i1);
        for j1=1:1:FIMS.GF.num 
            Gamma(j1,i1)=quad2d(@(g,D1) GFPDF_Twomey_matrix_double_integrand ...
                (FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                FIMS.Zpbinc,ZpDMA,FIMS.GF.chg,g,D1,T_FIMS,T_DMA,P_FIMS,P_DMA,j1), ...
                GFlobnd,GFupbnd,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
        end
    end
    c_real = interp1(FIMS.GF.gNorm, FIMS.GF.c_sim, GFbinc);
    FIMS.GF.cPDF=Twomey_inv(Gamma,FIMS.GF.R/max(FIMS.GF.R));
    FIMS.GF.cPDF=FIMS.GF.cPDF(:)./(sum(FIMS.GF.cPDF(:))*(GFbinc(2)-GFbinc(1)));
    g=GFbinc;
    FIMS.GF.R_inv_tm=Gamma*FIMS.GF.cPDF;
    R_inv = FIMS.GF.R_inv_tm(:)/max(FIMS.GF.R_inv_tm(:));
    cPDF = FIMS.GF.cPDF(:);
    FIMS.GF.residue_tm=sum(power((FIMS.GF.R/max(FIMS.GF.R)-FIMS.GF.R_inv_tm/max(FIMS.GF.R_inv_tm)),2));
    t_run = toc;
    if action.show_fitting
        figure; 
        subplot(2,1,1); plot(FIMS.GF.Dpbinc,FIMS.GF.R/max(FIMS.GF.R),'o',FIMS.GF.Dpbinc,FIMS.GF.R_inv_tm/max(FIMS.GF.R_inv_tm),'.-');
        title(['Twomey: ','RH = ',num2str(FIMS.GF.RHmix),'%, Dry size = ', num2str(FIMS.GF.DpDMA*1e9), 'nm']);
        xlabel('Dp (nm)');ylabel('Normalized R')
        subplot(2,1,2); plot(g,cPDF);hold on;plot(g,c_real);legend('{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
        xlabel('g'); ylabel('c(g,Dp1)');
    end
end

%%% lognormal fitting
if option.lognormal_fit
    p = varargin{1};
%     p = 2;
%     if i == 1; p=input('set lognormal num = '); end
    lb=[]; ub=[]; x0=[]; 
    lb_set=[0.0,0.9,1.01]; % fcond, G, sigma
    ub_set=[1.0,2.0,1.50];
    x0_set=[0.5,1.4,1.05];
%     lb_set=[0.0,0.8,1.00]; % lower limit of fcond, G, sigma
%     ub_set=[1.0,2.0,1.5]; % upper limit of fcond, G, sigma
%     x0_set=[1.0,1.1,1.01]; % initial values of fcond, G, sigma
    for j=1:1:p
        lb=horzcat(lb,lb_set);
        ub=horzcat(ub,ub_set);
        x0=horzcat(x0,x0_set);
    end
    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
        'MaxFunctionEvaluations',1000);  
    tic
    fit_to=FIMS.GF.R;
    x=lsqnonlin(@(x) DMAFIMS_log_lsqset(x,fit_to,FIMS.Zpbinc,ZpDMA,factor,T_FIMS,T_DMA,P_FIMS,P_DMA),x0,lb,ub,options);
    FIMS.GF.meas_norm=fit_to/max(fit_to);
    R_F=zeros(1,FIMS.GF.num);
    fcond=zeros(1,p);
    G=zeros(1,p);
    sigma=zeros(1,p);
    for k=1:1:p
        fcond(k)=x(1+(k-1)*3);
        G(k)=x(2+(k-1)*3);
        sigma(k)=x(3+(k-1)*3);
    end
    fcond=fcond/(sum(fcond));
    for j1=1:1:FIMS.GF.num    
        Zpc=FIMS.Zpbinc(j1);
        Dplobnd_FIMS=Zp2Dp(Zpc*1.2,T_FIMS+273.15,P_FIMS,FIMS.GF.chg)*1e9;
        Dpupbnd_FIMS=Zp2Dp(Zpc*0.8,T_FIMS+273.15,P_FIMS,FIMS.GF.chg)*1e9;
        FIMS.GF.R_inv_ml(j1)=FIMS.GF.factor*FIMS.GF.a/FIMS.GF.DMA_beta*quad2d(@(D2,D1) GFPDF_log_double_integrand(FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                FIMS.GF.Zpbinc,ZpDMA,FIMS.GF.chg,D2,D1,fcond,G,sigma,T_FIMS,T_DMA,P_FIMS,P_DMA,j1),Dplobnd_FIMS,Dpupbnd_FIMS,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
    end
    gNorm=[0.8:0.01:2.0];
    FIMS.GF.cPDF=GFPDF_log(gNorm,1,fcond,G,sigma)./gNorm;
    FIMS.GF.residue_ml=sum(power((FIMS.GF.R/max(FIMS.GF.R)-FIMS.GF.R_inv_ml/max(FIMS.GF.R_inv_ml)),2));
    g=gNorm;
    R_inv = FIMS.GF.R_inv_ml(:)/max(FIMS.GF.R_inv_ml(:));
    cPDF = (FIMS.GF.cPDF(:)./(sum(FIMS.GF.cPDF(:))*(g(2)-g(1)))).';
    t_run = toc;
    if action.show_fitting
        figure; 
        subplot(2,1,1); plot(FIMS.GF.Dpbinc,FIMS.GF.R/max(FIMS.GF.R),'o',FIMS.GF.Dpbinc,FIMS.GF.R_inv_ml/max(FIMS.GF.R_inv_ml),'.-');
        title(['ML: ','RH = ',num2str(FIMS.GF.RHmix),'%, Dry size = ', num2str(FIMS.GF.DpDMA*1e9), 'nm']);
        xlabel('Dp (nm)');ylabel('Normalized R')
        subplot(2,1,2); plot(g,cPDF);hold on; plot(g,c_real);legend('{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
        xlabel('g'); ylabel('c(g,Dp1)');
    end
end
%%% piecewise linear fitting   
if option.piecewise_fit
%     if i == 1; q=input('set piecewise segments = '); end
    q = varargin{1};
%     q = 20;
    range=linspace(0.8,2.0,q+2);
    lb=zeros(1,q);
    ub=30*ones(1,q);
    x0=1*ones(1,q);
    x0=interp1(FIMS.GF.Dpbinc*1e-9/FIMS.GF.DpDMA,FIMS.GF.R/max(FIMS.GF.R),range(2:end-1));
    x0(isnan(x0))=0;
    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
        'MaxFunctionEvaluations',1000);
    tic
    fit_to=FIMS.GF.R;
    x=lsqnonlin(@(x) DMAFIMS_pie_lsqset(x,fit_to,FIMS.Zpbinc,ZpDMA,factor,T_FIMS,T_DMA,P_FIMS,P_DMA,range),x0,lb,ub,options);
    FIMS.GF.meas_norm=fit_to/max(fit_to);
    R_F=zeros(1,FIMS.GF.num);

    for j1=1:1:FIMS.GF.num    
        Zpc=FIMS.Zpbinc(j1);
        Dplobnd_FIMS=Zp2Dp(Zpc*1.2,T_FIMS+273.15,P_FIMS,FIMS.GF.chg)*1e9;
        Dpupbnd_FIMS=Zp2Dp(Zpc*0.8,T_FIMS+273.15,P_FIMS,FIMS.GF.chg)*1e9;
        FIMS.GF.R_inv_pl(j1)=FIMS.GF.factor*FIMS.GF.a/FIMS.GF.DMA_beta*quad2d(@(D2,D1) GFPDF_pie_double_integrand(FIMS.normZp_full_range,FIMS.DMAomega,FIMS.omega_norm, ...
                    FIMS.GF.Zpbinc,ZpDMA,FIMS.GF.chg,D2,D1,x,range,T_FIMS,T_DMA,P_FIMS,P_DMA,j1),Dplobnd_FIMS,Dpupbnd_FIMS,FIMS.GF.Dplobnd_DMA,FIMS.GF.Dpupbnd_DMA);
    end
    norm_x=x/(sum(x)*(range(2)-range(1))); %normalized by the integral area
    gNorm=[0.8:0.01:2.0];
    FIMS.GF.cPDF=GFPDF_pie(gNorm,1,norm_x,range);
    FIMS.GF.residue_pl=sum(power((FIMS.GF.R/max(FIMS.GF.R)-FIMS.GF.R_inv_pl/max(FIMS.GF.R_inv_pl)),2));
    g=gNorm;
    R_inv = FIMS.GF.R_inv_pl(:)/max(FIMS.GF.R_inv_pl(:));
    cPDF = (FIMS.GF.cPDF(:)./(sum(FIMS.GF.cPDF(:))*(g(2)-g(1)))).';
    t_run = toc;
    if action.show_fitting
        figure; 
        subplot(2,1,1); plot(FIMS.GF.Dpbinc,FIMS.GF.R/max(FIMS.GF.R),'o',FIMS.GF.Dpbinc,FIMS.GF.R_inv_pl/max(FIMS.GF.R_inv_pl),'.-');
        title(['PL: ','RH = ',num2str(FIMS.GF.RHmix),'%, Dry size = ', num2str(FIMS.GF.DpDMA*1e9), 'nm']);
        xlabel('Dp (nm)');ylabel('Normalized R')
        subplot(2,1,2); plot(g,cPDF);hold on; plot(g,c_real);legend('{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
        xlabel('g'); ylabel('c(g,Dp1)');
    end
end

end