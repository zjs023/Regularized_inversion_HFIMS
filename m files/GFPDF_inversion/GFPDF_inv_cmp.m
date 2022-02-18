function GFPDF_inv_cmp(option)

    option.residue_cmp = 1;

    action.show_residue_cmp = 1;
    action.show_residue_cmp_stats = 1;


    if option.residue_cmp
        option.lognormal_fit = 0;
        option.piecewise_fit = 0;
        option.twomey_inversion = 0; 
        option.tikhonov0_inversion = 0;
        option.tikhonov1_inversion = 0;
        option.tikhonov2_inversion = 0;
        option.lsq_inversion = 0;
        
        option.lsq_inversion = 1;
        GFPDF_lsq = GFPDF_Poisson_cmp(option);
        option.lsq_inversion = 0;
        
        option.tikhonov0_inversion = 1; 
        GFPDF_tk0 = GFPDF_Poisson_cmp(option);
        option.tikhonov0_inversion = 0;
   
        option.tikhonov1_inversion = 1; 
        GFPDF_tk1 = GFPDF_Poisson_cmp(option);
        option.tikhonov1_inversion = 0; 
        
        option.tikhonov2_inversion = 1;
        GFPDF_tk2 = GFPDF_Poisson_cmp(option);
        option.tikhonov2_inversion = 0;
  
        option.twomey_inversion = 1;
        GFPDF_tm = GFPDF_Poisson_cmp(option);
        option.twomey_inversion = 0;

        option.lognormal_fit = 1;
        GFPDF_ml = GFPDF_Poisson_cmp(option);
        option.lognormal_fit = 0;

        option.piecewise_fit = 1;
        GFPDF_pl = GFPDF_Poisson_cmp(option);
        option.piecewise_fit = 0;
    end

    if option.residue_cmp && action.show_residue_cmp
        figure;
        bar_series = [GFPDF_lsq.residue_R_mean'; GFPDF_tk0.residue_R_mean'; GFPDF_tk1.residue_R_mean'; GFPDF_tk2.residue_R_mean'; GFPDF_tm.residue_R_mean'];
        bar_error = [GFPDF_lsq.residue_R_std'; GFPDF_tk0.residue_R_std'; GFPDF_tk1.residue_R_std'; GFPDF_tk2.residue_R_std'; GFPDF_tm.residue_R_std'];
        plot_hist_with_errbar(bar_series,bar_error)
        legend('GFPDF-1','GFPDF-2','GFPDF-3');
        xticklabels({'LSQ','Tik 0^{th}','Tik 1^{st}','Tik 2^{nd}','Twomey'}); ylabel('{\it\chi}^{2}');

        figure;
        bar_series = [GFPDF_lsq.residue_c_mean'; GFPDF_tk0.residue_c_mean'; GFPDF_tk1.residue_c_mean'; GFPDF_tk2.residue_c_mean'; GFPDF_tm.residue_c_mean'];
        bar_error = [GFPDF_lsq.residue_c_std'; GFPDF_tk0.residue_c_std'; GFPDF_tk1.residue_c_std'; GFPDF_tk2.residue_c_std'; GFPDF_tm.residue_c_std'];
        plot_hist_with_errbar(bar_series,bar_error)
        legend('GFPDF-1','GFPDF-2','GFPDF-3');
        xticklabels({'LSQ','Tik 0^{th}','Tik 1^{st}','Tik 2^{nd}','Twomey'}); ylabel('\gamma^{2}');
        
        figure;
        bar_series = [GFPDF_lsq.roughness_mean'; GFPDF_tk0.roughness_mean'; GFPDF_tk1.roughness_mean'; GFPDF_tk2.roughness_mean'; GFPDF_tm.roughness_mean'];
        bar_error = [GFPDF_lsq.roughness_std'; GFPDF_tk0.roughness_std'; GFPDF_tk1.roughness_std'; GFPDF_tk2.roughness_std'; GFPDF_tm.roughness_std'];
        plot_hist_with_errbar(bar_series,bar_error)
        legend('GFPDF-1','GFPDF-2','GFPDF-3');
        xticklabels({'LSQ','Tik 0^{th}','Tik 1^{st}','Tik 2^{nd}','Twomey'}); ylabel('\it\xi');
        
        figure;
        bar_series = [GFPDF_lsq.t_run_mean'./GFPDF_lsq.t_run_mean'; GFPDF_tk0.t_run_mean'./GFPDF_lsq.t_run_mean'; GFPDF_tk1.t_run_mean'./GFPDF_lsq.t_run_mean'; GFPDF_tk2.t_run_mean'; GFPDF_tm.t_run_mean'./GFPDF_lsq.t_run_mean'];
        bar_error = [GFPDF_lsq.t_run_std'; GFPDF_tk0.t_run_std'; GFPDF_tk1.t_run_std'; GFPDF_tk2.t_run_std'; GFPDF_tm.t_run_std'];
        plot_hist_with_errbar(bar_series,bar_error)
        legend('GFPDF-1','GFPDF-2','GFPDF-3');
        xticklabels({'LSQ','Tik 0^{th}','Tik 1^{st}','Tik 2^{nd}','Twomey'}); ylabel('Computing time [s]');
    end
    if option.residue_cmp && action.show_residue_cmp_stats
%         figure;boxplot([GFPDF_lsq.residue_R',GFPDF_tk0.residue_R',GFPDF_tk1.residue_R',GFPDF_tk2.residue_R',GFPDF_tm.residue_R'],'symbol', '');ylabel('{\it\chi}^{2}');
%         figure;boxplot([GFPDF_lsq.residue_c',GFPDF_tk0.residue_c',GFPDF_tk1.residue_c',GFPDF_tk2.residue_c',GFPDF_tm.residue_c'],'symbol', '');ylabel('\gamma^{2}');
%         figure;boxplot([GFPDF_lsq.roughness',GFPDF_tk0.roughness',GFPDF_tk1.roughness',GFPDF_tk2.roughness',GFPDF_tm.roughness'],'symbol', '');ylabel('\it\xi');
        figure;i=3;boxplot([GFPDF_lsq.residue_R{i}',GFPDF_tk0.residue_R{i}',GFPDF_tk1.residue_R{i}',GFPDF_tk2.residue_R{i}',GFPDF_tm.residue_R{i}'], 'Colors', 'r', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5]+0.15, 'Width', 0.3); set(gca,'xtick',[]);
        hold on;i=2;boxplot([GFPDF_lsq.residue_R{i}',GFPDF_tk0.residue_R{i}',GFPDF_tk1.residue_R{i}',GFPDF_tk2.residue_R{i}',GFPDF_tm.residue_R{i}'], 'Colors', 'b', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5], 'Width', 0.3); set(gca,'xtick',[]);
        hold on;i=1;boxplot([GFPDF_lsq.residue_R{i}',GFPDF_tk0.residue_R{i}',GFPDF_tk1.residue_R{i}',GFPDF_tk2.residue_R{i}',GFPDF_tm.residue_R{i}'], 'Colors', 'k', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5]-0.15, 'Width', 0.3); set(gca,'xtick',[]); 
        ylabel('{\it\chi}^{2}'); xticks([1 2 3 4 5]); xticklabels({'LSQ','Tik 0^{th}','Tik 1^{st}','Tik 2^{nd}','Twomey'}); set(gca, 'YGrid', 'on', 'XGrid', 'off'); set(gca,'TickLabelInterpreter', 'tex');
        
        figure;i=3;boxplot([GFPDF_lsq.residue_c{i}',GFPDF_tk0.residue_c{i}',GFPDF_tk1.residue_c{i}',GFPDF_tk2.residue_c{i}',GFPDF_tm.residue_c{i}'], 'Colors', 'r', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5]+0.15, 'Width', 0.3); set(gca,'xtick',[]);
        hold on;i=2;boxplot([GFPDF_lsq.residue_c{i}',GFPDF_tk0.residue_c{i}',GFPDF_tk1.residue_c{i}',GFPDF_tk2.residue_c{i}',GFPDF_tm.residue_c{i}'], 'Colors', 'b', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5], 'Width', 0.3); set(gca,'xtick',[]);
        hold on;i=1;boxplot([GFPDF_lsq.residue_c{i}',GFPDF_tk0.residue_c{i}',GFPDF_tk1.residue_c{i}',GFPDF_tk2.residue_c{i}',GFPDF_tm.residue_c{i}'], 'Colors', 'k', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5]-0.15, 'Width', 0.3); set(gca,'xtick',[]); 
        ylabel('\gamma^{2}'); xticks([1 2 3 4 5]); xticklabels({'LSQ','Tik 0^{th}','Tik 1^{st}','Tik 2^{nd}','Twomey'}); set(gca, 'YGrid', 'on', 'XGrid', 'off'); set(gca,'TickLabelInterpreter', 'tex');
        
        figure;i=3;boxplot([GFPDF_lsq.roughness{i}',GFPDF_tk0.roughness{i}',GFPDF_tk1.roughness{i}',GFPDF_tk2.roughness{i}',GFPDF_tm.roughness{i}'], 'Colors', 'r', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5]+0.15, 'Width', 0.3); set(gca,'xtick',[]);
        hold on;i=2;boxplot([GFPDF_lsq.roughness{i}',GFPDF_tk0.roughness{i}',GFPDF_tk1.roughness{i}',GFPDF_tk2.roughness{i}',GFPDF_tm.roughness{i}'], 'Colors', 'b', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5], 'Width', 0.3); set(gca,'xtick',[]);
        hold on;i=1;boxplot([GFPDF_lsq.roughness{i}',GFPDF_tk0.roughness{i}',GFPDF_tk1.roughness{i}',GFPDF_tk2.roughness{i}',GFPDF_tm.roughness{i}'], 'Colors', 'k', 'Symbol', '', 'PlotStyle','compact','BoxStyle', 'filled','Position', [1:5]-0.15, 'Width', 0.3); set(gca,'xtick',[]); 
        ylabel('\it\xi'); xticks([1 2 3 4 5]); xticklabels({'LSQ','Tik 0^{th}','Tik 1^{st}','Tik 2^{nd}','Twomey'}); set(gca, 'YGrid', 'on', 'XGrid', 'off'); set(gca,'TickLabelInterpreter', 'tex');
    end
end
function plot_hist_with_errbar(bar_series,bar_error)
    bar(bar_series, 'grouped'); hold on;
    ngroups = size(bar_series, 1);
    nbars = size(bar_series, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%         errorbar(x, bar_series(:,i), bar_error(:,i), 'k', 'linestyle', 'none');
    end
    hold off
end
