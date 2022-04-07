function GFPDF = GFPDF_Poisson_cmp(option)

if nargin < 1
    option.lognormal_fit = 0;
    option.piecewise_fit = 0;
    option.twomey_inversion = 0; 
    option.tikhonov0_inversion = 0;
    option.tikhonov1_inversion = 0;
    option.tikhonov2_inversion = 0;
    option.lsq_inversion = 0;
end

action.generate_new_data = 1;
action.show_fitting = 1;
action.show_L_curve = 0;

global FIMS
flt_num = 332;
image_bins = 30;
N_tot = 100;
% GFbins = [10,15,20,25,30,35,40,45,50];
GFbins = 20;
% MLmodes = [1,2,3,4];
MLmodes = [1,2,3];
% PLbins = [10,15,20,25,30];
PLbins = [15,15,15];
N_trial=10; % run 10 trials and take average

var2str = @(~) inputname(1);
var_str={var2str(GFbins),var2str(MLmodes),var2str(PLbins)};
if option.lsq_inversion
    var_name = var_str{1};
end
if option.twomey_inversion
    var_name = var_str{1};
end
if option.tikhonov0_inversion
    var_name = var_str{1};
end
if option.tikhonov1_inversion
    var_name = var_str{1};
end
if option.tikhonov2_inversion
    var_name = var_str{1};
end
if option.lognormal_fit
    var_name = var_str{2};
end
if option.piecewise_fit
    var_name = var_str{3};
end
var=eval(var_name);

first_run = 0;
load('.\data\FIMS\FIMS_f332_r01.mat');

for i=1:1:3         % ith pre-defined GF
	if strcmp(var_name,'MLmodes') 
		var=eval(var_name);
		var = var(i);
	end
	if strcmp(var_name,'PLbins') 
		var=eval(var_name);
		var = var(i);
	end
	gNorm=[0.8:0.01:2.0];
	[cPDF_sim,R_sim_raw]=Forward_Response(i);  % ith pre-defined GF
	R_sim_raw=R_sim_raw/sum(R_sim_raw)*N_tot;
	FIMS.GF.R_sim = R_sim_raw;
	FIMS.GF.gNorm = gNorm;
	FIMS.GF.c_sim = cPDF_sim;
	for k=1:length(var) % kth GF bin number/ ML mode number/ PL seg bin number
		R_sim_tot=[];
		R_inv_tot=[];
		cPDF_inv_tot=[];
		t_run_tot=[];
		fprintf('\nEvaluating %d trails.\n',N_trial);
		for t=1:N_trial
			fprintf('%d.\n',t);
			if action.generate_new_data
				sig_gaus = 0.05;
				R_raw_Gaus = R_sim_raw.*(1+sig_gaus*randn(size(R_sim_raw)));
				R_poiss = poissrnd(R_raw_Gaus);
				R_synth(t,:) = R_poiss;

				sig_poiss = sqrt(R_poiss+power(sig_gaus*R_poiss,2));
				FIMS.GF.Sigma = sig_poiss;     
			else
				R_synth=FIMS.GF.R_synth{i};
			end
			FIMS.GF.R=R_synth(t,:);
			R_sim_tot = horzcat(R_sim_tot, FIMS.GF.R(:)./max(FIMS.GF.R));
			[g,cPDF_inv,R_inv,t_run]=DMAFIMS_GFPDF_Poisson(FIMS,option,var(k));
			cPDF_inv = interp1(g,cPDF_inv',gNorm)';
			R_inv_tot = horzcat(R_inv_tot, R_inv);
			cPDF_inv_tot = horzcat(cPDF_inv_tot,cPDF_inv);
			t_run_tot = horzcat(t_run_tot,t_run);
			if 0
				figure;
				subplot(2,1,1); plot(FIMS.GF.Dpbinc,FIMS.GF.R(:)./max(FIMS.GF.R),'o',FIMS.GF.Dpbinc,R_inv/max(R_inv),'.-');
				title('Residue=',nansum(power((FIMS.GF.R(:)./max(FIMS.GF.R)-R_inv/max(R_inv)),2)));
				xlabel('{\itD}p (nm)');ylabel('Normalized {\itR}');legend('{\itR}_{sim}','{\itR}_{inv}');
				subplot(2,1,2);plot(gNorm,cPDF_sim,gNorm,cPDF_inv');
				title('GFPDF err=',nansum(power((cPDF_sim-cPDF_inv'),2))/length(cPDF_inv'));
				xlabel('\itg'); ylabel('{\itc}({\itg},{\itD}p_{1})');legend('{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
				pause
			end
		end
		FIMS.GF.R_synth{i}=R_synth;
		R_sim=nanmean(R_sim_tot,2);
		R_inv=nanmean(R_inv_tot,2);
		cPDF_inv=nanmean(cPDF_inv_tot,2)';
		cPDF_sim_tot = repmat(cPDF_sim',1,size(cPDF_inv_tot,2));
		residue_R=nansum(power((R_sim_tot-R_inv_tot),2));
		residue_c=nansum(power((cPDF_sim_tot-cPDF_inv_tot),2))/length(cPDF_inv);
		roughness=nansum(abs(2*cPDF_inv_tot(2:end-1,:)-cPDF_inv_tot(1:end-2,:)-cPDF_inv_tot(3:end,:)));
		avg.residue_R(i)=nansum(power((R_sim-R_inv),2));
		avg.residue_c(i)=nansum(power((cPDF_sim-cPDF_inv),2))/length(cPDF_inv);
		avg.roughness(i)=nansum(abs(2*cPDF_inv(2:end-1)-cPDF_inv(1:end-2)-cPDF_inv(3:end)));

		GFPDF.avg=avg;
		GFPDF.residue_R{i}=residue_R;
		GFPDF.residue_c{i}=residue_c;
		GFPDF.roughness{i}=roughness;
		GFPDF.residue_R_mean(i,k)=nanmean(residue_R);
		GFPDF.residue_R_std(i,k)=nanstd(residue_R);
		GFPDF.residue_c_mean(i,k)=nanmean(residue_c);
		GFPDF.residue_c_std(i,k)=nanstd(residue_c);
		GFPDF.roughness_mean(i,k)=nanmean(roughness);
		GFPDF.roughness_std(i,k)=nanstd(roughness);
		GFPDF.t_run_mean(i,k)=nanmean(t_run_tot);
		GFPDF.t_run_std(i,k)=nanstd(t_run_tot);
		if action.show_fitting
			if 0
				figure;
				for i=1:size(R_inv_tot,2) 
					subplot(2,1,1); plot(FIMS.GF.Dpbinc,R_sim_tot(:,i)/max(R_sim_tot(:,i)),'o',FIMS.GF.Dpbinc,R_inv_tot(:,i)/max(R_inv_tot(:,i)),'.-');
					title('Residue=',nansum(power((R_sim_tot(:,i)-R_inv_tot(:,i)),2)));
					subplot(2,1,2);plot(gNorm,cPDF_sim_tot(:,i),gNorm,cPDF_inv_tot(:,i));
					title('GFPDF err=',nansum(power((cPDF_sim_tot(:,i)-cPDF_inv_tot(:,i)),2))/length(cPDF_inv_tot(:,i)));
					pause
				end
			end

			figure; 
			subplot(2,1,1); plot(FIMS.GF.Dpbinc,R_sim/max(R_sim),'o',FIMS.GF.Dpbinc,R_inv/max(R_inv),'.-');
			title(['Image bins = ',num2str(image_bins),', ',var_name, ' = ',num2str(var(k)),', Dry size = ', num2str(FIMS.GF.DpDMA*1e9), 'nm']);
			xlabel('{\itD}p (nm)');ylabel('Normalized {\itR}');legend('{\itR}_{sim}','{\itR}_{inv}');
			subplot(2,1,2); 
			confplot(gNorm,cPDF_inv,nanstd(cPDF_inv_tot,0,2)',nanstd(cPDF_inv_tot,0,2)','Color',[1 0 0],'LineWidth',2);
			hold on; plot(gNorm,cPDF_sim);
			xlabel('\itg'); ylabel('{\itc}({\itg},{\itD}p_{1})');legend('{\itc}_{PDF\_inv} \pm 1*\sigma','{\itc}_{PDF\_inv}','{\itc}_{PDF\_sim}');
		end
	end
        
    GFPDF.residue_R_mean(GFPDF.residue_R_mean==0)=nan;
    GFPDF.residue_R_std(GFPDF.residue_R_std==0)=nan;
    GFPDF.residue_c_mean(GFPDF.residue_c_mean==0)=nan;
    GFPDF.residue_c_std(GFPDF.residue_c_std==0)=nan;
    GFPDF.roughness_mean(GFPDF.roughness_mean==0)=nan;
    GFPDF.roughness_std(GFPDF.roughness_std==0)=nan;
    
    if action.show_L_curve
        figure; hold on;
%         edges = 0:0.6:2.4;
%         [residue_c_grp,residue_c_grp_bnd]=discretize(GFPDF.residue_c_mean(i,:),edges);
        [residue_c_grp,residue_c_grp_bnd]=discretize(GFPDF.residue_c_mean(i,:),4);
        bubsizes = unique(residue_c_grp)';
        legentry=cell(size(bubsizes));
        for idx = 1:numel(bubsizes)
            bubleg(idx) = plot(nan, nan, 'ko','MarkerSize', sqrt(25*bubsizes(idx)), 'MarkerFaceColor', 'k');
            legentry{idx} = strcat(num2str(residue_c_grp_bnd(bubsizes(idx))),'~',num2str(residue_c_grp_bnd(bubsizes(idx)+1)));
        end    
        scatter(GFPDF.residue_R_mean(i,:),GFPDF.roughness_mean(i,:),100.*residue_c_grp,var,'filled'); 
        hold on;
        errorbar(GFPDF.residue_R_mean(i,:),GFPDF.roughness_mean(i,:),GFPDF.roughness_std(i,:),GFPDF.roughness_std(i,:),...
            GFPDF.residue_R_std(i,:),GFPDF.residue_R_std(i,:),'LineStyle','none');
        leg=legend(legentry);title(leg,'\gamma^{2}')
        cbh=colorbar; colormap(hsv(length(var))); ylabel(cbh,var_name);
        set(cbh,'YTick',var);
        set(cbh,'YTickLabel',num2cell(var));
        xlabel('{\it\chi}^{2}');ylabel('\it\xi');
        title(['Pre-defined GF-PDF ',num2str(i)]);
    end
end


           
            