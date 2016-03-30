%Explore the bioactivity (ccq75) of all DMSO's in A2 across cell lines

%% Replicate reproducibility for DMSO's across cell lines
dmso = sig_info('{"pert_iname":"DMSO"}');

core_lines = parse_grp('/cmap/data/vdb/spaces/lincs_core_lines.grp');

figure;
for ii = 1:numel(core_lines)
    idx = find(strcmp({dmso.cell_id},core_lines(ii)));
    data = {dmso(idx).distil_cc_q75};
    
    %Convert to matrix and remove -666
    data = cell2mat(data);
    data(data == -666) = [];
    
    [~,f,x,cdf] = kde(data,2^6,0,1);
    plot(x,cdf, ...
        'DisplayName', strcat(core_lines{ii},'_n=',num2str(numel(data))))
    hold on

end

l = legend('show');
set(l,'Interpreter','none')
xlabel('cc_q75','Interpreter','none')
ylabel('cdf')
title(sprintf('Replicate Reproducibility for DMSO \n A2 signatures'))

%% Compare reproducibility of DMSO's to cps.
cp = sig_info('{"pert_type":"trt_cp"}','fields',{'cell_id','distil_cc_q75'});
untrt = sig_info('{"pert_type":"ctl_untrt"}');

dmso_ccq75 = cell2mat({dmso.distil_cc_q75});
dmso_ccq75(dmso_ccq75 <= 0) = [];
cp_ccq75 = cell2mat({cp.distil_cc_q75});
cp_ccq75(cp_ccq75 <= 0) = [];
untrt_ccq75 = cell2mat({untrt.distil_cc_q75});
untrt_ccq75(untrt_ccq75 <= 0) = [];

figure;
mesh_size = 2^6;
[~,f_dmso,x_dmso,cdf_dmso] = kde(dmso_ccq75,mesh_size,0,1);
[~,f_cp,x_cp,cdf_cp] = kde(cp_ccq75,mesh_size,0,1);
[~,f_untrt,x_untrt,cdf_untrt] = kde(untrt_ccq75,mesh_size,0,1);
plot(x_dmso,f_dmso,'DisplayName','DMSO')
hold on
plot(x_cp,f_cp,'DisplayName','cp')
hold on
plot(x_untrt,f_untrt,'DisplayName','untrt')
xlabel('cc_q75','Interpreter','none')
ylabel('density')
legend('show')
title('CP vs UNTRT vs DMSO reproducibility')

figure;
plot(x_dmso,cdf_dmso,'DisplayName','DMSO')
hold on
plot(x_cp,cdf_cp,'DisplayName','cp')
hold on
plot(x_cp,cdf_cp,'DisplayName','untrt')
xlabel('cc_q75','Interpreter','none')
ylabel('cdf')
legend('show')
title('CP vs UNTRT vs DMSO reproducibility')
