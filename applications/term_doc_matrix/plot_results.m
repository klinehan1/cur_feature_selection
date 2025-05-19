% plot cur_perf_dtm.m results

load accuracy_results.mat
load time_results.mat

%%

plot(nc_nr,sf_err, '-', 'Marker', ".", 'MarkerSize', 17, 'LineWidth',2);
hold on; grid on;
errorbar(nc_nr, ls_r_err, ls_r_sd, '-', 'LineWidth',2);
plot(nc_nr,ls_d_err, '-d', 'LineWidth',2);
plot(nc_nr,deim_err, '-o', 'LineWidth',2);
plot(nc_nr,qr_err, '-^', 'LineWidth',2);
plot(nc_nr,svd_err, '-s', 'LineWidth',2);
legend('sf','ls-r','ls-d','deim','qr', 'svd', 'Location','southwest');
title('CUR Accuracy on Document Term Matrix');
xlabel("Number of Selected Columns/Rows");
ylabel("Relative Error");
ax = gca;
ax.Children = ax.Children([6 1 2 3 4 5]);

print('dtm_acc.eps', '-depsc')

%%

semilogy(nc_nr,sf_t, '-', 'Marker', ".", 'MarkerSize', 17, 'LineWidth',2);
hold on; grid on;
errorbar(nc_nr, ls_r_t, ls_r_t_sd, '-', 'LineWidth',2);
semilogy(nc_nr,ls_d_t, '-d', 'LineWidth',2);
semilogy(nc_nr,deim_t, '-o', 'LineWidth',2);
semilogy(nc_nr,qr_t, '-^', 'LineWidth',2);
semilogy(nc_nr,svd_t, '-s', 'LineWidth',2);
legend('sf','ls-r','ls-d','deim','qr', 'svd', 'Location', 'southeast','NumColumns',2);
title('CUR Time on Document Term Matrix');
xlabel("Number of Selected Columns/Rows");
ylabel("Time (seconds)");
ax = gca;
ax.Children = ax.Children([6 1 2 3 4 5]);

print('dtm_time.eps', '-depsc')