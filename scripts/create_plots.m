% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function create_plots(test)
% PLOTTING: Plot convergence trends for statistics and sparsity boxplots

% load parameters
d        = test.d;
order    = test.order;
func_str = test.func_str;
N        = test.N;

% load data
load(['../results/' func_str '_d' num2str(d) '_ord' num2str(order)]);

% load sample data
load(['../results/' func_str '_d' num2str(d) '_samples'], 'MC_out', 'QMC_out');

%% Plot results

% Declare mean and standard error functions
Mout  = @(c, field) cell2mat(cellfun(@(x) x.(field), c, 'uniformoutput',false));
pmean = @(c, field) mean(Mout(c, field),2);
pste  = @(c, field) 1.96*std(Mout(c, field),[],2)/sqrt(size(c,2));

% define colors
gColor = [85;170;170]/255;
bColor = [60;60;230]/255;
rColor = [200;0;0]/255;
pColor = [170;0;170]/255;
yColor = [225;125;0]/255;

% plot error in mean
figure
hold on
grid on
errorbar(N, pmean(OMP_out,'MeanE'),  pste(OMP_out,'MeanE'),  '-o', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP');
errorbar(N, pmean(OMPN_out,'MeanE'), pste(OMPN_out,'MeanE'), '-d', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMPN');
errorbar(N, pmean(BPDN_out,'MeanE'), pste(BPDN_out,'MeanE'), '-o', 'Color', bColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BPDN');
errorbar(N, pmean(RGA_out,'MeanE'),  pste(RGA_out,'MeanE'),  '-o', 'Color', rColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'RGA');
errorbar(N, pmean(MC_out,'MeanE'),   pste(MC_out,'MeanE'),   '-.o', 'Color', pColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'MC');
errorbar(N, pmean(QMC_out,'MeanE'),  pste(QMC_out,'MeanE'),  '-.o', 'Color', yColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'QMC');
xlabel('Number of Samples, N')
ylabel('Relative Mean Error')
legend('show', 'location', 'southwest')
xlim([0,max(N)])
set(gca,'YScale','log')
set(gca,'YMinorGrid','Off')
hold off
print('-depsc',['../results/' func_str '_MeanErr_vs_N'])

% plot error in standard deviation
figure
hold on
grid on
errorbar(N, pmean(OMP_out,'StdE'),  pste(OMP_out,'StdE'),  '-o', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP');
errorbar(N, pmean(OMPN_out,'StdE'), pste(OMPN_out,'StdE'), '-d', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMPN');
errorbar(N, pmean(BPDN_out,'StdE'), pste(BPDN_out,'StdE'), '-o', 'Color', bColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BPDN');
errorbar(N, pmean(RGA_out,'StdE'),  pste(RGA_out,'StdE'),  '-o', 'Color', rColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'RGA');
errorbar(N, pmean(MC_out,'StdE'),   pste(MC_out,'StdE'),   '-.o', 'Color', pColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'MC');
errorbar(N, pmean(QMC_out,'StdE'),  pste(QMC_out,'StdE'),  '-.o', 'Color', yColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'QMC');
xlabel('Number of Samples, N')
ylabel('Relative Standard Deviation Error')
legend('show', 'location', 'southwest')
xlim([0,max(N)])
set(gca,'YScale','log')
set(gca,'YMinorGrid','Off')
hold off
print('-depsc',['../results/' func_str '_StdErr_vs_N'])

% plot test error
figure
hold on
grid on
errorbar(N, pmean(OMP_out,'TestE'),  pste(OMP_out,'TestE'),  '-o', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP');
errorbar(N, pmean(OMPN_out,'TestE'), pste(OMPN_out,'TestE'), '-d', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMPN');
errorbar(N, pmean(BPDN_out,'TestE'), pste(BPDN_out,'TestE'), '-o', 'Color', bColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BPDN');
errorbar(N, pmean(RGA_out,'TestE'),  pste(RGA_out,'TestE'),  '-o', 'Color', rColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'RGA');
xlabel('Number of Samples, N')
ylabel('Test Set Error')
legend('show', 'location', 'southwest')
xlim([0,max(N)])
set(gca,'YScale','log')
set(gca,'YMinorGrid','Off')
hold off
print('-depsc',['../results/' func_str '_TestErr_vs_N'])

% plot sparsity
figure
hold on
grid on
errorbar(N, pmean(OMP_out,'spars'),  pste(OMP_out,'spars'),  '-o', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP');
errorbar(N, pmean(OMPN_out,'spars'), pste(OMPN_out,'spars'), '-d', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMPN');
errorbar(N, pmean(BPDN_out,'spars'), pste(BPDN_out,'spars'), '-o', 'Color', bColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BPDN');
errorbar(N, pmean(RGA_out,'spars'),  pste(RGA_out,'spars'),  '-o', 'Color', rColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'RGA');
xlabel('Number of Samples, N')
ylabel('$\|\mathbf{c}\|_{0}$')
legend('show', 'location', 'southeast')
xlim([0,max(N)])
set(gca,'YMinorGrid','Off')
hold off
print('-depsc',['../results/' func_str '_Spars_vs_N'])

% plot runtime
figure
hold on
grid on
errorbar(N, pmean(OMP_out,'time'),  pste(OMP_out,'time'),  '-o', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP');
errorbar(N, pmean(OMPN_out,'time'), pste(OMPN_out,'time'), '-d', 'Color', gColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMPN');
errorbar(N, pmean(BPDN_out,'time'), pste(BPDN_out,'time'), '-o', 'Color', bColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BPDN');
errorbar(N, pmean(RGA_out,'time'),  pste(RGA_out,'time'),  '-o', 'Color', rColor, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'RGA');
xlabel('Number of Samples, N')
ylabel('Runtime, T')
legend('show', 'location', 'southeast')
xlim([0,max(N)])
set(gca,'YScale','log')
set(gca,'YMinorGrid','Off')
hold off
print('-depsc',['../results/' func_str '_Runtime_vs_N'])

end

% -- END OF FILE --
