% plotting analytical results presented in Figure 3 of the manuscript- 
% "Coupling between DNA replication, segregation and the onset of 
% constriction in Escherichia coli". Note that this program plots for 
% 1 growth medium on a single run. For plotting all 4 growth media shown 
% in Figure 3, change the "med" parameter in the program manually and 
% execute each individual section. To overlay different growth media 
% statistics on same figures, have the right figure active and in the
% plotting sections comment out the figure command. 

clear

% Initializations
x_all = 1:-0.1:0; % fraction of C from termination at which constriction 
% control is found. 1 corresponds to initiation and 0 to termination. 
[num,txt,raw] = xlsread('D:\cell division\exp data\jaan mannik new data\jaan paper figures\Figure 2 complete\params.xlsx');
% The above excel file stores different parameters such as mean C period,
% standard deviation of C period, mean time of initiation relative to birth
% etc., used in this program. These parameters are derived from
% experimental data.
med= 1; % media chosen. 1 corresponds to slowest growth medium. In 
% Figure 3, we plot results from 4 different media. Hence this paramater 
% can go from 1 to 4. 
exp_t= 3; % in minutes, the time for which DnaN marker is hypothesized to 
% remain bound to DNA after replication termination.

% Outputs
cvtnt= zeros(1, length(x_all)); % CV of Tn-Trt. (Tn= timing of onset of constriction, Trt= timing of termination)
corttn= zeros(1, length(x_all)); % Correlation coefficient (CC) between (Trt, Tn).
Pttn= zeros(1, length(x_all)); % Slope of best linear fit of Tn vs Trt plot.

% Output statistics obtained from experimental data (dotted line in Figure
% 3B-3D)
cv_exp = num(med,12); % CV of Tn-Trt.
r_exp= num(med, 15); % Correlation coefficient (CC) between (Trt, Tn).
P_exp= num(med, 18); % Slope of best linear fit of Tn vs Trt plot.

for i=1:length(x_all)
    
mu_C = num(med,2)-exp_t; % mean C period
sig_C = sqrt((num(med,3)*num(med,2))^2-exp_t^2); % standard deviation of C
mu_ti= num(med,6); % mean Tri value
sig_ti= num(med,7); % standard deviation of Tri
x=x_all(i); % current value of x_all
inv_r = (num(med,5)+exp_t+mu_C*(x)); % rate of Poisson process (r)

cvtnt(i) = sqrt(inv_r^2+sig_C^2*x^2+exp_t^2)/(inv_r-(mu_C*x+exp_t)); % CV using Equation 1 of manuscript
sig_rt= sqrt(sig_ti^2+sig_C^2+exp_t^2); sig_n = sqrt(sig_ti^2+sig_C^2*(1-x)^2+inv_r^2); % variances using Equations 14, 15 of manuscript
co_rn= sig_ti^2+ (1-x)*sig_C^2 ; % covariance using Equation 13 of manuscript
corttn(i) = co_rn/(sig_rt*sig_n); % CC using Equation 12 of manuscript
Pttn(i) = co_rn/sig_rt^2; % Slope using Equation 16 of manuscript

end

%%
% figures for report

% plots CC of (Trt, Tn) pair for different values of x. Dotted line
% represents the experimental value of CC of (Trt, Tn) for the particular
% growth media.

figure
hold on
s=plot(x_all, corttn,'-bs','MarkerSize',4, 'LineWidth',2); 
xlabel('x')
ylabel('R(Trt, Tn)')
box on
hold on
pl= plot(x_all, ones(1,length(x_all))*r_exp, 'b--', 'LineWidth',2);
set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%%

% plots slope of best linear fit of Tn vs Trt plot for different values of 
% x. Dotted line represents the experimental value of same for the 
% particular growth media.

figure
hold on
s=plot(x_all, Pttn,'-bs','MarkerSize',4, 'LineWidth',2); 
xlabel('x')
ylabel('Slope(Trt, Tn)')
box on
hold on
pl= plot(x_all, ones(1,length(x_all))*P_exp, 'b--', 'LineWidth',2);
set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%%

% plots CV(Tn-Trt) for different values of x. Dotted line represents the 
% experimental value of the CV for the particular growth media.

figure
hold on
s=plot(x_all, cvtnt,'-bs','MarkerSize',4 ,'LineWidth',2);
xlabel('x')
ylabel('CV(Tn-Trt)')
box on
hold on
pl= plot(x_all, ones(1,length(x_all))*cv_exp, 'b--', 'LineWidth',2);
set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

