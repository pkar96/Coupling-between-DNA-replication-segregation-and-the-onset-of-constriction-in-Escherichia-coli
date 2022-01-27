% calls codes which simulate various cell cycle models such as adder per
% origin model, Cooper Helmstetter model, parallel adder model and the
% independent adder model. Change the function called name in this section
% to simulate a different model.

clear

no_runs= 6; % The simulations are run no_runs times

% Output obtained
times=NaN(no_runs,6000); % timing of simulation relative to thymine shift time
cd_pop=NaN(no_runs,6000); % C+D timing of the dividing population at time "times" 
td_pop=NaN(no_runs,6000); % Td timing of the dividing population at time "times" 
lb_pop=NaN(no_runs,6000); % Length at birth of the dividing population at time "times" 
ld_pop=NaN(no_runs,6000); % Length at division of the dividing population at time "times" 
li_pop=NaN(no_runs,6000); % Length at initiation per origin of the dividing population at time "times" 
li_t_pop=NaN(no_runs,6000); % Length at initiation of the dividing population at time "times" 
rate_pop=NaN(no_runs,6000); % growth rates of the dividing population at time "times" 

% runs parallel simulations of the models. The first parameter in the 
% called function chooses the E. coli strain to be simulated. The
% parameters used in the cell cycle models are chosen from experiments
% conducted on that strain. 1= STK37 strain
parfor i=1:no_runs
    [times(i,:), cd_pop(i,:), td_pop(i,:), lb_pop(i,:), ld_pop(i,:), li_pop(i,:), li_t_pop(i,:), rate_pop(i,:)] = indep_adder_change_c(1, i);
end

%%

% stores the output obtained from the run in an excel file
fileName= 'ao_model.xlsx';
xlswrite(fileName,times,'times');
xlswrite(fileName,cd_pop,'C+D');
xlswrite(fileName,td_pop,'Td');
xlswrite(fileName,lb_pop,'Lb');
xlswrite(fileName,ld_pop,'Ld');
xlswrite(fileName,li_pop,'Del_i');
xlswrite(fileName,li_t_pop,'Li');
xlswrite(fileName,rate_pop,'Rate');

%%
% Plot figures similar to that shown in Figure 5

% Extract the data from excel files stored in previous section
[numx,txtx,rawx] = xlsread('ao_model.xlsx','times'); 
[numy,txty,rawy] = xlsread('ao_model.xlsx','Td');

no_runs=20; % same as that in first section
x_tot=-300:25:600; % binned edges
y_tot= zeros(no_runs,length(x_tot)-1); % averaged cell characteristics
for i=1:no_runs
x= numx(i,:); y= numy(i,:);
[bin,da,yfit, P, err]=binning_with_error_1(y,x, x_tot); % da stores the average data
y_tot(i,:) = da;
end

% plots the figures to generate plots similar to Figure 5
figure
fp=plot((x_tot(1:end-1)+x_tot(2:end))/2,mean(y_tot),'-r');
fp.LineWidth=1.5;
ylabel('T_d (min)');
xlabel('Time (min)');
box on
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,642])
xlim([-300,600])
