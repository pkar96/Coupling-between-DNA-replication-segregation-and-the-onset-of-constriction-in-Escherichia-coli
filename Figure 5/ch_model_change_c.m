% The Cooper Helmstetter model is simulated for a population of cells. When
% a cell divides information about the cell cycle such as C+D time,
% generation time, length at birth, length at division, length at
% initiation per origin, length at initiation and the growth rate for the 
% cell cycle is stored. In the simulation at simTime=ch_cd, the C period 
% changes due to thymine shift. This is captured by dcd minutes increase in
% the parameter cd (captures C+D time). For more information about model, 
% see "Modeling of thymine shift experiments" section in the manuscript.

function [times, cd_pop, td_pop, lb_pop, ld_pop, li_pop, li_t_pop, rate_pop] = ch_model_change_c(num_in, iter_num)
iter_num

[num,txt,raw] = xlsread('D:\cell division\exp data\jaan mannik new data\jaan paper figures\Figure 5 new\params_thy.xlsx');
% The excel file above contains parameters used in the simulations. The values are derived from experimental data.

cells=[];
% Input parameters
v_inp=1;
tau= num(num_in,1); % generation time in minutes
gr= log(2)/tau; % growth rate in min^-1
ngen=12; % number of generation
cd= num(num_in,4);  % C+D period mean in minutes
bii=num(num_in,10)*2;   % change in length from initiation to initiation
kp_cells= 5*tau; % time till which both daughter cells are noted in the experiment
ch_cd= 7*tau; % the C-period changes by dcd minutes at this time in the experiment
dcd= 40.1999; % change C-period by this many minutes after ch_cd

%-----Noise parameters
cvl=num(num_in,3)*num(num_in,2);  % standard deviation of growth rate
cvt=num(num_in,12)*bii; % standard deviation of size additive initiation adder noise
sigcd=num(num_in,5)*cd;  % standard deviation of time additive C+D noise

%-----Outputs
cd_pop=NaN(1,6000); % C+D timing of the population
td_pop=NaN(1,6000); % Generation time of the population
lb_pop=NaN(1,6000); % Length at birth of the population
ld_pop=NaN(1,6000); % Length at division of the population
li_pop=NaN(1,6000); % Length at initiation per origin of the population
li_t_pop=NaN(1,6000); % Length at initiation of the population
rate_pop=NaN(1,6000); % rates of the population
times=NaN(1,6000); % timing of the simulation

lmin=gr/(5);

for no_in=1:15
x= Cell(v_inp); % Obtains an object Cell with defined cell characteristics
x.rate= gr;

%-----For keeping track of events
ii_in = bii+ randn()*cvt;
while ii_in<0
    ii_in=bii+ randn()*cvt;
end
x.vNextInit = x.oris*(ii_in);

x.tNextDivs(1)=cd+randn()*sigcd;
while x.tNextDivs(1)<0
    x.tNextDivs(1)=cd+randn()*sigcd;
end

cells=[cells x];
end

%-----Advance in time for ngen generations
gens = ngen*tau; % total time in mins for which the simulation runs
tStep = 0.01; % time increment in units of mins
simTime = 0; % simulation time progress
post_ch= 0; % counter tells if ch_cd time has passed from beginning of simulation
cnt=1; % array index counter to store cell characteristics information after cells have divided

while simTime < gens
    %-----Step
    simTime = simTime + tStep;
    
    % check for time at which thymine shift occurs
    if simTime>ch_cd && post_ch==0
        cd=cd+dcd;
        post_ch= 1;
    end
    
    for ct=1:length(cells)
    x=cells(ct);
    x.t = x.t + tStep;
    x.v = x.v*exp(tStep*x.rate); % grow exponentially
    x.tNextDivs= x.tNextDivs-tStep;
    end
    %-----Perform events
    
    cellsNew=[];
    for ct=1:length(cells)
    x=cells(ct);
    %-----Divide
    %Upon accumulating enough length
    if ~isempty(x.tNextDivs) && x.tNextDivs(1)<=0
        %-----Record event
        x.td_cell = x.t-x.tLastDiv;
        x.cd_cell = x.t - x.tOfInits(1);
        x.lb_cell = x.vb;
        x.ld_cell = x.v;
        x.li_cell = x.vOfInits(1)/x.oOfInits(1)*2;
        x.li_t_cell = x.vOfInits(1);
        x.rate_cell= x.rate;
        
        x.vOfInits(1)=[];
        x.tOfInits(1)=[];
        x.oOfInits(1)=[];
        
        %-----Update cell
        
        x.tLastDiv = x.t;
        x.vd = x.v;
        x.v = x.vd*0.5;
        x.vb = x.v;
        x.rec_data=1;
        x.rate = max(gr + randn()*cvl,lmin);
        x.oris = x.oris/2;
        x.tNextDivs(1)=[];
        %----------Update vNextInit
        x.vNextInit = x.vNextInit/2;
        
        %----------other daughter
        if simTime<=kp_cells
        y=copy(x);
        y.td_cell= NaN;
        y.cd_cell= NaN;
        y.lb_cell = NaN;
        y.ld_cell = NaN;
        y.li_cell = NaN;
        y.li_t_cell = NaN;
        y.rate_cell = NaN;
        y.rec_data=0;
        y.rate= max(gr + randn()*cvl,lmin);
        cellsNew =[cellsNew y];
        end
    end
    
    %-----Initiate
    %Upon accumulating enough length
    if x.v - x.vNextInit >= 0
        %-----Record event
        x.tLastInit = x.t;
        x.vi = x.v;
        x.vOfInits= [x.vOfInits x.v];
        x.tOfInits= [x.tOfInits x.t];
        x.oOfInits= [x.oOfInits 2*x.oris];
        
        %-----Update cell
        x.oris = 2*x.oris; % double number of origins
        ii_in = bii+ randn()*cvt;
        while ii_in<0
            ii_in=bii+ randn()*cvt;
        end
        x.vNextInit = x.oris*(ii_in);
        xcd=cd+randn()*sigcd;
        while xcd<0
            xcd=cd+randn()*sigcd;
        end
        x.tNextDivs= [x.tNextDivs xcd];
    end
    end
    cells = [cells cellsNew];
    
    % record output for dividing cells
    ind_in = find([cells.rec_data]==1);
    cd_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).cd_cell];
    td_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).td_cell];
    lb_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).lb_cell];
    ld_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).ld_cell];
    li_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).li_cell];
    li_t_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).li_t_cell];
    rate_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).rate_cell];
    times(cnt:cnt+length(ind_in)-1)= (simTime-ch_cd)*ones(1,length(ind_in));
    cnt=cnt+length(ind_in);
    for cnt_ind=1:length(ind_in)
        cells(ind_in(cnt_ind)).rec_data=0;
    end
    
end
