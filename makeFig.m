%%% THIS CODE CALLS MSO_dae.m %%%

% Computes Vm of MSO neuron model and Ve in extracellular "virtual cylinder"
% Accompanies the manuscript [submitted to J. Neuroscience]:
% "A model of the medial superior olive explains spatiotemporal features of local field potentials"
% JH Goldwyn, M Mc Laughlin, E Verschooten, PX Joris, J Rinzel
%
% MSO model neuron was introduced in:
% "Control of submillisecond synaptic timing in binaural coincidence detectors by Kv1 channels"
% Paul J Mathews, Pablo E Jercog,	John Rinzel, Luisa L Scott, Nace L Golding
% Nature Neuroscience 13 601-609 (2010)

% Simulation code by Joshua H Goldwyn
% Submitted to ModelDB 1/14/13 by Joshua H Goldwyn [jgoldwyn@nyu.edu]

close all
clear all

% Reproduce figures from manuscript [1=Yes, 0=No]
MakeFig4  = 1; 
MakeFig8  = 1;
MakeFig11 = 1;


%%%     MAKE FIGURE 8 -- MONOLATERAL 1kHZ EXCITATION  %%%
if MakeFig4
    
    figure(4), clf
    set(gcf,'position',[1           1        1436         805])
    
    %%% MONOLATERAL 1kHz EXCITATION %%%
    tEnd = 7.;             % simulation duration [ms]
    stimType = 'left';     % monolateral excitation
    gE = 10;               % excitatory conductance [mS / cm2]
    tauE = 0.2;            % excitatory time constant (alpha function) [ms]
    csynE = [2 22];        % location of excitation (compartment number)
    gI = 0;                % inhibitory conductance [mS / cm2]
    tauI = [.4 2];         % inhibitory time constants (double exponential function) [ms]
    csynI = [12];          % location of inhibition (compartment number)
    synFreq = [1000 1001]; % EPSP frequency (Hz) for each dendrite. inhibition freq is first entry
    synDelay = [.0 .0];    % Delay of EPSP onset times in each dendrite [ms]
    inhibDelay = 0;        % Delay of inhibition relative to excitation in first entry of synDelay
    FreezeKLT = 0;         % Whether to Freeze KLT conductance at rest (0=No)
    rB = 11;               % radius of extracellular virtual cylinder (must be larger than soma radius = 10) [micro m]

    % Run model
    out = MSO_dae(tEnd, stimType, gE, tauE, csynE, gI, tauI, csynI, synFreq, synDelay, inhibDelay, FreezeKLT, rB);
    [nt,nx] = size(out.Ve);
    Ve = out.Ve; % Extracellular potential [mV]
    Vm = out.Vm; % Membrane potential [mV]
    Isyn = 1e3*(out.Isyn.*repmat(out.ParamStruct.Surface',nt,1)); % synaptic current [nA]
    Im =1e3* (out.Im.*repmat(out.ParamStruct.Surface',nt,1));  % net membrane current [nA]
    x = out.x; % spatial location of compartments [micro m, 0 is soma center]
    t = out.t; % time (ms)
    ParamStruct = out.ParamStruct; % Input parameters
    dx = ParamStruct.dx*1e4; % distance between compartments [micro m]
    dt = t(2)-t(1); % time step [ms]
    dxG = out.ParamStruct.dxG*1e4; % Distance to ground [microm m]

    % Time and space grid for surface plots of Vm
    [tSurf, xSurf] = meshgrid([t t(end)+dt], [-160 ;-160+cumsum(dx)]); tSurf = tSurf'; xSurf = xSurf';

    % Time and space grid for surface plots of extracellular domain [extended to 0mV at ground via linear decay]
    [tSurfE, xSurfE] = meshgrid([t t(end)+dt], [[-160-dxG:10:-170]' ; -160 ;-160+cumsum(dx) ; [170:10:160+dxG]']); tSurfE = tSurfE'; xSurfE = xSurfE';

    % Extracellular voltage extended to 0mV at ground via linear decay
    for i=1:nt
        VeE(i,:) = [linspace(0,Ve(i,1),100) Ve(i,:) linspace(Ve(i,end),0,100)] ;
    end

    m = mean(VeE(4001:end,:),1); % Mean of "ongoing" response at each spatial location
    VeE0 = VeE -repmat(m,nt,1); % "DC removed" responses

    FontSize = 24;
    
    % Panel A: Vm
    subplot(2,3,1)
    surf(tSurf,xSurf,zeros(size(tSurf)),Vm); shading flat; view([0, 90])
    axis([4 7 min(min(xSurf)) max(max(xSurf))])
    caxis([-80,-40]);
    set(gca,'XTick',0:tEnd,'YTick',[-160 0 160],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-80 -40],'XTick',[-80:20:-4],'FontSize',FontSize);
    xlabel('Time (ms)', 'FontSize',FontSize)
    ylabel('Distance (\mum)', 'FontSize',FontSize)
    ti = title('Vm (mV)','FontSize',FontSize);
    set(gca,'position',[.13 .6 .2 .2])
    set(ti,'position',get(ti,'position')+[0 160 0])



    % Panel B: Ve 
    subplot(2,3,2)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE); shading flat; view([0, 90])
    caxis(.325*[-1 1])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.325 .325],'XTick',[-.3 0 .3],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve (mV)','FontSize',FontSize);
    set(gca,'position',[.43 .6 .2 .2])
    set(ti,'position',get(ti,'position')+[0 1200 0])


    % Panel C -- Ve [DC removed]
    subplot(2,3,3)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE0); shading flat; view([0, 90])
    caxis(.2*[-1 1])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.2 .2],'XTick',[-.2 0 .2],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve [No DC] (mV)','FontSize',FontSize);
    set(gca,'position',[.76 .6 .2 .2])
    set(ti,'position',get(ti,'position')+[0 1200 0])



    % Panel D : Membrane current 
    subplot(2,3,4)
    surf(tSurf,xSurf,zeros(size(tSurf)),Im); shading flat; view([0, 90])
    axis([4 7 min(min(xSurf)) max(max(xSurf))])
    caxis([-.75 .75]);
    set(gca,'FontSize',FontSize,'Xtick',0:tEnd)
    cb = colorbar('location','northoutside');
    ti = title('I memb (nA)','FontSize',FontSize);
    set(gca,'YTick',[-160 0 160])
    set(cb,'xlim',[ -.75 .15],'XTick',[ -.6:.3:0],'YTick',[-160 0 160],'FontSize',FontSize);
    set(ti,'position',get(ti,'position')+[0 160 0])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ylabel('Distance (\mum)', 'FontSize',FontSize)
    set(gca,'position',[.13 .1 .2 .2])


    % Panel E: Membrane current without synapse current
    subplot(2,3,5)
    surf(tSurf,xSurf,zeros(size(tSurf)),Im-Isyn); shading flat; view([0, 90])
    axis([4 7 min(min(xSurf)) max(max(xSurf))])
    caxis([-.121,.121]);
    set(gca,'FontSize',FontSize,'Xtick',4:tEnd,'YTick',[-160 0 160])
    xlabel('Time (ms)', 'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    ti = title('I memb w/out  I syn (nA)','FontSize',FontSize);
    set(cb,'xlim',[ -.01 .121],'XTick',[ 0 .05 .1],'FontSize',FontSize);
    set(ti,'position',get(ti,'position')+[0 160 0])
    set(gca,'position',[.43 .1 .2 .2])



    % Panel F: Membrane current by region
    % compartments of cell regions
    iDend1 = 1:10; iSoma = 11:13; iDend2 = 14:23;

    subplot(2,3,6), hold all
    plot(t,sum(Isyn'),'k','linewidth',2)
    plot(t,[sum(Im(:,iDend1)') ;sum(Im(:,iSoma)') ; sum(Im(:,iDend2)')],'linewidth',2)
    leg = legend({'Synapse','Near Dend.','Soma','Far Dend.'},'location','northoutside','fontsize',16);
    set(gca,'XTick',0:tEnd,'YTick',-.8:.4:.4,'FontSize',FontSize)
    xlabel('Time (ms)','FontSize',FontSize)
    ylabel('Current (nA)','FontSize',FontSize)
    axis([4 7 -.85 .4])
    legend('boxoff')
    set(gca,'position',[.76 .1 .2 .2])


end



%%%     MAKE FIGURE 8 -- BILATERAL 1kHz EXCITATION  %%%
if MakeFig8
    
    figure(8)
    set(gcf,'position',[1           1        1436         805])
    
    %%% Out of phase (.5ms delay) %%%
    tEnd = 7.;             % simulation duration [ms]
    stimType = 'both';     % bilateral excitation
    gE = 10;               % excitatory conductance [mS / cm2]
    tauE = 0.2;            % excitatory time constant (alpha function) [ms]
    csynE = [2 22];        % location of excitation (compartment number)
    gI = 0;                % inhibitory conductance [mS / cm2]
    tauI = [.4 2];         % inhibitory time constants (double exponential function) [ms]
    csynI = [12];          % location of inhibition (compartment number)
    synFreq = [1000 1000]; % EPSP frequency (Hz) for each dendrite. inhibition freq is first entry
    synDelay = [0  .5];    % Delay of EPSP onset times in each dendrite [ms]
    inhibDelay = 0;        % Delay of inhibition relative to excitation in first entry of synDelay
    FreezeKLT = 0;         % Whether to Freeze KLT conductance at rest (0=No)
    rB = 11;               % radius of extracellular virtual cylinder (must be larger than soma radius = 10) [micro m]

    % Run model
    out = MSO_dae(tEnd, stimType, gE, tauE, csynE, gI, tauI, csynI, synFreq, synDelay, inhibDelay, FreezeKLT, rB);
    [nt,nx] = size(out.Ve);
    Ve = out.Ve; % Extracellular potential [mV]
    Vm = out.Vm; % Membrane potential [mV]
    Isyn = 1e3*(out.Isyn.*repmat(out.ParamStruct.Surface',nt,1)); % synaptic current [nA]
    Im =1e3* (out.Im.*repmat(out.ParamStruct.Surface',nt,1));  % net membrane current [nA]
    x = out.x; % spatial location of compartments [micro m, 0 is soma center]
    t = out.t; % time (ms)
    ParamStruct = out.ParamStruct; % Input parameters
    dx = ParamStruct.dx*1e4; % distance between compartments [micro m]
    dt = t(2)-t(1); % time step [ms]
    dxG = out.ParamStruct.dxG*1e4; % Distance to ground [microm m]

    % Time and space grid for surface plots of Vm
    [tSurf, xSurf] = meshgrid([t t(end)+dt], [-160 ;-160+cumsum(dx)]); tSurf = tSurf'; xSurf = xSurf';

    % Time and space grid for surface plots of extracellular domain [extended to 0mV at ground via linear decay]
    [tSurfE, xSurfE] = meshgrid([t t(end)+dt], [[-160-dxG:10:-170]' ; -160 ;-160+cumsum(dx) ; [170:10:160+dxG]']); tSurfE = tSurfE'; xSurfE = xSurfE';

    % Extracellular voltage extended to 0mV at ground via linear decay
    for i=1:nt
        VeE(i,:) = [linspace(0,Ve(i,1),100) Ve(i,:) linspace(Ve(i,end),0,100)] ;
    end

    m = mean(VeE(4001:end,:),1); % Mean of "ongoing" response at each spatial location
    VeE0 = VeE -repmat(m,nt,1); % "DC removed" responses

    FontSize = 24;
    
    % Panel A: Vm
    subplot(2,3,1)
    surf(tSurf,xSurf,zeros(size(tSurf)),Vm); shading flat; view([0, 90])
    axis([4 7 min(min(xSurf)) max(max(xSurf))])
    caxis([-80,-40]);
    set(gca,'XTick',0:tEnd,'YTick',[-160 0 160],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-80 -40],'XTick',[-80:20:-4],'FontSize',FontSize);
    xlabel('Time (ms)', 'FontSize',FontSize)
    ylabel('Distance (\mum)', 'FontSize',FontSize)
    ti = title('Vm (mV)','FontSize',FontSize);
    set(gca,'position',[.13 .6 .2 .2])
    set(ti,'position',get(ti,'position')+[0 160 0])



    % Panel B: Ve 
    subplot(2,3,2)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE); shading flat; view([0, 90])
    caxis(.35*[-1 1])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.35 .35],'XTick',[-.3 0 .3],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve (mV)','FontSize',FontSize);
    set(gca,'position',[.43 .6 .2 .2])
    set(ti,'position',get(ti,'position')+[0 1200 0])


    % Panel C -- Ve [DC removed]
    subplot(2,3,3)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE0); shading flat; view([0, 90])
    caxis([-.21 .21])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.21 .21],'XTick',[-.2 0 .2],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve [No DC] (mV)','FontSize',FontSize);
    set(gca,'position',[.76 .6 .2 .2])
    set(ti,'position',get(ti,'position')+[0 1200 0])



    %%% In phase (0ms delay) %%%
    tEnd = 7.;             % simulation duration [ms]
    stimType = 'both';     % bilateral excitation
    gE = 10;               % excitatory conductance [mS / cm2]
    tauE = 0.2;            % excitatory time constant (alpha function) [ms]
    csynE = [2 22];        % location of excitation (compartment number)
    gI = 0;                % inhibitory conductance [mS / cm2]
    tauI = [.4 2];         % inhibitory time constants (double exponential function) [ms]
    csynI = [12];          % location of inhibition (compartment number)
    synFreq = [1000 1000]; % EPSP frequency (Hz) for each dendrite. inhibition freq is first entry
    synDelay = [0 .0];    % Delay of EPSP onset times in each dendrite [ms]
    inhibDelay = 0;        % Delay of inhibition relative to excitation in first entry of synDelay
    FreezeKLT = 0;         % Whether to Freeze KLT conductance at rest (0=No)
    rB = 11;               % radius of extracellular virtual cylinder (must be larger than soma radius = 10) [micro m]

    % Run model
    out = MSO_dae(tEnd, stimType, gE, tauE, csynE, gI, tauI, csynI, synFreq, synDelay, inhibDelay, FreezeKLT, rB);
    [nt,nx] = size(out.Ve);
    Ve = out.Ve; % Extracellular potential [mV]
    Vm = out.Vm; % Membrane potential [mV]
    Isyn = 1e3*(out.Isyn.*repmat(out.ParamStruct.Surface',nt,1)); % synaptic current [nA]
    Im =1e3* (out.Im.*repmat(out.ParamStruct.Surface',nt,1));  % net membrane current [nA]
    x = out.x; % spatial location of compartments [micro m, 0 is soma center]
    t = out.t; % time (ms)
    ParamStruct = out.ParamStruct; % Input parameters
    dx = ParamStruct.dx*1e4; % distance between compartments [micro m]
    dt = t(2)-t(1); % time step [ms]
    dxG = out.ParamStruct.dxG*1e4; % Distance to ground [microm m]

    % Time and space grid for surface plots of Vm
    [tSurf, xSurf] = meshgrid([t t(end)+dt], [-160 ;-160+cumsum(dx)]); tSurf = tSurf'; xSurf = xSurf';

    % Time and space grid for surface plots of extracellular domain [extended to 0mV at ground via linear decay]
    [tSurfE, xSurfE] = meshgrid([t t(end)+dt], [[-160-dxG:10:-170]' ; -160 ;-160+cumsum(dx) ; [170:10:160+dxG]']); tSurfE = tSurfE'; xSurfE = xSurfE';

    % Extracellular voltage extended to 0mV at ground via linear decay
    for i=1:nt
        VeE(i,:) = [linspace(0,Ve(i,1),100) Ve(i,:) linspace(Ve(i,end),0,100)] ;
    end

    m = mean(VeE(4001:end,:),1); % Mean of "ongoing" response at each spatial location
    VeE0 = VeE -repmat(m,nt,1); % "DC removed" responses

    %Panel D: Vm
    subplot(2,3,4)
    surf(tSurf,xSurf,zeros(size(tSurf)),Vm); shading flat; view([0, 90])
    axis([4 7 min(min(xSurf)) max(max(xSurf))])
    caxis([-80,-40]);
    set(gca,'XTick',0:tEnd,'YTick',[-160 0 160],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-80 -40],'XTick',[-80:20:-4],'FontSize',FontSize);
    xlabel('Time (ms)', 'FontSize',FontSize)
    ylabel('Distance (\mum)', 'FontSize',FontSize)
    ti = title('Vm (mV)','FontSize',FontSize);
    set(gca,'position',[.13 .1 .2 .2])
    set(ti,'position',get(ti,'position')+[0 160 0])



    % Panel E: Ve 
    subplot(2,3,5)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE); shading flat; view([0, 90])
    caxis([-.5 .5])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.5 .5],'XTick',[-.5 0 .5],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve (mV)','FontSize',FontSize);
    set(gca,'position',[.43 .1 .2 .2])
    set(ti,'position',get(ti,'position')+[0 1200 0])


    % Panel F -- Ve [DC removed]
    subplot(2,3,6)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE0); shading flat; view([0, 90])
    caxis([-.21 .21])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.21 .21],'XTick',[-.2 0 .2],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve [No DC] (mV)','FontSize',FontSize);
    set(gca,'position',[.76 .1 .2 .2])
    set(ti,'position',get(ti,'position')+[0 1200 0])


end


%%%     MAKE FIGURE 11 -- MONOLATERAL 1kHz EXCITATION + SOMA INHIBIITON %%%
if MakeFig11
    
    figure(11), clf
    set(gcf,'position',[1           300        1436         400])
    
    %%% MONOLATERAL 1kHz EXCITATION with soma inhibition %%%
    tEnd = 7.;             % simulation duration [ms]
    stimType = 'left';     % monolateral excitation
    gE = 10;               % excitatory conductance [mS / cm2]
    tauE = 0.2;            % excitatory time constant (alpha function) [ms]
    csynE = [2 22];        % location of excitation (compartment number)
    gI = 4;                % inhibitory conductance [mS / cm2]
    tauI = [.4 2];         % inhibitory time constants (double exponential function) [ms]
    csynI = [12];          % location of inhibition (compartment number)
    synFreq = [1000 1001]; % EPSP frequency (Hz) for each dendrite. inhibition freq is first entry
    synDelay = [.0 .0];    % Delay of EPSP onset times in each dendrite [ms]
    inhibDelay = -.35;        % Delay of inhibition relative to excitation in first entry of synDelay
    FreezeKLT = 0;         % Whether to Freeze KLT conductance at rest (0=No)
    rB = 11;               % radius of extracellular virtual cylinder (must be larger than soma radius = 10) [micro m]

    % Run model
    out = MSO_dae(tEnd, stimType, gE, tauE, csynE, gI, tauI, csynI, synFreq, synDelay, inhibDelay, FreezeKLT, rB);
    [nt,nx] = size(out.Ve);
    Ve = out.Ve; % Extracellular potential [mV]
    Vm = out.Vm; % Membrane potential [mV]
    Isyn = 1e3*(out.Isyn.*repmat(out.ParamStruct.Surface',nt,1)); % synaptic current [nA]
    Im =1e3* (out.Im.*repmat(out.ParamStruct.Surface',nt,1));  % net membrane current [nA]
    x = out.x; % spatial location of compartments [micro m, 0 is soma center]
    t = out.t; % time (ms)
    ParamStruct = out.ParamStruct; % Input parameters
    dx = ParamStruct.dx*1e4; % distance between compartments [micro m]
    dt = t(2)-t(1); % time step [ms]
    dxG = out.ParamStruct.dxG*1e4; % Distance to ground [microm m]

    % Time and space grid for surface plots of Vm
    [tSurf, xSurf] = meshgrid([t t(end)+dt], [-160 ;-160+cumsum(dx)]); tSurf = tSurf'; xSurf = xSurf';

    % Time and space grid for surface plots of extracellular domain [extended to 0mV at ground via linear decay]
    [tSurfE, xSurfE] = meshgrid([t t(end)+dt], [[-160-dxG:10:-170]' ; -160 ;-160+cumsum(dx) ; [170:10:160+dxG]']); tSurfE = tSurfE'; xSurfE = xSurfE';

    % Extracellular voltage extended to 0mV at ground via linear decay
    for i=1:nt
        VeE(i,:) = [linspace(0,Ve(i,1),100) Ve(i,:) linspace(Ve(i,end),0,100)] ;
    end
    VeE_withI = VeE;
    
    m = mean(VeE(4001:end,:),1); % Mean of "ongoing" response at each spatial location
    VeE0 = VeE -repmat(m,nt,1); % "DC removed" responses

    FontSize = 24;
    
    % Panel A: Vm
    subplot(1,4,1)
    surf(tSurf,xSurf,zeros(size(tSurf)),Vm); shading flat; view([0, 90])
    axis([4 7 min(min(xSurf)) max(max(xSurf))])
    caxis([-80,-40]);
    set(gca,'XTick',0:tEnd,'YTick',[-160 0 160],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-80 -40],'XTick',[-80:20:-4],'FontSize',FontSize);
    xlabel('Time (ms)', 'FontSize',FontSize)
    ylabel('Distance (\mum)', 'FontSize',FontSize)
    ti = title('Vm (mV)','FontSize',FontSize);
    set(gca,'position',[.08 .2 .16 .5])
    set(ti,'position',get(ti,'position')+[0 100 0])
    
    % Panel B: Ve 
    subplot(1,4,2)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE); shading flat; view([0, 90])
    caxis([-.61 .61])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.61 .61],'XTick',[-.6 0 .6],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve (mV)','FontSize',FontSize);
    set(gca,'position',[.32 .2 .16 .5])
    set(ti,'position',get(ti,'position')+[0 700 0])


    % Panel C -- Ve [DC removed]
    subplot(1,4,3)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE0); shading flat; view([0, 90])
    caxis([-.2 .2])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.2 .2],'XTick',[-.2 0 .2],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve [No DC] (mV)','FontSize',FontSize);
    set(gca,'position',[.56 .2 .16 .5])
    set(ti,'position',get(ti,'position')+[0 700 0])
    
    
    %%% MONOLATERAL 1kHz EXCITATION without soma inhibition %%%
    tEnd = 7.;             % simulation duration [ms]
    stimType = 'left';     % monolateral excitation
    gE = 10;               % excitatory conductance [mS / cm2]
    tauE = 0.2;            % excitatory time constant (alpha function) [ms]
    csynE = [2 22];        % location of excitation (compartment number)
    gI = 0;                % inhibitory conductance [mS / cm2]
    tauI = [.4 2];         % inhibitory time constants (double exponential function) [ms]
    csynI = [12];          % location of inhibition (compartment number)
    synFreq = [1000 1001]; % EPSP frequency (Hz) for each dendrite. inhibition freq is first entry
    synDelay = [.0 .0];    % Delay of EPSP onset times in each dendrite [ms]
    inhibDelay = -.35;        % Delay of inhibition relative to excitation in first entry of synDelay
    FreezeKLT = 0;         % Whether to Freeze KLT conductance at rest (0=No)
    rB = 11;               % radius of extracellular virtual cylinder (must be larger than soma radius = 10) [micro m]

    % Run model
    out = MSO_dae(tEnd, stimType, gE, tauE, csynE, gI, tauI, csynI, synFreq, synDelay, inhibDelay, FreezeKLT, rB);
    [nt,nx] = size(out.Ve);
    Ve = out.Ve; % Extracellular potential [mV]
    Vm = out.Vm; % Membrane potential [mV]
    Isyn = 1e3*(out.Isyn.*repmat(out.ParamStruct.Surface',nt,1)); % synaptic current [nA]
    Im =1e3* (out.Im.*repmat(out.ParamStruct.Surface',nt,1));  % net membrane current [nA]
    x = out.x; % spatial location of compartments [micro m, 0 is soma center]
    t = out.t; % time (ms)
    ParamStruct = out.ParamStruct; % Input parameters
    dx = ParamStruct.dx*1e4; % distance between compartments [micro m]
    dt = t(2)-t(1); % time step [ms]
    dxG = out.ParamStruct.dxG*1e4; % Distance to ground [microm m]

    % Time and space grid for surface plots of Vm
    [tSurf, xSurf] = meshgrid([t t(end)+dt], [-160 ;-160+cumsum(dx)]); tSurf = tSurf'; xSurf = xSurf';

    % Time and space grid for surface plots of extracellular domain [extended to 0mV at ground via linear decay]
    [tSurfE, xSurfE] = meshgrid([t t(end)+dt], [[-160-dxG:10:-170]' ; -160 ;-160+cumsum(dx) ; [170:10:160+dxG]']); tSurfE = tSurfE'; xSurfE = xSurfE';

    % Extracellular voltage extended to 0mV at ground via linear decay
    for i=1:nt
        VeE(i,:) = [linspace(0,Ve(i,1),100) Ve(i,:) linspace(Ve(i,end),0,100)] ;
    end
    

    % Panel D -- [Ve w/inhibition - Ve w/out inhibition]
    subplot(1,4,4)
    surf(tSurfE,xSurfE,zeros(size(tSurfE)),VeE_withI - VeE); shading flat; view([0, 90])
    caxis([-.34 .34])
    set(gca,'XTick',0:tEnd,'YTick',[-1000 0 1000],'FontSize',FontSize)
    cb = colorbar('location','northoutside');
    set(cb,'xlim',[-.34 .34],'XTick',[-.3 0 .3],'FontSize',FontSize);
    axis([4 7 min(min(xSurfE)) max(max(xSurfE))])
    xlabel('Time (ms)', 'FontSize',FontSize)
    ti = title('Ve Difference (mV)','FontSize',FontSize);
    set(gca,'position',[.8 .2 .16 .5])
    set(ti,'position',get(ti,'position')+[0 700 0])

    
end


