%%% THIS CODE IS CALLED BY runNeurophonic.m %%%

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

function OutStruct = MSO_dae(tEnd, stimType, gE, tauE, csynE, gI, tauI, csynI, synFreq, synDelay, inhibDelay, FreezeKLT, rB, varargin)

    % INPUT PARAMETERRS
    optargin = size(varargin,2);
    stdargin = nargin - optargin;
    if stdargin==0
        tEnd = 5;              % Simulation time
        stimType = 'left';     % stimType left, right, or both
        gE = 20;               % gE [mS / cm2]
        tauE = 0.2;            % tauE [ms, alpha function]
        gI = 0;                % gI [mS / cm2]
        tauI = [0.4 2.];       % tauI [ms, double exponential]
        synFreq = [1000 1000 ];% synFreq EPSG frequency (Hz) on each dendrite [2 elements, inhibition uses first]
        synDelay = [0.2 0.2];  % synDelay [2 elements] delay start of EPSG train on left and right dendrite inhibition uses first element 
        inhibDelay = -0.35;    % inhibDelay -- additional amount to delay inhibition [relative to left dendrite]
        FreezeKLT = 0;         % [0 or 1]
        csynE = [2 23-1];      % location (compartment #) of excitatory inputs
        csynI = 12;            % location (compartment #) of inhibitory input
        rB = 11;               % micro m.  must be bigger than radius of soma=10
    end
    
    % Simulation time
    ParamStruct.t0   = 0;
    ParamStruct.tEnd = tEnd;  % ms
    ParamStruct.dt   = 1e-3;  % ms
    
    % Neuron Parameter Values
    ParamStruct.Cm       = 0.9;  % micro F / cm2
    ParamStruct.Ri       = 200.; % axial resistivity [Ohm cm]
    ParamStruct.Gleak    = 0.3;  % mS / cm2
    ParamStruct.Vleak    = -60.; % mV used for dendrites and soma
    ParamStruct.VK       = -106.;% mV
    ParamStruct.Vh       = -43.; % mV
    ParamStruct.Vrest    = -59.6;% mV
    
    % Soma Parameters
    ParamStruct.nS       = 3;     % # soma compartments
    ParamStruct.GKLVAS   = 17.;   % mS / cm2
    ParamStruct.GhS      = 0.86;  % mS / cm2 
    ParamStruct.diamS    = 20E-4; % soma diameter (cm)
    ParamStruct.lS       = 20.E-4;% soma length (cm)
    ParamStruct.dxS      = ParamStruct.lS/ParamStruct.nS;
    ParamStruct.SurfaceS = pi * ParamStruct.diamS * ParamStruct.dxS;
    ParamStruct.CrossS   = pi * (ParamStruct.diamS/2.).^2. ;
    
    % Dendrite Parameters
    ParamStruct.nD       = 10; % # dendrite compartments
    ParamStruct.GKLVAD   = ParamStruct.GKLVAS*.21 ;% [step model]
    ParamStruct.GhD      = ParamStruct.GhS*.21; %  [step model]
    ParamStruct.diamD    = 3.5E-4;   % dendrite diameter (cm)
    ParamStruct.lD       = 150.E-4;  % dendrite length (cm)
    ParamStruct.dxD      = ParamStruct.lD/ParamStruct.nD;
    ParamStruct.SurfaceD = pi * ParamStruct.diamD * ParamStruct.dxD ;
    ParamStruct.CrossD   = pi * (ParamStruct.diamD/2.).^2.; 

    % Put together to make neuron
    ParamStruct.nC       = 2*ParamStruct.nD + ParamStruct.nS;
    ParamStruct.GKLVA    = [ParamStruct.GKLVAD+zeros(1,ParamStruct.nD)    , ParamStruct.GKLVAS+zeros(1,ParamStruct.nS)    , ParamStruct.GKLVAD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.Gh       = [ParamStruct.GhD+zeros(1,ParamStruct.nD)       , ParamStruct.GhS+zeros(1,ParamStruct.nS)       , ParamStruct.GhD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.dx       = [ParamStruct.dxD+zeros(1,ParamStruct.nD)       , ParamStruct.dxS+zeros(1,ParamStruct.nS)       , ParamStruct.dxD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.Cross    = [ParamStruct.CrossD+zeros(1,ParamStruct.nD)    , ParamStruct.CrossS+zeros(1,ParamStruct.nS)    , ParamStruct.CrossD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.Surface  = [ParamStruct.SurfaceD+zeros(1,ParamStruct.nD)  , ParamStruct.SurfaceS+zeros(1,ParamStruct.nS)  , ParamStruct.SurfaceD+zeros(1,ParamStruct.nD) ]';

    % Extracellular parameters
    ParamStruct.dxG      = 0.1; % Distance to ground (cm)
    ParamStruct.RE       = 1.5*ParamStruct.Ri;  % resistance Ohm cm 
    ParamStruct.rB       = rB*1e-4; % Radius of "bounding volume" (convert um to cm)
    ParamStruct.ECross   = pi*ParamStruct.rB^2 - ParamStruct.Cross;
    ParamStruct.ECross0  = pi*ParamStruct.rB^2;
    
    % Spatial coordinates [micro m]
    x= zeros(size(ParamStruct.dx));
    x(1) = ParamStruct.dx(1)*1E4/2-(ParamStruct.lD+ParamStruct.lS/2)*1E4;
    for i=2:size(ParamStruct.dx)
     x(i) = x(i-1)+.5*(ParamStruct.dx(i-1)+ParamStruct.dx(i))*1E4;
    end
    
    %Excitatory Synapse Parameters
    ParamStruct.gsynE = gE;
    ParamStruct.tausynE = tauE; % ms
    ParamStruct.VsynE = 0; 
    if strcmp(stimType,'left')
        tsynE = [synDelay(1):1000/synFreq(1):ParamStruct.tEnd]'; % time of synaptic events
        nsynE = length(tsynE); % Number synapse events in time
        idsynE = repmat(csynE(1),nsynE,1);
        ParamStruct.synE = [idsynE reshape(tsynE,nsynE,1)];
    elseif strcmp(stimType,'right')
        tsynE = [synDelay(2):1000/synFreq(2):ParamStruct.tEnd]'; % time of synaptic events
        nsynE = length(tsynE); % Number synapse events in time
        idsynE = repmat(csynE(2),nsynE,1);
        ParamStruct.synE = [idsynE reshape(tsynE,nsynE,1)];
    elseif strcmp(stimType,'both')
        % First synapse
            tsynE = [synDelay(1):1000/synFreq(1):ParamStruct.tEnd]'; % time of synaptic events
            nsynE = length(tsynE); % Number synapse events in time
            idsynE = repmat(csynE(1),nsynE,1);
            ParamStruct.synE = [idsynE reshape(tsynE,nsynE,1)];
        % Second synapse
            tsynE = [synDelay(2):1000/(synFreq(2)):ParamStruct.tEnd]'; % time of synaptic events
            nsynE = length(tsynE); % Number synapse events in time
            idsynE = repmat(csynE(2),nsynE,1); 
            ParamStruct.synE = [ParamStruct.synE ; [idsynE reshape(tsynE,nsynE,1)]];
    end
    
    %Inhibitory Synapse Parameters
    ParamStruct.gsynI = gI;
    ParamStruct.tausynI = tauI; % ms
    ParamStruct.VsynI = -90; 
    ParamStruct.synI = [];
    if gI>0 
        tsynI = [synDelay(1)+inhibDelay:1000/synFreq(2):ParamStruct.tEnd]'; % time of synaptic events
        tsynI = tsynI(tsynI>0);
        nsynI = length(tsynI); % Number synapse events in time
        idsynI = repmat(csynI,nsynI,1);
        ParamStruct.synI = [idsynI reshape(tsynI,nsynI,1)];
    end

    % IKLT 
    ParamStruct.FreezeKLT = FreezeKLT;
    
    %%%%%%%%%% RUN SOLVER %%%%%%%%%%
    [t,y] = DAESolveIt(ParamStruct);

    nt = length(t);
    dt = ParamStruct.dt;

    % Intracellular voltage
    Vi = y(1:ParamStruct.nC,:)';
    
    % Extracellular Voltage
    Ve = y(ParamStruct.nC+1:2*ParamStruct.nC,:)';
    
    % Membrane potenial
    Vm = Vi-Ve;
    
    % Membrane Current densities
    gate = gate_func;
    m = y(2*ParamStruct.nC+1:3*ParamStruct.nC,:)';
    h = y(3*ParamStruct.nC+1:4*ParamStruct.nC,:)';

    % Currents 
    if ParamStruct.FreezeKLT==0
        IKLT = repmat(ParamStruct.GKLVA',nt,1) .* m.^4. .* h .* (Vm - ParamStruct.VK); % dynamic
    else
        IKLT = repmat(ParamStruct.GKLVA',nt,1) .* repmat(m(1,:),nt,1).^4. .* repmat(h(1,:),nt,1) .* (Vm - ParamStruct.VK);  % Freeze
    end
    Ih   = repmat(ParamStruct.Gh',nt,1) .* (Vm - ParamStruct.Vh);
    Ilk  = ParamStruct.Gleak .* (Vm - ParamStruct.Vleak);
     
    % Synapses
    IsynE = zeros(size(IKLT));
    [nsynE,~] = size(ParamStruct.synE);
    IsynI = zeros(size(IKLT));
    [nsynI,~] = size(ParamStruct.synI);
    for i=1:nt; 
       for j=1:nsynE
           jE = ParamStruct.synE(j,1);
       IsynE(i,jE) = IsynE(i,jE) + ParamStruct.gsynE*(Vm(i,jE)-ParamStruct.VsynE).*alpha(t(i),ParamStruct.tausynE,ParamStruct.synE(j,2));
       end
       for j=1:nsynI
           jI = ParamStruct.synI(j,1);
       IsynI(i,jI) = IsynI(i,jI) + ParamStruct.gsynI*(Vm(i,jI)-ParamStruct.VsynI).*dblExp(t(i),ParamStruct.tausynI,ParamStruct.synI(j,2));
       end
    end
    Isyn = IsynE + IsynI;
     
    % Cable current  [1000/R * V * CrossArea/dx = 1000/(Ohm cm) * mV *cm = 1000mV/Ohm = muA] 
    % Divide by area (cm2) below to get density muA/cm2
    nC = ParamStruct.nC;
    IlongI     = zeros(size(Vi));
    IlongI(:,1)  =  (1000./ParamStruct.Ri) *  ((Vi(:,2)-Vi(:,1))) * ParamStruct.CrossD / ParamStruct.dxD;  % Sealed end boundary
    IlongI(:,nC) =  (1000./ParamStruct.Ri) *  ((Vi(:,ParamStruct.nC-1)-Vi(:,ParamStruct.nC))) * ParamStruct.CrossD / ParamStruct.dxD ; % Sealed end boundary
    for j=2:nC-1
      IlongI(:,j) = (2.*1000./ParamStruct.Ri) * ( (Vi(:,j-1)-Vi(:,j))/(ParamStruct.dx(j-1)/ParamStruct.Cross(j-1)+ ParamStruct.dx(j)  /ParamStruct.Cross(j)) ...
                  + (Vi(:,j+1)-Vi(:,j))/(ParamStruct.dx(j)  /ParamStruct.Cross(j)  + ParamStruct.dx(j+1)/ParamStruct.Cross(j+1)) ) ; 
    end
    IlongI = IlongI./repmat(ParamStruct.Surface',nt,1); % Current Density muA/ cm2

    % Capacitance current
    Ic = -(IKLT + Ih + Ilk + Isyn) + IlongI;

    % Total membrane current density (should integrate to 0 wrt x when scaled by area)
    Im = IlongI;
        
    % Output variables
    OutStruct.t = t;
    OutStruct.x = x;
    OutStruct.Vm = Vm;
    OutStruct.Vi = Vi;
    OutStruct.Ve = Ve;
    OutStruct.Im = Im;
    OutStruct.IKLT = IKLT;
    OutStruct.Ih = Ih;
    OutStruct.Ilk = Ilk;
    OutStruct.IsynE = IsynE;
    OutStruct.IsynI = IsynI;
    OutStruct.Isyn = Isyn;
    OutStruct.Ic = Ic;
    OutStruct.ParamStruct = ParamStruct;

function out  = gate_func()
	out.minf   = @(V) 1.0 ./(1.+exp((V-(-57.34))/(-11.7)));
    out.taum   = @(V) ((21.5 ./ (6.*exp((V+60.)/7.) + 24.*exp(-(V+60.)./50.6))) + 0.35); %((21.5 ./ (6.*exp((V+60.)/7.) + 24.*exp(-(V+60.)./50.6))) + 0.35);
    out.hinf   = @(V) (1.-.27) ./ (1 + exp((V-(-67))/6.16)) + .27;
    out.tauh   = @(V) (170./(4.9*exp((V+60.)/10.) + exp(-(V+70.)/8.)) + 10.7);
    out.ninf   = @(V) (1+exp(-(V+15)/5)).^-0.5;
    out.taun   = @(V) 0.17*100 * (11*exp((V+60)/24) + 21*exp(-(V+60)/23)  ).^(-1) + 0.7;
    out.pinf   = @(V) 1./ (1+exp(-(V+23)/6));
    out.taup   = @(V) 0.17*100 * ( 4*exp((V+60)/32) + 5*exp(-(V+60)/22)  ).^(-1) + 5; 
    
function out = alpha(t,tau,t0)
    out = (t>t0)*(t-t0)/tau * exp(1-(t-t0)/tau);

function out = dblExp(t,tau,t0)  
   a = tau(1)*tau(2)*log(tau(1)/tau(2))/(tau(1)-tau(2));
   m = -exp(-(a)/tau(1))+exp(-(a)/tau(2));
   out = (1/m)*(t>t0)*(-exp(-(t-t0)/tau(1))+exp(-(t-t0)/tau(2)));
   
function [t,y] = DAESolveIt(ParamStruct)

    % Gating functions
    gate = gate_func;

    % Spatial coordinates
    x= zeros(size(ParamStruct.dx));
    x(1) = ParamStruct.dx(1)*1E4/2-(ParamStruct.lD+ParamStruct.lS/2)*1E4;
    for i=2:size(ParamStruct.dx)
      x(i) = x(i-1)+.5*(ParamStruct.dx(i-1)+ParamStruct.dx(i))*1E4;
    end

    % Differential variables : VI, VE, m, h
    nVar = 4;
    id = ones(nVar*ParamStruct.nC,1);

    % IDA solver options
    options = IDASetOptions('RelTol',1.e-6,...
                            'AbsTol',1.e-6,...
                            'VariableTypes',id,...
                            ...%  'MaxStep', 1e-8,...
                            ...%'MaxNumSteps', 100,...
                            'UserData', ParamStruct); % Structure of parameter values

    % Initialize variables
    Vrest0 = ParamStruct.Vrest*ones(2*ParamStruct.nD+ParamStruct.nS,1);%-58; % Resting potential (mV)
    Vi0 = Vrest0;
    Ve0 = zeros(size(Vrest0));
    m0 = gate.minf(Vrest0);
    h0 = gate.hinf(Vrest0);
    y0 = [Vi0 ; Ve0 ; m0 ; h0];
    yp0 = zeros(size(y0));

    % Initialize solver
    IDAInit(@res_dae,ParamStruct.t0,y0,yp0,options);
    [status, y0_mod, yp0_mod] = IDACalcIC(ParamStruct.tEnd, 'FindAlgebraic');

    % % Change Parameter Values  
    % IDASet('UserData',ParamStruct);

    % Run IDA solver
    tic;
    [status,t, y] = IDASolve([ParamStruct.dt:ParamStruct.dt:ParamStruct.tEnd],'Normal');
    display(['Compute time = ', num2str(toc),'sec'])

    t = [ParamStruct.t0 t];
    y = [y0 y];

    % Deallocate IDAS memory
    IDAFree


function [res, flag, new_data] = res_dae(t,y,yp,ParamStruct) 
    % Derivatives: Vmp: Cm*Vmp - (Ileak + IKLVA + Ih + Isyn)) - IlongI =0
    % Derivatives: Vmp: Cm*Vmp - (Ileak + IKLVA + Ih + Isyn)) + IlongE =0
    %              xp - (xinf-x)/taux

    % Variables [Dendrite Soma Dendrite]
    Vi = y(1:ParamStruct.nC);       
    Ve = y(ParamStruct.nC+1:2*ParamStruct.nC);       
    m =  y(2*ParamStruct.nC+1:3*ParamStruct.nC);
    h =  y(3*ParamStruct.nC+1:4*ParamStruct.nC);

    % Derivatives
    Vip =  yp(1:ParamStruct.nC);        
    Vep =  yp(ParamStruct.nC+1:2*ParamStruct.nC);        
    mp  =  yp(2*ParamStruct.nC+1:3*ParamStruct.nC);
    hp  =  yp(3*ParamStruct.nC+1:4*ParamStruct.nC);

    % Membrane potential
    Vm  = Vi  - Ve;
    Vmp = Vip - Vep;

    % Current density (mA / cm^2)
    if ParamStruct.FreezeKLT==0
        IKLT = ParamStruct.GKLVA .* m.^4. .* h .* (Vm - ParamStruct.VK); % dynamic
    else
        minfvrest = 1.0 ./(1.+exp((ParamStruct.Vrest-(-57.34))/(-11.7)));
        hinfvrest = (1.-.27) ./ (1 + exp((ParamStruct.Vrest-(-67))/6.16)) + .27;
        IKLT = ParamStruct.GKLVA .* minfvrest.^4. .* hinfvrest .* (Vm - ParamStruct.VK);  % Freeze
    end
    Ih   = ParamStruct.Gh    .* (Vm - ParamStruct.Vh);
    Ilk  = ParamStruct.Gleak .* (Vm - ParamStruct.Vleak);

    nC = ParamStruct.nC;
    IlongI     = zeros(size(Ilk)); %mu A
    IlongI(1)  =  (1000./ParamStruct.Ri) *  ((Vi(2)-Vi(1))) * ParamStruct.CrossD / ParamStruct.dxD;  % Sealed end boundary
    IlongI(nC) =  (1000./ParamStruct.Ri) *  ((Vi(ParamStruct.nC-1)-Vi(ParamStruct.nC))) * ParamStruct.CrossD / ParamStruct.dxD ; % Sealed end boundary

    for j=2:nC-1
      IlongI(j) = (2.*1000./ParamStruct.Ri) * ( (Vi(j-1)-Vi(j))/(ParamStruct.dx(j-1)/ParamStruct.Cross(j-1)+ ParamStruct.dx(j)  /ParamStruct.Cross(j)) ...
                      + (Vi(j+1)-Vi(j))/(ParamStruct.dx(j)  /ParamStruct.Cross(j)  + ParamStruct.dx(j+1)/ParamStruct.Cross(j+1)) ) ; 
    end

    % mu A
    IlongE      = zeros(size(Ilk));
    IlongE(1)   =  (2*1000.) * ( -Ve(1)/(2*ParamStruct.dxG/(ParamStruct.ECross0) + (ParamStruct.dxD)/ParamStruct.ECross(1))*(1/ParamStruct.RE) ...
                                + (Ve(2)-Ve(1))/(ParamStruct.dxD/ParamStruct.ECross(1) + ParamStruct.dxD/ParamStruct.ECross(2))*(1/ParamStruct.RE) ) ; 
    IlongE(ParamStruct.nC) =  (2*1000.) * ( -Ve(ParamStruct.nC)/(2*ParamStruct.dxG/(ParamStruct.ECross0) + (ParamStruct.dxD)/ParamStruct.ECross(ParamStruct.nC))*(1/ParamStruct.RE) ...
                                + (Ve(ParamStruct.nC-1)-Ve(ParamStruct.nC))/(ParamStruct.dxD/ParamStruct.ECross(ParamStruct.nC) + ParamStruct.dxD/ParamStruct.ECross(ParamStruct.nC-1))*(1/ParamStruct.RE) ) ;
    for j=2:ParamStruct.nC-1
        IlongE(j) = (2.*1000./ParamStruct.RE) * ( (Ve(j-1)-Ve(j))/(ParamStruct.dx(j-1)/ParamStruct.ECross(j-1)+ ParamStruct.dx(j)/ParamStruct.ECross(j)) ...
                                    + (Ve(j+1)-Ve(j))/(ParamStruct.dx(j)/ParamStruct.ECross(j)  + ParamStruct.dx(j+1)/ParamStruct.ECross(j+1)) ) ; 
    end
    
    % Synapses
    IsynE = zeros(size(IKLT));
    [nsynE,~] = size(ParamStruct.synE);    
    IsynI = zeros(size(IKLT));
    [nsynI,~] = size(ParamStruct.synI);    
    activeSynE = find(ParamStruct.synE(:,2)<t  & (t-ParamStruct.synE(:,2))<10*ParamStruct.tausynE  );
    if nsynI>0
        activeSynI = find(ParamStruct.synI(:,2)<t  & (t-ParamStruct.synI(:,2))<10*ParamStruct.tausynI(2)  );
    else activeSynI = [];
    end
    for j=1:length(activeSynE)
      i =ParamStruct.synE(activeSynE(j),1);
      IsynE(i) = IsynE(i) + ParamStruct.gsynE*(Vm(i)-ParamStruct.VsynE).*alpha(t,ParamStruct.tausynE,ParamStruct.synE(activeSynE(j),2));
    end
    for j=1:length(activeSynI)
      i = ParamStruct.synI(activeSynI(j),1); 
      IsynI(i) = IsynI(i) + ParamStruct.gsynI*(Vm(i)-ParamStruct.VsynI).*dblExp(t,ParamStruct.tausynI,ParamStruct.synI(activeSynI(j),2));
    end
    Isyn = IsynE + IsynI;

    % Membrane current density
    Im = ParamStruct.Cm * Vmp  + IKLT + Ih + Ilk + Isyn ;

    % residuals for DAE solver
    res = [ParamStruct.Surface.*Im - IlongI 
           ParamStruct.Surface.*Im + IlongE
           mp   - ( (1.0 ./(1.+exp((Vm-(-57.34))/(-11.7))))  - m)   ./  (((21.5 ./ (6.*exp((Vm+60.)/7.) + 24.*exp(-(Vm+60.)./50.6))) + 0.35))
           hp   - ( ((1.-.27) ./ (1 + exp((Vm-(-67))/6.16)) + .27)   - h)   ./  (170./(4.9*exp((Vm+60.)/10.) + exp(-(Vm+70.)/8.)) + 10.7)];

    flag = 0;
    new_data = [];