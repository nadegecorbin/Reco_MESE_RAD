function [F0,Fn,Zn,F] = EPG_MESE(theta,ESP,T1,T2,Nexc,TR, varargin)
%   Single pool EPG (classic version) for TSE sequences
%   [F0,Fn,Zn,F] = EPG_MESE(theta,ESP,T1,T2,Nexc,TR,varargin)
%
%   arguments:
%               theta:      vector of flip angles (rad) including excitation
%               ESP:        Echo spacing, ms
%               T1:         T1, ms
%               T2:         T2, ms
%               Nexc :       Nb of repetions
%               TR :         Repetition time [ms]

%
%   optional arguments (use string then value as next argument)
%
%               kmax:       maximum EPG order to include. Can be used to
%                           accelerate calculation.
%                           Setting kmax=inf ensures ALL pathways are
%                           computed
%               diff:       structure with fields:
%                           G    - Gradient amplitude(s)
%                           tau  - Gradient durations(s)
%                           D    - Diffusion coeff m^2/s (i.e. expect 10^-9)
%
%               zinit:      User specified initial state of Z0
%
%   Outputs:
%               F0:         signal (F0 state) at each echo time
%               Fn:         full EPG diagram for all transverse states
%               Zn:         full EPG diagram for all longitudinal states
%               F:          full state matrix. Each column is arranged as
%                           [F0 F0* Z0 F1 F-1* Z1 F2 F-2* Z2 ...] etc
% Nad?ge Corbin Modified 2021 early
%   Shaihan Malik 2017-07-21


%% Extra variables

for ii=1:length(varargin)
    
    % Kmax = this is the maximum EPG 'order' to consider
    % If this is infinity then don't do any pruning
    if strcmpi(varargin{ii},'kmax')
        kmax = varargin{ii+1};
    end
    
    % Diffusion - structure contains, G, tau, D
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
    end
    
    % Zinit - allow the initial state to vary, but only allow the Z0 term
    % to be different, assume all other terms are zero
    if strcmpi(varargin{ii},'zinit')
        zinit=varargin{ii+1};
    end
    
end


%% Calculation is faster when considering only appropriate EPG orders

%%% The maximum order varies through the sequence. This can be used to speed up the calculation
np = length(theta)*Nexc; % number of pulses, this includes exciatation
nPerExc = length(theta);

% if not defined, assume want max
if ~exist('kmax','var')
    kmax = 2*(np - 1);  % kmax at last echo time, this is twice # refocusing pulses
end

if isinf(kmax)
    % this flags that we don't want any pruning of pathways
    allpathways = true;
    kmax = 2*(np - 1); % this is maximum value
else
    allpathways = false;
end

%%% Variable pathways
if allpathways
    kmax_per_pulse = 2*(1:np);
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
else
    kmax_per_pulse = 2*[1:ceil(np/2) (floor(np/2)):-1:1]+1;
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
    
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
end

%%% Number of states is 3x(kmax +1) -- +1 for the zero order
N=3*(kmax+1);

%%
%%% Sort out RF pulse amplitude and phases
% split magnitude and phase
alpha = abs(theta);phi=angle(theta);

% add CPMG phase
phi(2:end) = phi(2:end) + pi/2;


%% Relaxation and shift matrices


%%% Build Shift matrix, S
S = EPG_shift_matrices(kmax);
S = sparse(S);

% Relaxation over half the echo spacing
E1 = exp(-0.5*ESP/T1);
E2 = exp(-0.5*ESP/T2);
E = diag([E2 E2 E1]);


%%% regrowth
b = zeros([N 1]);
b(3) = 1-E1;%<--- just applies to Z0

% Relaxation until TR
TRFill=TR-ESP*(nPerExc-1);
E1TR = exp(-TRFill/T1);
E2TR = exp(-TRFill/T2);
ETR = diag([E2TR E2TR E1TR]);

%%% regrowth over TR
bTR = zeros([N 1]);
bTR(3) = 1-E1TR;%<--- just applies to Z0

%%% Add in diffusion at this point
if exist('diff','var')
    E = E_diff(E,diff,kmax,N);
    ETR = E_diff(ETR,diff,kmax,N);
else
    % If no diffusion, E is the same for all EPG orders
    E = spdiags(repmat([E2 E2 E1],[1 kmax+1])',0,N,N);
    ETR = spdiags(repmat([E2TR E2TR E1TR],[1 kmax+1])',0,N,N);
end


%%% Composite relax-shift
ES=E*S;
ES=sparse(ES);

ETR=sparse(ETR);

%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

% store the indices of the top 3x3 corner, this helps build_T
i1 = [];
for ii=1:3
    i1 = cat(2,i1,sub2ind(size(T),1:3,ii*ones(1,3)));
end



%% F matrix - records the state at each echo time
F = zeros([N np-Nexc]); % there are np-Nexc echoes because of excitation pulses

%%% Initial State - can be specified by user
FF = zeros([N 1]);
if ~exist('zinit','var')
    FF(3)=1;
else
    FF(3)=zinit;
end


%% Now handle excitation
ind=0;
for exc=1:Nexc
    ind=ind+1;
    
    A = RF_rot(alpha(1),phi(1));
    build_T(A);
    
    kidx = 1:3*kmax_per_pulse(ind);
    FF(kidx) = T(kidx,kidx)*FF(kidx); %<---- state straight after excitation,
    
    
    %%% Now simulate the dephase gradient & evolution, half the readout
    FF(kidx) = ES(kidx,kidx)*FF(kidx)+b(kidx);
    
    
    %% Now simulate the refocusing pulses
    
    for jj=2:length(theta)
        ind=ind+1;
        A = RF_rot(alpha(jj),phi(jj));
        build_T(A);%<- replicate A to make large T matrix
        
        % variable maximum EPG order - accelerate calculation
        kidx = 1:3*kmax_per_pulse(ind);
        
        % Apply RF pulse to current state
        FF(kidx)=T(kidx,kidx)*FF(kidx);
        
        % Now evolve for half echo spacing, store this as the echo
        F(kidx,(jj-1)+(exc-1)*(nPerExc-1)) = ES(kidx,kidx)*FF(kidx)+b(kidx);
        % Deal with complex conjugate after shift
        F(1,(jj-1)+(exc-1)*(nPerExc-1))=conj(F(1,(jj-1)+(exc-1)*(nPerExc-1))); %<---- F0 comes from F-1 so conjugate
        
        if jj==length(theta)
            % Finally, evolve again up to next RF pulse
            FF(kidx) = ETR(kidx,kidx)*F(kidx,(jj-1)+(exc-1)*(nPerExc-1))+bTR(kidx);
            FF(1)=conj(FF(1)); %<---- F0 comes from F-1 so conjugate
            
        else
            
            % Finally, evolve again up to next RF pulse
            FF(kidx) = ES(kidx,kidx)*F(kidx,(jj-1)+(exc-1)*(nPerExc-1))+b(kidx);
            FF(1)=conj(FF(1)); %<---- F0 comes from F-1 so conjugate
        end
        
    end
end
%%% Return signal
F0 = F(1,:)*1i;

%%% Construct Fn and Zn
idx=[fliplr(5:3:size(F,1)) 1 4:3:size(F,1)];
kvals = -kmax:kmax;
%%% Remove the lowest two negative states since these are never populated
%%% at echo time
idx(1:2)=[];
kvals(1:2)=[];

%%% Now reorder
Fn = F(idx,:);
%%% Conjugate
Fn(kvals<0,:)=conj(Fn(kvals<0,:));

%%% Similar for Zn
Zn = F(3:3:end,:);


%%% NORMAL EPG transition matrix as per Weigel et al JMR 2010 276-285
    function Tap = RF_rot(a,p)
        Tap = zeros([3 3]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(4) = conj(Tap(2));
        Tap(5) = Tap(1);
        Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(7) = -1i*exp(1i*p)*sin(a);
        Tap(8) = 1i*exp(-1i*p)*sin(a);
        Tap(9) = cos(a);
    end

    function build_T(AA)
        ksft = 3*(3*(kmax+1)+1);
        for i2=1:9
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end
