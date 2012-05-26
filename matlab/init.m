function init
%% init.h
% LINES TO CHECK: 
%%

clear
load('workspace');
% loaded at start of  profit.g  eql_ma.g  welf_ma.g  and other programs %
% pp = load(configfile); % load pp=^configfile;

it=pp(IND_TYPE);
et=pp(EQL_TYPE);
[prefix,eqltype,invtype]=acronym(et,it);    % what is acronym?
sprintf('[Model: %s finvestment in %s ]',eqltype,invtype);

%  Prefix is used in file names as follows:
%   Bertrand MPNE       		--> b,
%   	monopolist 				--> mb,
%   	social planner  		--> sb,
%   Cournot  MPNE         	--> c,
%     monopolist          	--> mc,
%     social planner      	--> sc,
%   Capacity constraints MPNE --> p,
%   	monopolist 				--> mp,
%   	socialplanner   		--> sp.

small = 1e-8;        % used by profit.g and welf_ma.g when restricting p<y when using log utility (i.e. when RLG_y > 0) %

kmax=pp(KMAX);              % max value of omega; min. value is 0 
x_entryl=pp(ENTRY_LOW);     % Sunk cost of entry - lowerbound 
x_entryh=pp(ENTRY_HIGH);    % Sunk cost of entry - upperbound 
phi=pp(SCRAP_VAL);          % scrap value 
phiH=pp(RLG_SCRAP_HI);      % if >0, upperbound of random scrap.  if 0, scrap is fixed at phi 

entry_k = pp(ENTRY_AT);     % omega at which firms enter 
rlnfirms = pp(MAX_FIRMS);   % Max # of firms 
stfirm = pp(START_FIRMS);   % Firm to start at 
Beta = pp(BETA);            % Discounting factor 
delta = pp(DELTA);          % prob. that outside world moves up 
% a = pp(INV_MULT); %       % From p(x) = ax / (1 + ax) 


RLG_w0     = 1;   % if 1 then the lowest state is a non-firm --> no share and static profit is computed by removing the firm 
sprintf(' ');
if RLG_w0 == 1	 
    sprintf('Hey: in  init.h  RLG_w0 = 1  -->  lowest w is a place-holder NOT a firm offering lowest quality   ');
else
    sprintf('Hey: in  init.h  RLG_w0 = 0  -->  lowest w is actually a firm, not a placeholder for an open spot ');
end
sprintf(' ');

RLG_no_force_exit = 0;
if (pp(ENTRY_HIGH) > 1e99 && pp(SCRAP_VAL) < -1e99) || ( 1 && pp(ENTRY_HIGH) < .0001 && pp(SCRAP_VAL) < .0001)
    RLG_no_force_exit = 1;  
end
if RLG_no_force_exit;  
    sprintf('Hey: firms cannot be bumped off lowest rung since RLG_no_force_exit = 1 in init.h');  
end

RLG_wscale = pp(RLG_WSCALE);
RLG_wshift = pp(RLG_WSHIFT);
RLG_wstar  = pp(WSTAR);
RLG_sh_cap = pp(RLG_SH_CAP);
RLG_y      = pp(RLG_y);      % if >0 then   u = w - ln( RLG_y - p )  instead of  u = w - p  
RLG_alpha  = pp(RLG_ALPHA);  %  price coeff 

maxp = 99999;  % HEY: Could have pp(RLG_MAXP) be negative when a factor of max duop price and positive when the maxp itself, as computed in profit.g

% RLG_leap   = pp(RLG_LEAP); %  % Prob entrant goes to KMAX instead of ENTRY_AT 
RLG_leap   = pp(34);        % until restart batch_ma.g jobs  since RLG_leap defined by pmg.h which is NOT run inside for loops in batch_*.g %


if RLG_y>0 && RLG_wstar>0    
    sprintf(' ');
    sprintf(' should not use both RLG_y and RLG_wstar  ');
    sprintf(' ');
end
if RLG_y>0 && RLG_y<pp(MC)
    sprintf(' ');
    sprintf(' must have RLG_y > MC else profits = zero ');
    sprintf(' ');
end

RLG_out = pp(RLG_OUTGOOD);
if RLG_out
    sprintf('Outside good is turned ON');
else
    sprintf('Outside good is turned OFF.  Monopolist price is RLG_y (if RLG_y>0) or RLG_MAXP * leader price in max differentiated duopoly');
end

% aeff = (pp(INV_MULT)/pp(KMAX)) .* rev(seqa(1,1,pp(KMAX)+1)); 
if pp(RLG_INV)>0
  if  0  && BatchJob<3600   % NEED TO CHECK THIS regarding the 0 part
    % aeff = pp(INV_MULT)+ pp(RLG_INV)*rev(seqa(0,1,pp(KMAX)+1)) ; %   % ^2 % 
    aeff = pp(INV_MULT)*(1+pp(RLG_INV)*rev(seqa(0,1,pp(KMAX)+1)));     % proportional % % WHAT IS SEQA
    sprintf('Using linear spillover with coeff of %f  , yielding aeff = %f',fpp(RLG_INV),aeff);
  else
    aeff = ones(pp(KMAX),1)*pp(INV_MULT)*[(1+pp(RLG_INV)),pp(INV_MULT)];     % all laggards at same aeff, which is some factor of frontier aeff %
    sprintf('Using common spillover for all laggards, yielding aeff = %f',aeff);
    if delta>0 || RLG_wstar>0
        sprintf('HEY: probably should NOT use this form of spillover (set in init.h), since leader often NOT at frontier');  
    end 
  end
  % NOTE: something must be wrong in planner (and possibly monop) solution       %
  %       since Nash yields HIGHER CS+PS when using wstar = 0 (ie. GG unbounded) %
else
  aeff = pp(INV_MULT)*ones(pp(KMAX)+1,1);
end

RLG_minx1 = 0;      % ax = (1+ax) m  --> ax(1-m) = m -->  x = m/a(1-m) % 
if pp(RLG_MIN_INNOV1)>0  
    RLG_minx1 = pp(RLG_MIN_INNOV1)/((1-pp(RLG_MIN_INNOV1))*aeff(rlnfirms+1))
end


% outwidth 256;
save('init');
save('workspace');
