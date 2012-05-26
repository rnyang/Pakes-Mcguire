function welf_ma
%% welf_ma.g
% LINES TO CHECK: 
%%

% This version is markwds.prg
% Written by: Gautam Gowrisankaran
% April 25, 1993
% This program generates the welfare output programs used for the
% Markov-perfect game. %

% new; %   % new; must be first line in program, else it terminates the program.  Hence, my idea to retain interactive PM code by running  new;  only when BatchJob==0 does not work %
    load('workspace');
    if BatchJob==0  
        load pmg.mat;   
    end

    % clrscr();     
    disp '**** Computing Welfare Statistics ****';

    load init.mat;
    M  = pp(MKT_SIZE);
    mc = pp(MC);

    rand('seed', 28010);

    wstart = zeros(rlnfirms,1);

    if BatchJob==0
      disp "Enter initial efficiency level(s)";
      disp "Valid range for parameter values is 0 - " compact(kmax);
      i=1;
      while i<=rlnfirms
         wdef=iif(i==1, entry_k+2, 0);
         prompt=['Firm ',compact(i)];
         wstart(i)=getint(prompt, 0, kmax, wdef);
         i=i+1;
      end
      numtimes=getint('Number of periods to simulate (1-1K, default 100)',1,1000,100);
      numruns=getint('Number of runs to simulate for (1-100K, default 1000)',1,100000,1000);
    else
      wstart(1) = kmax;                   % AFTER load isentry(), wstart(2) = pp(ENTRY_AT) if isentry>0 with this monop  %
      if pp(WSTAR)>0 || delta>0,         wstart(1) = min(kmax|pp(ENTRY_AT)+2);  end   % alt: min(kmax-1|(ceil((6+pp(WSTAR)-pp(RLG_WSHIFT))/pp(RLG_WSCALE))))  where WSTAR is on post shift/scale grid, but wstart is on 0,1,2...Kmax grid %
      if pp(MAX_FIRMS)==2 || RLG_y==0,   wstart(2) = pp(ENTRY_AT);               end   % fixed DUOP, so start with DUOP, || linear model(since want to avoid monopoly) %
      if (pp(ENTRY_HIGH)<1e-8 || pp(ENTRY_HIGH)>1e100) && pp(WSTAR)<=0,      wstart = kmax+zeros(rlnfirms,1);  end    % NO ENTRY/EXIT --> start with all firms at frontier %
      numtimes= 100;   if pp(WSTAR)>0 || BatchJob<-1,   numtimes= 500;  end  % increase numtimes if WSTAR to remove uncertainty of whether wstart(1) is near wfinal(1).  If WSTAR<=0  wstart(1) = wfinal(1) = kmax %
      numruns=  1000;
    end

    % jobstart = 5; %  % Start of job creation statistics.  What is this variable used for? %

    % Set up binomial coefficients for decoding/encoding of n-tuples %
    if pp(WSTAR)==0, kmax = kmax+1;  end         % need tables 1 beyond kmax %
    binom = eye(rlnfirms+kmax+1);
    binom = [zeros(rlnfirms+kmax+1,1),binom];
    i=2;
    while i <= rlnfirms+kmax+1
      binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
      i=i+1;
    end

    if pp(WSTAR)==0, kmax = kmax-1;  end

    wmax = binom(rlnfirms+kmax+1,kmax+2);     % Number of possible industry structures %

    %
%     FixedX = 3.5;
%     while FixedX < 5 % WHERE DOES THE LOOP END?
%     FixedX = FixedX + .2;
%     x = FixedX + 0*x;  p = aeff*x;  p = p./(1+p);
%     printd(' ', 'hardcoding fixed policy of x=  ', x(1), '  yielding ax/(1+ax)=  ', p(1) , ' ');
    %

    % Load in all the data stored by the equilibrium generation program %
    % This data is: v (value), x (investment), p (probability of state rising), isentry %

    filename = [prefix , 'markov.' , num2str(rlnfirms) , 'ot'];

    % filename = "BJ/" $+ filename; %      % to load files copied to some directory, e.g. BJ/, to avoid recomputing eqm %

    % load bigread() = ^filename;
    bigread = load(filename);
    v = (reshape(bigread(1:wmax*rlnfirms),wmax,rlnfirms));   % value functions   %
    bigread = bigread(wmax*rlnfirms+1:rows(bigread));
    x = (reshape(bigread(1:wmax*rlnfirms),wmax,rlnfirms));   % investment policy %
    bigread = bigread(wmax*rlnfirms+1:rows(bigread));
    p = (reshape(bigread(1:wmax*rlnfirms),wmax,rlnfirms));   % prob(own w up)    %
    if phiH ~= 0                            % random scrap --> load exit policy %
      bigread = bigread(wmax*rlnfirms+1:rows(bigread));
      isexit = (reshape(bigread(1:wmax*rlnfirms),wmax,rlnfirms));    % prob exit %
    end
    bigread = bigread(wmax*rlnfirms+1:rows(bigread));
    isentry = bigread(1:wmax);                                      % prob entry %
    bigread = bigread(wmax+1:rows(bigread));
    v_entry = bigread(1:wmax);                                  % value of entry %

    % Load in all data from the static profit calculation %
    % The data is: firm profits, consumer surplus, market shares at each state,
    % price/cost margins, one-firm concentration ratios. %
    % filename = prefix $+ "cons." $+ num2str(rlnfirms,"%*.*lf",1,0) $+ "f";

    % filename = "BJ/" $+ filename;  %

    bigread = load(filename); % CHECK LATER

    i = wmax;  if pp(WSTAR)==0,   wmax = binom(rlnfirms+2+kmax,kmax+3);  end    % RLG_wstar==0 extended kmax by one in profit.g %

    profit   = reshape(bigread(1:wmax*rlnfirms),wmax,rlnfirms);
    bigread  = bigread(wmax*rlnfirms+1:rows(bigread));
    csurplus = reshape(bigread(1:wmax),wmax,1);
    bigread  = bigread(wmax+1:rows(bigread));
    share    = reshape(bigread(1:wmax*rlnfirms),wmax,rlnfirms);
    bigread  = bigread(wmax*rlnfirms+1:rows(bigread));
    pmcmargm = reshape(bigread(1:wmax),wmax,1);
    bigread  = bigread(wmax+1:rows(bigread));
    concentm = reshape(bigread(1:wmax),wmax,1);
    bigread  = bigread(wmax+1:rows(bigread));
    price    = reshape(bigread(1:wmax*rlnfirms),wmax,rlnfirms);

    wmax = i;  % restore wmax %

    if BatchJob>0 && isentry(encode(wstart))>0,    wstart(2) = pp(ENTRY_AT);    end
    wstart = sort(wstart,1);
    wstart = wstart(end:-1:1);

%     format /rd 3,0;
    disp(outline('Simulation'));
    disp 'Starting conditions: (" wstart' ")';
    % format /rd 6,2;
    disp 'Number of simulation periods: '; disp( compact(numtimes));
    disp '; Number of runs: '; disp(compact(numruns));

    if pp(WSTAR)==0  && 0 
      kmax = kmax-round(2+pp(KMAX)/6*(1-delta));
      printd(' ', '** Adjusting kmax down by ', round(2+pp(KMAX)/6*(1-delta)), ' (to" kmax ") in welf_ma since using WSTAR = 0 is more accurate when simulate with kmax lower than used to get investment policy', ' ');
    end
    wmax = binom(rlnfirms+kmax+1,kmax+2);

    consurp  = zeros(numruns,1); % (Total) consumer surplus %
    prodsurp = zeros(numruns,1); % Producer surplus %

    prodsurp2= zeros(numruns,1); % PS with unbounded outside utility using reoptimized prices  %
    consurp2 = zeros(numruns,1); % CS with unbounded outside utility using reoptimized prices  %

    totinnov = zeros(numruns,1);  % other RLG stats   %
    totinnov1= zeros(numruns,1);  % leader innovation %
    totinnov2= zeros(numruns,1);  % firm 2 innovation %
    totinnovT= zeros(numruns,1);  % frontier innov T  %
    totinvest= zeros(numruns,1);  % industry invest   %
    totinves1= zeros(numruns,1);  % leader   invest   %
    totinves0= zeros(numruns,1);  % industry invest 0 && no entry || exit --> Absorbing State %
    totinve0n= zeros(numruns,1);  % n firms invest 0  %
    totprof1 = zeros(numruns,1);  % leader profit excluding investment %
    totprof2 = zeros(numruns,1);  % firm 2 profit excluding investment %
    totnfirms= zeros(numruns,1);
    tot_firms= zeros(numruns,rlnfirms+1);  % last col is freq of empty industry %
    totentry = zeros(numruns,1);
    totshare0= zeros(numruns,1);
    totshare1= zeros(numruns,1);
    totshare2= zeros(numruns,1);
    totprice1= zeros(numruns,1);
    totprice2= zeros(numruns,1);
    totw2    = zeros(numruns,1);
    totw2T   = zeros(numruns,1);
    totshar1T= zeros(numruns,1);
    totmarkup= zeros(numruns,1);
    totwstar = zeros(numruns,1);
    totwsame = zeros(numruns,1);
    totLRstay= zeros(numruns,1); % a firm at Lowest Rung stays with some probability %
    totLRx   = zeros(numruns,1); % a firm at Lowest Rung invests > 0 %
    totpdiff = zeros(numruns,1); % static p + RLGforced  -  newton p %

    totexit  = zeros(numruns,1); % total number of firms that exit, && therefore contribute to totv0_vx %
    totv0    = zeros(numruns,1); % initial v %
    totvx    = zeros(numruns,1); % realized profits %
    totv0i   = zeros(numruns,1); % initial v, initial leading firm only %
    totvxi   = zeros(numruns,1); % realized profits, initial leading firm only %
    totlspi  = zeros(numruns,1); % lifespan, initial leading firm only %
    totxT0   = zeros(numruns,1); % # runs ending in a state with zero investment (absorbing state if delta=0) %
    discEU_T = zeros(numruns,1); % discounted EU from final period (check if small as needed for run to be long enough && CS to be well-defined) %
    discCS_T = zeros(numruns,1); % discounted CS from final period (check if small as needed for run to be long enough && CS to be well-defined) %
    hitkmax  = zeros(numruns,1); % # times lead firm at kmax %
    maxk     = 0;                % maximum k achieved across all runs %
    x1absorb = 0;                % leader investment in absorbing state in which leader invests to deter laggards' investments %

    % from ds_ma.g %
    active   = zeros(rlnfirms+1,1);  % No. of periods with n firms active. %
    exitors  = 0;  % No. of exitors %
    entrants = 0;  % No. of entrants %
    entexit  = 0;  % No. of periods with both entry && exit %
    lifedis  = 0;  % Distribution of firm lifespans %
    valuedis = 0;  % Distribution of firm's total profits %

    RLGforced = 0;
    RLGinnov  = 0;

    simfile = ['BJ' , num2str(BatchJob) , '_sim'];
%     output file = ^simfile reset;   screen off;    output on;
    file = fopen(strcat('simfile','.mat'),'wt');
    fprintf(file,  [nr(1), t(2),  cum_innov(3),  pr_innov(4),  pr_entry(5),   w,  thisexit,   ax/(1+ax),   share,  price,  profit]);
%     screen on;     output off;
    fclose(file);

    nr = 1;
    while nr <= numruns
      RLGforced = 0;
      RLGinnov  = 0;
      discEU    = zeros(numtimes,1);      % discounted EU : used to compute CS when RLG_out = 0 && RLG_y > 0 %
      discCS    = zeros(numtimes,1);      % discounted CS : used to compute CS when RLG_out = 0 && RLG_y > 0 %

      firm_id  = 1;
      if rlnfirms>1
        firm_id  = [1;zeros(rlnfirms-1,1)];       % track initial value && realized values of 1st leading firm %
      end

      lifemx   = zeros(rlnfirms,1);  % No. of periods that firm has been active for %
      valuemx  = zeros(rlnfirms,1);  % Total profits of each firm to date %
      valuem0  = zeros(rlnfirms,1);  % value of each firm's initial state %
      wthis    = wstart;
      wthisprev= wstart;
      lifemx   = lifemx + (wthis(wthis > 0));

      ii = lifemx(lifemx==1);  % initialize with firms that already exist prior to 1st simulated period %
      if sum(ii)>0;  
        valuem0 = valuem0.*(1-ii) + ii.*(v(encode(wthis),:)');
      end

      % "nr " nr "  valuem0 " valuem0'  "  firm_id " firm_id' %

      t = 0;
      while t < numtimes

        codew2 = encode(wthis);    % Get probabilities of investment causing a rise in eff level & actual investment && value function for this state %
        prob = p(codew2,:)';       % prob(own w up) %
        xx = x(codew2,:)';         % Investment decisions made prior to observing entry/exit outcomes %
        vv = v(codew2,:)';         % Profit should also be based on codew2 since also before outcomes %

        if (t+1==numtimes) && (mod(nr , (numruns/10)) == 0)     % disp last state for 10 of the runs -- shows progress in completing numtimes, && gives limited sample of distribution of terminal state %
          disp "Runs computed: "; disp(padr(nr,5,0)); 
    %       format /m1 /rd 3,0;  
          printd('  w_T= ', wthis');
    %       format /m1 /rd 8,4;  
          printd('  invest= ', xx', '  innov= ', prob',  '  entry= ', isentry(codew2),  '  lifemx= ', lifemx', '  valuemx= ', valuemx', '  valuem0= ', valuem0',  '  RLGforced= ', padr(RLGforced,3,0),  '  discEU avg~final= ', meanc(discEU), padr(discEU(numtimes),10,8));   
        end

        if sum( (wthis==1) & xx(xx>1e-8) )>0                     totLRx(nr)    = totLRx(nr)   +1;   end
        temp = isexit(codew2,:)';
        if sum( (wthis==1) & temp(temp<1-1e-8))>0     totLRstay(nr) = totLRstay(nr)+1;   end

        ii = sum(wthis>0);    
        if ii==0  
          tot_firms(nr,rlnfirms+1)= tot_firms(nr,rlnfirms+1) + 1;
          if sum(tot_firms(1:nr,rlnfirms+1)) < 5
        printd(' \n',  'HEY: no firms are in the industry !!   t= ',  t  ,' wthis= ', wthis' ,' wthisprev= ' ,wthisprev' ,' at row ', encode(wthisprev),  ' in policy functions');
          end
        else  
          tot_firms(nr,ii)= tot_firms(nr,ii) + 1;                  % based on pre- entry/exit/investment   %
          totnfirms(nr)   = totnfirms(nr) + ii;                    % before exit & entry since profit & CS %
          totinve0n(nr)   = totinve0n(nr) + sum(xx(xx(1:ii)==0) );  % pre-exit to track deadbeat firms      %    % since totinve0n focuses on laggards, ignore RLG_minx1 as the "zero" for frontier firms %
        end

        totinves1(nr)= totinves1(nr) + xx(1);            % if xx(1)==0;  "Hey: leader invests 0 in welf_ma.g";  end  %
        totinnov1(nr)= totinnov1(nr) + prob(1);
        totinnov2(nr)= totinnov2(nr) + prob(2);
        totprof1(nr) = totprof1(nr)  + profit(codew2,1);
        totprof2(nr) = totprof2(nr)  + profit(codew2,2);

        totw2(nr) = totw2(nr) + wthis(1)-wthis(2);       % leader advantage %

        maxk = max([maxk;wthis(1)]);
        if wthis(1)==kmax,  hitkmax(nr)  = hitkmax(nr)+1;   end
        if RLG_wstar>0
            temp = (wthis.*RLG_wscale+RLG_wshift);
            totwstar(nr) = totwstar(nr) + sum(  temp(temp > RLG_wstar));  
        end
        if rlnfirms>1,      totwsame(nr) = totwsame(nr) + (wthis(1)==wthis(2));  end

        % Find out which firms want to leave. %
        wtrans = zeros(rlnfirms,1);
        if phiH==0                                           % fixed scrap %
          ii = (min(vv) == phi)*(minindc(vv)-1) + (min(vv) > phi)*rlnfirms;
          if ii > 0,  wtrans(1:ii) = wthis(1:ii);   end    % ii = #firms remaining %
        else                                                 % random scrap %
          prexit = rand(rlnfirms,1);
          wtrans = wthis .* (prexit(prexit > isexit(codew2,:)'));    % exit if prexit <= isexit %
        end

        % Now figure out exit. Must capture people who exit voluntarily,
    %     as well as firms whose efficiency has been driven down to zero. %
        lifemx = lifemx + (wthis > 0);
        thisexit = (wtrans == 0) & (lifemx > 0);    % checking lifemx ensures the exitor is a real firm, not a vacancy %
        lifemx = lifemx - (wthis > 0);

        if t<0
            printd('t= ', t, ' thisexit ', thisexit',  '  isexit ',  isexit(codew2,:),  '  wthis ', wthis',   '  wtrans ', wtrans',  ' firm_id ', firm_id',  ' totexit ', totexit(nr));
        end

        phis = zeros(rlnfirms,1);       % vector of realized phi values.  Want simulated phis, not E(phi|exit) %
        if sum(thisexit) > 0
          % "Exit in period " t "firms" selif(firm_id, thisexit); %
          ii = indexcat(thisexit,1);
          if phiH==0
              phis(ii) = phi*ones(rows(ii),1);             % ii = exitor indices %
          else
              phis(ii) = phiH - (phiH-phi)*prexit(ii);     % prexit are the rand() draws compared to isexit() above %
          end

          xx   =   xx .* (1-thisexit);    % zero out exitors xx && prob since those values are conditional on NOT exiting %
          prob = prob .* (1-thisexit);    % only relevant for random scrap case, but harmless to do for all cases %
        end


        % codew = encode(wtrans); %      % HEY: why does PM original code use wtrans (i.e., AFTER EXIT market structure)?? %
        codew = codew2;                  % eql_ma.g uses profit(w,j) where w is BEFORE exit, so RLG uses codew2 NOT encode(wtrans) %

        sigma = share(codew,:)';       % may get overwritten later %

        % Now, tally the statistics.  RLG moved this block from post-entry determination to here (since profits are BEFORE entry/exit/invest outcomes %

        consurp(nr) = consurp(nr) + (Beta^t) * csurplus(codew);                                      % if RLG_y > 0 will convert to true CS after last period  %
        consurp2(nr)= consurp2(nr)+ (Beta^t) *(csurplus(codew)+M*RLGforced*RLG_wscale/RLG_alpha);    % add extra CS (or M/alpha scaled EU) if grid has shifted %

        discEU(t+1) =               (Beta^t) *(csurplus(codew)*RLG_alpha/M + RLGforced*RLG_wscale);  % enables checking whether EU && CS from last period is small after being discounted to 0 %
        if RLG_y
            discCS(t+1) = (Beta^t) * M * ( exp( (discEU(t+1)/(Beta^t))/RLG_alpha ) - RLG_y );
        else
            discCS(t+1) = (Beta^t) *(csurplus(codew)+M*RLGforced*RLG_wscale/RLG_alpha);
        end

        i = sum(wthis(wthis>0));       % #active firms  %
        if i==0                  % EMPTY industry %
          prodsurp(nr) = prodsurp(nr)  + (Beta^t)*(sum(phis) - sum(xx));
          prodsurp2(nr)= prodsurp2(nr) + (Beta^t)*(sum(phis) - sum(xx));
          valuemx = valuemx + (Beta^(lifemx-1)).*( -xx);  % NOTE: firm that exists in wstart has lifemx = 1 when t=0, so use lifemx-1 %

          totinves0(nr)= totinves0(nr) + (isentry(codew)==0);   % absorbing state %
          totinvest(nr)= totinvest(nr) + sum(xx);              % consurp2, totmarkup, totshare0, are unchanged %
          totshare0(nr)= totshare0(nr) + 1;

        else     % industry NOT EMPTY %
          www   = (RLGforced + wthis(1:i))*RLG_wscale + RLG_wshift;    % inside qualities relative to TRUE outside quality %
          sigma = share(codew,1:i)';   % by NOT recomputing prices (if unbounded) we are assuming KMAX high enough that prices independent of RLGforced %
          ppp   = price(codew,1:i)';   % (i.e., outside good is inconsequential if more than 1 firm) %

          if  0  && RLG_wstar==0 && RLG_y==0     % unbounded state-space && linear utility -->  want to RECOMPUTE price via Newton -->  BUT NOT USING unbounded in eql_ma.g anymore, so do not do this anymore  %
        % I ran this code && confirmed approximation used below for monopolist with linear utility is accurate.  Get same answer, but slower to re-newton, so SKIP %
        if i==1                                % note, but no alarm:  RLG_MAXP in profit.g may have constrained monopolist's price %
          ppp = price(codew,1) + RLGforced*RLG_wscale/RLG_alpha ;
        else
          www = www - max([0;min(www)-30]);     % avoids singular deriv in newton, && prices should still be correct %
          ppp = price(codew,1:i)';              % initialize ppp search at offline static price %
          ppp = newton2(ppp,@cfunk2);
        end
        sigma = exp(www - RLG_alpha*ppp);
        sigma = sigma./( RLG_out + sum(sigma) );

        totpdiff(nr)= totpdiff(nr) + meanc( price(codew,1:i)' - ppp );

        if RLGforced==0 && max(abs( price(codew,1:i)' - ppp ))>0.015
%           format /m1 /rds 4,2;  
          pritd(  't= ', padr(t,4,0), '.  RLGforced = 0 so should get  price(codew) = ppp, but wtrans=', wtrans', ' has price(): ', price(codew,1:i), '  ppp: ',  ppp');
        end

          else  % use prices ppp && shares sigma as computed by profit.g %

        % an unbounded monopolist with linear utility consumers raises price to maintain same utility %
            if RLG_y==0 && RLG_wstar==0 && rlnfirms==1,   ppp = price(codew,1)+RLGforced*RLG_wscale/RLG_alpha;    end

        www = M * (ppp-mc).*sigma - pp(FIXED_COST);     % current period profit %
        prodsurp2(nr)= prodsurp2(nr) + (Beta^t)*(sum( www )             + sum(phis) - sum(xx)) ;      % entryfee subtracted later %
        prodsurp(nr) = prodsurp(nr)  + (Beta^t)*( sum(profit(codew,:)') + sum(phis) - sum(xx)) ;

        if RLG_wstar == 0
            valuemx(1:i) = valuemx(1:i) + (Beta^(lifemx(1:i)-1)).*(www - xx(1:i));     % NOTE: firm that exists in wstart has lifemx = 1 when t=0, so use lifemx-1 %
        else
            valuemx = valuemx + (Beta^(lifemx-1)).*( profit(codew,:)' - xx );          % NOTE: firm that exists in wstart has lifemx = 1 when t=0, so use lifemx-1.  Enables monopolist valuemx to match its prodsurp %
        end
        totinvest(nr)= totinvest(nr) + sum(xx);
        if delta==0
          totinves0(nr)= totinves0(nr) +(sum(xx)==0 && isentry(codew)==0 &&  sum(isexit(codew,:)'.*wthis)==0 );   % Absorbing State %

          if RLG_no_force_exit && xx(1)>0 && max(wthis(2:rlnfirms))<=1   %  Can also have absorbing state with Sole Frontier x > 0 && all other firms at 1 || out && all entry, exit = 0, IF RLG_no_force_exit = 0 %
            totinves0(nr)= totinves0(nr) + (sum(xx(2:rlnfirms))==0 && isentry(codew)==0 &&  sum(isexit(codew,:)'.*wthis)==0 );   % Absorbing State, but with investment %
            x1absorb = xx(1);
          end
        end

        totmarkup(nr)= totmarkup(nr) + sum(ppp.*sigma) / mc / sum(sigma);
        totshare0(nr)= totshare0(nr) + 1-sum(sigma);
        totshare1(nr)= totshare1(nr) + sigma(1);
        totprice1(nr)= totprice1(nr) + ppp(1);
        if i>1 
          totshare2(nr)= totshare2(nr) + sigma(2);  
          totprice2(nr)= totprice2(nr) + ppp(2);  
        end
        % "nr~t= " nr t "  wthis " wthis' "  wtrans " wtrans' "  m0 " valuem0' "  mx " valuemx'  "  ii " ii'  "  v " v(encode(wthis),:) "  profit " profit(codew,:) "  xx " xx' ;  %

          end

        end  % end of block moved by RLG to precede exit instead of follow entry stuff %


        if sum(thisexit) > 0

          valval = valuemx(thisexit)+(Beta^(lifemx(thisexit)).*phis(thisexit));   % RLG added profit() here since exit occurs AFTER profits %

          if ( firm_id' * thisexit )>0               % the leading firm exited %
        i = indexcat(firm_id,1);
        totlspi(nr)= lifemx(i);                   % mean totv0i = mean totvxi  if value function correct %
        totv0i(nr) = valuem0(i);
        totvxi(nr) = valuemx(i) + Beta^lifemx(i)*phis(i);
        % "t = " t " firm_id " firm_id' "  valuem0 " valuem0'  "  lifemx " lifemx'  "  firm_id==1 "  (firm_id.==1)'  " thisexit " thisexit'; %
        firm_id = zeros(rlnfirms,1);              % no longer have the initial leading firm %
          end

          % lifedis  = lifedis|selif(lifemx,thisexit); %
          % valuedis = valuedis|valval; %
          totexit(nr)= totexit(nr) + sum(thisexit);
          totv0(nr)  = totv0(nr) + sum( valuem0(thisexit)  );
          totvx(nr)  = totvx(nr) + sum( valval );
          lifemx     = lifemx  .* (1-thisexit);
          valuemx    = valuemx .* (1-thisexit);
          valuem0    = valuem0 .* (1-thisexit);
        end

        entrypr = rand(1,1);           % original PM used  encode(wtrans), which is AFTER exit occurs, in place of codew below %
        if entrypr < isentry(codew)   % with random scrap, entry/exit DECISIONS both depend on wthis NOT wtrans (i.e., BEFORE either entry || exit is observed) %
          % prob(rlnfirms)  = 0; %

            % when nfirms > 2, only do leap if TWO open spots so can still check frequency of full industry for need to increase nfirms %

          doleap = (RLG_leap>0 && rlnfirms>1);   if doleap && rlnfirms>2,   doleap = wthis(rlnfirms-1)==0 ;   end
          if doleap && rand(1,1) < RLG_leap  
        wtrans(rlnfirms)  = kmax;
          else                  
        wtrans(rlnfirms)  = entry_k;
          end
          entryfee          = x_entryl + entrypr * (x_entryh - x_entryl);      % uses realized entryfee, NOT E(fee|enter) % 
          totentry(nr)      = totentry(nr)+1;
          valuemx(rlnfirms) = 0;
          valuem0(rlnfirms) = v_entry(codew);
          % "nr~t= " nr t "  wthis " wthis' "  wtrans " wtrans' "  m0 " valuem0' "  mx " valuemx' "  v_entry " v_entry(i);  %

          prodsurp(nr)  = prodsurp(nr)  - (Beta^t)*entryfee;    % PM subtracted entryfee earlier %
          prodsurp2(nr) = prodsurp2(nr) - (Beta^t)*entryfee;    % RLG moved here so that profit earned before entry/exit/invest %
        end    % entry %

        lifemx = lifemx + (wtrans > 0);  % moved here from before profit tallied by RLG since want lifemx = 0 for new entrant %

        % RLG moved wnext stuff AFTER computing surpluses since do not want to update RLGforced until after computed CS %
        wnext = wtrans+(prob(prob>=rand(rlnfirms,1)));
        if rand(1,1) <= delta 

          if RLG_no_force_exit;       wnext = wnext-1 + (wnext==1);          % Never bump firms off lowest rung of ladder %
          else                       wnext = wnext-1;
          end

          RLGinnov = RLGinnov + 1;
        end

        if RLG_wstar <= 0 && sum(wnext==kmax+1)>0;  % RLG statespace && frontier moving beyond kmax, so shift ALL down %

          if RLG_no_force_exit;       wnext = wnext-1 + (wnext==1);          % Never bump firms off lowest rung of ladder %
          else                       wnext = wnext-1;
          end

          RLGforced = RLGforced+1;
        end

        wnext = max(([wnext,zeros(rlnfirms,1)])');
        wnext = min(([wnext,kmax*ones(rlnfirms,1)])');  % RLG added %

        if nr<101;  % save to file for plotting by  simplots.m :  nr  t  cum_innov  pr_innov  pr_entry  w  thisexit  ax/1+ax  share  price  profit %
%           output file = ^simfile;   screen off;    output on;           % pr_innov uses investment ignoring exit,   fine since leader rarely exits %
          file = fopen(strcat(simfile,'.mat'),'at');
%           format /m1 /rds 1,0;  
          fprintf(file, sprintf('%1.0f %1.0f %1.0f', nr,  t+1,  RLGforced+RLGinnov+wthis(1)-wstart(1)));
%           format /m1 /rds 6,4;  
          fprintf(file, sprintf('%6.4f %6.4f', (1-(1-p(codew,1))^sum(wthis(wthis==wthis(1)))),  isentry(codew))) ;
%           format /m1 /rds 1,0;  
          fprintf(file, sprintf('%1.0f %1.0f', wthis', thisexit'));
%           format /m1 /rds 6,4;  
          fprintf(file, sprintf('%6.4f %6.4f %6.4f %6.4f', p(codew,:),  share(codew,:),  price(codew,:),  profit(codew,:))) ;
%           screen on;     output off;
          fclose(file);
        end

        % Now re-sort all firm level data, to reflect the fact that firms must be in descending order next year %
        wthisprev = wthis;
        temp = sort([wnext,lifemx,valuemx,valuem0,firm_id],1);
        temp = temp(end:-1:1);
        wthis = temp(:,1);  lifemx = temp(:,2);  valuemx = temp(:,3);  valuem0 = temp(:,4);  firm_id = temp(:,5);

        % if wthis(2)==1;  " t= " t " wthis " wthis';  end %

        t = t+1;
      end

      % "firm_id " firm_id' " lifemx %

      % get realized + cont.value of firms active in last simulated period %
      thisexit = (lifemx>=1);  % RLG confirmed valuemx = valuem0 exact when wstart is an absorbing state (so discounting powers correct) %
      if sum(thisexit)>0;      % RLG also notes that mean totvx usually > mean totv0 ... probably due to value driven by outliers ? %
        temp = lifemx-1;
        valval   = valuemx(thisexit) + Beta^( temp(thisexit)).*v(encode(wthis)', thisexit);
        totexit(nr)  = totexit(nr)  + sum(thisexit);
        totv0(nr) = totv0(nr) + sum( valuem0(thisexit)  );  % mean totv0 = mean totvx  if value function correct %
        totvx(nr) = totvx(nr) + sum( valval );
      end

      % "nr~t= " nr t "  firm_id " firm_id' " totlspi(nr) " totlspi(nr) "  wthis " wthis' "  wtrans " wtrans' "  m0 " valuem0' "  mx " valuemx'  "  thisexit " thisexit'  "  v " v(encode(wthis),:) "  valval " valval' "  lifemx= " lifemx'; %

      if sum( firm_id==1 );
        totv0i(nr) = valuem0(firm_id==1);
        totvxi(nr) = valuemx(firm_id==1) +  Beta^numtimes * v(encode(wthis)', firm_id==1);
        totlspi(nr)= lifemx(firm_id==1);  % will be numtimes+1 %
      end

      % totxT0(nr) = (sum(xx)-RLG_minx1*sum(wthis.==kmax))<1e-15; %   % investment in terminal state is 0 --> an absorbing state if delta = 0 && entry && exit also 0 %
      totxT0(nr) = totinves0(nr)>0;  % any periods hit absorbing state %

      totinnov(nr) = RLGforced+RLGinnov+wthis(1)-wstart(1);

      if RLG_y > 0
        % Conceptually 2 ways to compute dynamic CS when utility has log(y-p): see discussion in profit.g %
        % discCS(t+1)= (Beta^t) * M * exp( ln( sum( exp( (RLGforced + wthis(ii)).*RLG_wscale+RLG_wshift + RLG_alpha * ln( RLG_y - price(codew,ii)' ) )))/RLG_alpha); %
        % alpha*ln(y)/(1-Beta) = consurp2 = cumsum( discounted EU ) --> y as a per-period compensating income = exp( (1-Beta)*consurp2/alpha ) --> Y = y/(1-Beta) is the discounted compensating y.  %
        % below is only valid if periods per run is large, since using 1/(1-Beta) which is infinite sum %
        % meanc(discEU())*Beta^numtimes/(1-Beta)) assumes all periods BEYOND numtimes have this EU, which is still too low, but better than assuming EU = 0 after numtimes %

        www = (1-Beta)*consurp(nr)*RLG_alpha/M;                      % find steam of income that generates this per-period equivalent EU to lifetime discounted EU %
        consurp(nr) = M * ( exp( www/RLG_alpha ) - RLG_y) /(1-Beta);

        www = (1-Beta)*consurp2(nr)*RLG_alpha/M;
        consurp2(nr)= M * ( exp( www/RLG_alpha ) - RLG_y) /(1-Beta);

      end

      discEU_T(nr) = discEU(numtimes);
      discCS_T(nr) = discCS(numtimes);

      if sum(consurp<0)>0,   disp 'HEY: in welf_ma.g CS <0,      resetting to 0';  consurp = max([zeros(1,numruns);consurp']);    end
      if sum(consurp2<0)>0,  disp 'HEY: in welf_ma.g CS2<0,      resetting to 0';  consurp2= max([zeros(1,numruns);consurp2']);   end
      if sum(discCS_T<0)>0,  disp 'HEY: in welf_ma.g discCS_T<0, resetting to 0';  discCS_T= max([zeros(1,numruns);discCS_T']);   end


      totw2T(nr)    = wthisprev(1)-wthisprev(2);   % leader advantage in last period %
      totshar1T(nr) = sigma(1);                    % leader share in last period %
      totinnovT(nr) = 1-(1-prob(1))^sum(wthisprev==wthisprev(1));       % prob frontier innovates last period.  Assumes no firms at frontier will exit. %

      nr=nr+1;
    end

    printd(' \n', 'Total Enter ~ Exit (includes firms alive at end) ~ Exit-Entry-initfirms = ', padr(sum(totentry),1,0) ,   padr(sum(totexit),1,0), '  ', padr(sum(totexit-totentry)-numruns*sum(wstart(wstart>0)),1,0));

%     after_loops: %% USE WHILE
%     while after_loops == 1
%         after_loops = 0;

    totinnov = totinnov /numtimes ;
    totinnov1= totinnov1/numtimes ;
    totinnov2= totinnov2/numtimes ;  % uses 0 for states in which only 1 firm %
    totinvest= totinvest/numtimes ;
    totinves1= totinves1/numtimes ;
    totinves0= totinves0/numtimes ;
    totinve0n= totinve0n/numtimes ;
    totnfirms= totnfirms/numtimes ;
    tot_firms= tot_firms/numtimes ;
    totentry = totentry /numtimes ;
    totshare0= totshare0/numtimes ;
    totshare1= totshare1/numtimes ;
    totshare2= totshare2/numtimes ;  % uses 0 for states in which only 1 firm %
    totprof1 = totprof1 /numtimes ;
    totprof2 = totprof2 /numtimes ;  % uses 0 for states in which only 1 firm %
    totprice1= totprice1/numtimes ;
    totprice2= totprice2/numtimes ;  % uses 0 for states in which only 1 firm %
    totmarkup= totmarkup/numtimes ;
    totwstar = totwstar /numtimes ;
    totwsame = totwsame /numtimes ;
    totLRstay= totLRstay/numtimes ;
    totLRx   = totLRx   /numtimes ;
    totw2    = totw2    /numtimes ;  % uses 0 for states in which only 1 firm %
    totpdiff = totpdiff /numtimes ;
    hitkmax  = hitkmax  /numtimes ;

    % totexit  = totexit/numtimes; %
    totv0  = totv0./totexit;
    totvx  = totvx./totexit;

    disp ' \n';

    printd(' RLGforced = ', padr(RLGforced,9,0), ' times, || ', RLGforced/numtimes, ' share of simulated periods (last run)');  % only uses last run since RLGforced,RLGinnov reset each run to zero %
    printd(' RLGinnov  = ', padr(RLGinnov, 9,0), ' times, || ',  RLGinnov/numtimes, ' share of simulated periods (last run)');
    disp ' \n';  % format /m1 /rds 8,2;
    % "Compare relevant surplus number to initial state value of  " v(encode(wstart),:) "  (ignore unless entry/exit shutdown so firms long-lived)"; %
    disp ' \n';
    printd('Industry characterization: ( WSTAR= ', padr(RLG_wstar,1,0), ' )');
    disp '--------------------------';
    printd('Mean disc EU period T : ', mse(discEU_T));
    printd('Mean disc CS period T : ', mse(discCS_T));
    disp ' \n';
    printd('Mean consumer surplus : ', mse(consurp));
    printd('Mean producer surplus : ', mse(prodsurp));
    printd('Mean    total surplus : ', mse(consurp + prodsurp));
    disp ' \n';
    printd('Mean CS    (adj. p*)  : ', mse(consurp2));
    printd('Mean PS    (adj. p*)  : ', mse(prodsurp2));
    printd('Mean CS+PS (adj. p*)  : ', mse(consurp2+prodsurp2));
    ' \n';
    printd('Mean industry innovate: ', mse(totinnov ));
    printd('Mean  leader  innovate: ', mse(totinnov1));
    printd('Mean  firm 2  innovate: ', mse(totinnov2));
    printd('Mean industry invest  : ', mse(totinvest));
    printd('Mean  leader  invest  : ', mse(totinves1));
    printd('Mean industry invest0 : ', mse(totinves0));
    printd('Mean  last t  invest0 : ', mse(totxT0   ));
    printd('Mean  last t front.inn: ', mse(totinnovT));
    printd('Mean freq #firms= 0   : ', mse(tot_firms(:,rlnfirms+1)));
    for nr = 1:1:rlnfirms
        printd('Mean freq #firms= ', padl(nr,2,0), ' : ', mse(tot_firms(:,nr)));  
    end
    printd('Mean leader advantage : ', mse(totw2));
    printd('Mean leader advantageT: ', mse(totw2T));
    printd('Mean # firms active   : ', mse(totnfirms));
    printd('Mean # firms invest0  : ', mse(totinve0n));
    printd('Mean # top 2 same w   : ', mse(totwsame ));
    printd('Mean # firms > wstar  : ', mse(totwstar ));
    printd('Mean # wstar /nfirms  : ', mse(totwstar./totnfirms));
    printd('Mean freq of entry    : ', mse(totentry ));
    printd('Mean freq lowrung stay: ', mse(totLRstay));
    printd('Mean freq lowrung x>0 : ', mse(totLRx   ));
    printd('Mean leader at kmax   : ', mse(hitkmax));
    printd('  max k by any firm   : ', padl(maxk,3,0));
    ' \n';
    printd('Mean share out. good  : ', mse(totshare0));
    printd('Mean share leader     : ', mse(totshare1));
    printd('Mean share firm 2     : ', mse(totshare2));
    printd('Mean share leader % T : ', mse(totshar1T));
    printd('Mean industry markup  : ', mse(totmarkup));
    printd('Mean leader   markup  : ', mse(totprice1/mc));
    printd('Mean firm 2   markup  : ', mse(totprice2/mc));
    printd('Mean leader   profit  : ', mse(totprof1));
    printd('Mean firm 2   profit  : ', mse(totprof2));
    printd('Mean price(w,:) - p*  : ', mse(totpdiff));
    ' \n';
    printd('Mean initial  value   : ', mse(totv0));
    printd('Mean realized value   : ', mse(totvx));
    printd('Mean initial  v 1st   : ', mse(totv0i));
    printd('Mean realized v 1st   : ', mse(totvxi));
    printd('  sh 1st firm exits   : ', padl(meanc(totlspi(totlspi<=numtimes)),5,3));
    disp '--------------------------';
    see_all=0;   % see_all=getbool("Would you like to see surplus realizations"); %
    if see_all
%         screen off;
%         output file = tmp.out reset;
        file = fopen('tmp.out','wt');
        fprintf(file,sprintf('%8.2f     %8.2f     %8.2f', CS, PS, TS));
        fprintf(file,sprintf('%8.2f',sort([consurp,prodsurp,(consurp+prodsurp)],3)),'descend');
        fclose(file);
%         output off;
%         screen on;
        command=[editor,' tmp.out'];
        dos command;
    end

    % save ^configfile=pp; %            % may have overwritten pp(EQL_DONE) if encountered empty industry %

    % ' '; disp "Press any key to return to the descriptive statistics menu"; ' ';  wait; %

    if BatchJob==0,  dsshell;  end % run dsshell.g

end % welf_ma


function retp = encode(ntuple)
% This function takes a weakly descending n-tuple (n = rlnfirms), with %
% min. elt. 0, max. elt. kmax, && encodes it into an integer %

  code = 1; % Coding is from 1 to wmax %
  for i = 1:rlnfirms                    % see code in eql_ma_foc.g to avoid for-loop %
    digit = ntuple(i);
    code = code + binom(digit+rlnfirms+1-i,digit+1);
  end
  retp = code;
end


function cfunk2(p)
% used for quality competition profit function %
	if RLG_y>0,  
	  n = exp(www+RLG_alpha*(ln(RLG_y-p)-ln(RLG_y)));
	else
	  % if RLG_wstar>0;  n = eg(www).*exp(-RLG_alpha*p); else ... end    <-- not needed since will never call this when RLG_wstar>0 as do in profit.g %
	  n = exp(www-RLG_alpha*p);
	end
	
	sigma = n./(RLG_out + sum(n));
	% sigma(1) = max(RLG_sh_cap|sigma(1));  <--  HOW should FOC be modified to account for SHARE CAP ? % 

	% if nfirms==2 && w(1)==10;  disp "w " w'  " prices " p'  " mktshares " sigma'  " sumshares " sum(sigma)  " exp() " n' " cfunk " (-(p-mc).*(1-sigma)./(RLG_y-p) + 1)' " profit " (M*p.*sigma-M*mc*sigma)'; 
% 	end %
      
	if RLG_y>0
        retp(-RLG_alpha*(p-mc).*(1-sigma)./(RLG_y-p) + 1);    % -RLG_alpha%(RLG_y-p) is du/dp.    dsigma/dp = sigma*(1-sigma)*du/dp  where u = utility %
    else
        % retp(-(p-mc).*sigma.*(1-sigma) + sigma);  %         % Caplin-Nalebuff FOC vector = 0 at solution.  divide by sigma to get rid of one of them  %
        retp(-RLG_alpha*(p-mc).*(1-sigma) + 1);               % the -RLG_alpha* is  du/dp (ie, utility slope) %
	end
end


function retp = newton2( p, objfunk )
% Newton-Raphson search to find the root of the function objfunk. %
% LB <= p <= UB is enforced %

  epsilon = .0001;  tol = .0001;   np = rows(p);   minstep = 1;   % stepsize is random uniform between (minstep,1) %
  trymax  = 200*np;                                               % iterations before restart search at rand() %
  
  if pp(IND_TYPE)==QUALITY
    LB = zeros(np,1) + pp(MC);
    if pp(RLG_y) >0
        UB = zeros(np,1) + pp(RLG_y) - small;    % global small defined in init.h %
    else
        UB = zeros(np,1) + maxp;                  % global maxp set as pp(RLG_MAXP)*price of maximal leader in duopoly %
    end
  else
    disp 'HEY: in Newton() you might need to add UB && LB definitions for Cournot versions of static profit game';
    UB = zeros(np,1)+99999;
    LB = zeros(np,1);
  end
  
  maxiter = trymax*5;   iter = 0;  norm = tol+1;  deriv = zeros(np,np);  x = tol+1;
  while ( (norm > tol) || (max(abs(x))>tol/100) ) && (iter < maxiter);
    iter = iter+1;
    x    = objfunk(p);    % Calculate function at p     %
    for i = 1:np       % Calculate derivative matrix %
      dp = p;
      if UB(i)-p(i) < 0.0001
          dp(i) = p(i)-epsilon;    % decrease p to get deriv instead since too near upper bound %
      else
          dp(i) = p(i)+epsilon;
      end
      deriv(:,i) = (objfunk( dp ) - x)/(dp(i)-p(i));
    end

    trap 1;    pnew = inv(deriv);  trap 0;

    if mod(iter,trymax)<4 && iter>10,
%         format /m1 /rd 7,2; 
        printd('Newton at iter ', iter, ' with w= ', wthis', ' p= ', p', ' foc= '); 
%         format /m1 /re 8,1; 
        printd( x', ' deriv= ', vecr(deriv)', ' sigma= ', sigma');    end
    
    if scalerr(pnew) || mod(iter,trymax)==0
%       format /m1 /rd 7,2;
      if scalerr(pnew)
        iter = trymax*(1+floor(iter/trymax));   % reset iter to avoid quick restart due to iter%trymax check %
        printd('restarting newton() at rand() since singular deriv at iter ', iter, ' with w= ', wthis', ' p= ', p', ' foc= ');
%         format /m1 /re 8,1; 
        printd( x', ' deriv= ', vecr(deriv)', ' sigma= ', sigma');
      else
        printd('restarting newton() at rand() since not converging at iter ', iter, ' with w= ', wthis', ' p= ', p', ' foc= ');
%         format /m1 /re 8,1; 
        disp(x');
      end
      pnew = LB+sort(rand(np,1),1,'descend');
      minstep = max([.3,minstep-.1]);          % non-converg less likely if use smaller stepsize? %
    else
      if mod(iter,trymax)<5
          pnew = p - pnew * ((.1+mod(iter,trymax)/10)*rand(1,1)) *x;     % start slow %
      else
          pnew = p - pnew *   (minstep+(1-minstep)*rand(1,1)) *x;
      end
    end
    if rows(p)>1 && sigma(1)<.001,             pnew(1) = max(pnew(1)-1|pnew(2)+.1);  end
    if rows(p)>1 && mod(iter,trymax)<5 && RLG_y>0,   pnew(1) = min([pnew(1);(RLG_y-.5+mod(iter,trymax)/10)]);  end
    pnew = min([UB';pnew']);
    pnew = max([LB';pnew']);
    norm = max(abs(pnew - p));
    p = pnew;
  end
  if norm > tol || max(abs(x))>tol/100
%       format /m1 /rds 8,5;  
      printd('newton() failed:  norm= ', norm, '  foc = ', x',  '  returning p= ', pnew') ;  
  end
  retp = pnew;
end
