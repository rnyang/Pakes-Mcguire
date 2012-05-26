function eql_ma_foc
%% eql_ma_foc.g
% LINES TO CHECK: 40, 179 , 365, 1586
%%

% This is a slightly modified version of markov.g program that allows
% for reading parameters from the configuration file.
% Version 18 uses the Allw x nfirms coding method, and allows for starting
% at any number of firms.
% Written by: Gautam Gowrisankaran
% May 16, 1993
% 
%   Modified by Ron Goettler April 2011 to solve FOC within each state each iteration.
%   Modification consists of adding a function to return foc for  1) entry prob  2) each firm's x  3) each firm's exit prob.
%   The foc are obtained simply by subtracting the candidate values from their optimal values as determined within optimize().
%   
%   To run in iterative BR mode (i.e., standard PM), the code calls  get_foc()  once, but doesn't call the Newton solver.
%   One minor difference between this iterative BR and PM is that  isentry  is updated along with other policies, not beforehand as in PM.
%   Hence, RLG added  oldentry.  The variable  isentry  is retained, instead of being renamed to  newentry.
%   If needed, the old PM timing of updating  isentry  prior to updating investments/exit could be reinstated without too much coding.
%   But the foc approach requires solving for entry simultaneously with investment/exit, so cleaner to update entry similarly when not using foc.
%   
%   When both FOC and non-FOC converge, they appear to converge to the same point.  But sometimes FOC converged when non-FOC doesn't, and vice-versa.
%   if RLG_foc = 0, use FOC method in  eql_ma_foc.g (where this value is hardcoded)  ONLY if non-FOC fails to converge.   If >0, do not try non-FOC.  If <0, never use FOC.
%

% new; %   % new; must be first line in program, else it terminates the program.  Hence, my idea to retain interactive PM code by running  new;  only when BatchJob==0 does not work %
    load('workspace');
    if BatchJob==0;  
        load pmg.mat;
    end

    x9 = 0;  y9=0;  % globals for storing best-response functions when invesigating non-convergence of nfirm = 2 case %

    % clrscr();   
    disp '**** Computing Dynamic Equilibrium using  eql_ma_foc.g  ****';
    init
    load init.mat;

%     dlibrary ./calcval ; % CHECK THIS
%     startdir = cd;
%     cd ./calcval;

    ddebug = 0;

    % constants not modifiable by user %
    tol = 0.001;        % Tolerance for convergence of value function %

    foc_tol = 1e-6;     % Tolerance for convergence of FOC solution within each state %

    rand('seed', 58921);      % random restarts in newton3() 
%     stream = RandStream.getGlobalStream;
%     reset(stream,58921);
    
    chk_foc_multi_eq = 0;  % # of random restarts of within state foc search to attempt to find a different solution %

    RLG_wstar  = pp(WSTAR);

    RLG_foc = pp(RLG_FOC);       % see documentation of RLG_foc at end of preamble %

    RLG_damp = 0; disp 'weight on old values is RLG_damp  = '; disp(RLG_damp);

    dop = 0;     % do print flag for debugging %

    if 1, RLG_w0_exit = 1; 'enforcing exit at w = 0 (as done in original PM code to signal room for entrants) ';
    else RLG_w0_exit = 0; 'NOT enforcing exit at w = 0 (as done in original PM code to signal room for entrants.  Okay if no entry/exit) ';
    end

    if phiH==0, disp ' ';  ' HEY: must revise fixed scrap to address profits earned in exit period --> nval() for exitors depends on j in optimize !! ... deliberately forcing abort with index out of range';   pp(9999999);  end

    if pp(RLG_SCRAP_HI)-pp(SCRAP_VAL) < .2,  ' HEY: solving FOC in each state may be difficult since range of scrap values is less than .2'; disp ' ';  end


    clear('newvalue','newx','newexit','oldvalue','oldx','oldexit','isentry','v_entry');


    % Set up binomial coefficients for decoding/encoding of n-tuples %

    if RLG_wstar == 0,  kmax = kmax+1;  end      % need tables 1 beyond kmax %

    binom = eye(rlnfirms+kmax+1);
    binom = [zeros(rlnfirms+kmax+1,1),binom];
    for i = 2:rlnfirms+kmax+1
      binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
    end
    binomv    = reshape(binom',size(binom,1)*size(binom,2),1);      % for RLG version of encode() that avoids for loop & %
    binomcols = size(binom',1);                                     % avoids etable1 (huge if nfirms big).  ~10% slower  %
    maxfirms_qencode = 9;                                           % maximum nfirms for which etable1 constructed for use in qencode() %
    etable1 = 0;  multfac1 = 0;                                     % initialize here in case nfirms > maxfirms_qencode  %

    if RLG_wstar == 0,  kmax = kmax-1;  end      % RESTORE kmax %


    % RLG asks what is value of encfirm < rlnfirms ??  The execution time is same.  And eql_sa.g crashes if encfirm < rlnfirms %
    encfirm = rlnfirms;      % was hardcoded to something like 4,5, or 6.   Max. number of firms to encode in table %

    filenam1 = [prefix , 'markov.dp'];
    printd('Execution protocol --> ', filenam1);
%     output file = ^filenam1 reset;
    file = fopen(filenam1,'wt');
%     fclose(file);
%     file = fopen(filenam1,'at');
    nfirms = stfirm;        % max. # of firms at initial computation stage %
    nfirms = rlnfirms;      % hardcoded ignoring of option to use equilibrium from industry with fewer firms as starting point %
    oneton = 1:1:nfirms;
    oneto2nfirms1 = 1:1:2*nfirms+1;

    LB = zeros(1+2*nfirms,1);  % bounds on policies:  x in 1:nfirms, exit in nfirms+1:2*nfirms, entry in 1+2*nfirms %
    UB =  ones(1+2*nfirms,1);  UB(1:nfirms) = UB(1:nfirms) + 9999999999;  % investment has no upper-bound %


    % Next blocks are somewhat messy since PM allowed for possibility of using equilibrium with fewer firms as starting values  %
    % RLG ignores this option and also sets encfirm = nfirms.  Hence, multfac2 is not used and the opening block is skipped.    %
    % Much of code in next few blocks could therefore be removed, but is retained in case desire to use it again someday.       %
    % RLG also added saving and loading of  dtable  etable1  mask, since building them is SLOW (nested loops) when many firms.  %

    if   0   && nfirms > 1   % Skip this block when don't want to use model with fewer firms as starting values %
      nfirms = nfirms - 1;
      wmax = binom(nfirms+1+kmax,kmax+2);
      if 0
        % Read data: v (value), x (investment), p (probability of state rising), isentry, newexit %
        filename = [prefix , 'markov.' , padr(nfirms,1,0) , 'ot'];
        bigread = load(filename);
        newvalue = (reshape(bigread(1:wmax*nfirms),wmax,nfirms));
        bigread = bigread(wmax*nfirms+1:size(bigread,1));
        newx = (reshape(bigread(1:wmax*nfirms),wmax,nfirms));
        if phiH ~= 0
          bigread = bigread(wmax*nfirms+1:size(bigread,1));
          newexit = (reshape(bigread(1:wmax*nfirms),wmax,nfirms));
        end
        clear bigread;               % isentry not read in since it's computed first %
      else
        newx     = zeros(wmax,nfirms);
        newexit  = zeros(wmax,nfirms);
        newvalue = zeros(wmax,nfirms);
      end
      oneton =  1:1:nfirms;

      if nfirms >= encfirm
        multfac1 = (kmax+1)^(oneton(1:encfirm)-1);
        nfirms = encfirm;
        % Encode all numbers from 1 to kmax^nfirms %
        etable1 = zeros((kmax+1)^nfirms,1);
        for i = 1:1:size(etable1,1)     % was:  i=0;  while i < length(etable1); ... i= i+1; end %
          msk = []; % Build one mask %
          j=kmax+1; k=i-1;
          while j <= size(etable1,1)
            msk = [(mod(k , j)*(kmax+1)/j);msk];
            k = k - mod(k , j);
            j = j*(kmax+1);
          end
          temp = sort(msk(1:nfirms),1);
          etable1(i) = encode2(temp(end:-1:1));
        end
        nfirms = stfirm-1;
        if nfirms > encfirm;
          multfac1 = [zeros(nfirms-encfirm,1);multfac1];
        end
      end
      nfirms = stfirm;
    end

    while nfirms <= rlnfirms
      nfirms_oneton = nfirms-1:1:nfirms;     % used in encode3() and qencode() %

      filenam2 = [prefix , 'markov.' , num2str(nfirms) , 'ot'];
      d1 = date;
    %   format /rd 9,6;
      % Number of different combinations of competitors positions faced by a firm %
      wmax = binom(nfirms+1+kmax,kmax+2);
      fprintf(file,'Number of firms: %s; Total number of states: %s\n', compact(nfirms),compact(wmax));
      fprintf(file,'Initializing ...');
      filename = [prefix , 'pr.' , num2str(nfirms) , 'f'];

      i = wmax;
      if RLG_wstar==0,   i = binom(nfirms+2+kmax,kmax+3);  end    % RLG_wstar==0 extended kmax by one in profit.g %
      
%       profit(i,nfirms)=load(filename);
      profit(i,1:length(load(filename)))=load(filename)';

      if nfirms<maxfirms_qencode
        % check if can load  dtable etable1 and mask  since they take LONG time to build if nfirms high %
        filename = ['bmarkov.dtable.' , num2str(nfirms) , 'f.' , num2str(kmax) , 'k'];
%         if (not filesa(filename) $== "");
        filechar = ls;
        filesa = {};
        for i=1:size(filechar,1),filesa{end+1}=filechar(i,:);end
        if any(strcmp(filesa,filename)) % CHANGE TO FILE SEARCH IN DIRECTORY 
          oneton = 1:1:nfirms;          % these variables also created in the code below that builds  dtable etable1 and mask %
          multfac1 = (kmax+1)^(oneton-1);
          two_n = 2^(nfirms-1);
          
          dtable(nfirms,wmax)=load(filename);

          filename = ['bmarkov.etable1.' , num2str(nfirms) , 'f.' , num2str(kmax) , 'k'];
          etable1((kmax+1)^nfirms,1) = load(filename);

          filename = ['bmarkov.mask.' , num2str(nfirms) , 'f.' , num2str(kmax) , 'k'];
          mask(nfirms-1,two_n) = load(filename);

          fprintf(file,'Loaded dtable etable1 and mask instead of building them...\n');
        else   % rebuild dtable etable1 and mask %

          two_n = 2^(nfirms-1);
          clear dtable;
         if nfirms > 1
            oneton = 1:1:nfirms;
            % Build a mask of all binary numbers from 0 to two_n - 1 %
            mask = zeros(nfirms-1,two_n);
            for i = 1:two_n                 % was:  i=0;  while i < two_n; ... i= i+1; end %
                msk = []; % Build one mask %
                j=2; k=i-1;
                while j <= two_n
                  if mod(k , j) == 0
                      msk = [0;msk];
                  else
                      k = k - j/2; msk = [1;msk];
                  end
                  j = j*2;
                end
                mask(:,i) = msk(1:nfirms-1);
            end
        disp 'Mask done';
        end

          % Make a table for quick decoding %
          dtable = zeros(nfirms,wmax);
          for i = 1:1:wmax,   dtable(:,i) = decode(i);   end

          % Make a table for quick encoding %
          % Fill in multfac1, multfac2, for quick encoding %
          if nfirms <= encfirm
            multfac1 = (kmax+1)^(oneton-1);
            % Encode all numbers from 1 to kmax^nfirms %
            etable1 = zeros((kmax+1)^nfirms,1);
            for i = 1:1:size(etable1,1)     % was:  i=0;  while i < length(etable1); ... i= i+1; end %
                msk = []; % Build one mask %
                % j=kmax+1; %
                k = i-1;
                if mod(i,100000)==0,  printd('etable1 row ', i, ' of ', size(etable1,1));  end
                  % while j <= length(etable1); %       % loop will execute nfirms times %
                  for t0 = 1:1:nfirms
                    j = (kmax+1)^t0;
                    msk = [(mod(k , j)*(kmax+1)/j);msk];
                    k = k - mod(k , j);
                    % j = j*(kmax+1);  end %
                  end
                  temp = sort(msk(1:nfirms),1);
                  save('workspace');
                  etable1(i) = encode2(temp(end:-1:1));
            end
            disp 'etable1 done';
            if nfirms>1     % save file since creating dtable etable1 and mask are time-consuming %
              % write output %
              disp 'Saving dtable  etable1  and  mask ';
        % 	  format /m1 /rds 1,0;
    %           screen off;
              fclose(file);
              filename = ['bmarkov.dtable.' , num2str(nfirms) , 'f.' , num2str(kmax) , 'k'];
        % 	  output file = ^filename reset;
              file = fopen(filename,'wt');
              fprintf(file,num2str(dtable));
              fclose(file);

              filename = ['bmarkov.etable1.' , num2str(nfirms) , 'f.' , num2str(kmax) , 'k'];
        % 	  output file = ^filename reset;
              file = fopen(filename,'wt');
              fprintf(file,num2str(etable1));
              fclose(file);

              filename = ['bmarkov.mask.' , num2str(nfirms) , 'f.' , num2str(kmax) , 'k'];
        % 	  output file = ^filename reset;
              file = fopen(filename,'wt');
              fprintf(file,num2str(mask));
              fclose(file);

    %           output file = ^filenam1;
              file = fopen(filenam1,'at');
    %           screen on;
            end  % save  dtable etable1 and mask  files %

          else        % nfirms > encfirm.   NO LONGER EVER HERE since RLG set encfirm = nfirms %
            multfac1 = [0;multfac1];
            multfac2 = (kmax+1)^[(oneton(1:nfirms-encfirm)-1);zeros(encfirm,1)];
            %
            disp 'Multfac1 is '; disp(multfac1);
            disp 'Multfac2 is '; disp(multfac2);
            %
            etable2 = zeros((kmax+1)^(nfirms-encfirm),1);
            for i = 1:1:size(etable2,1)     % was:  i=0;  while i < length(etable2); ... i= i+1; end %
              msk = []; % Build one mask %
              j=kmax+1; k=i-1;
              while j <= size(etable2,1)
                msk = [(mod(k , j)*(kmax+1)/j);msk];
                k = k - mod(k , j);
                j = j*(kmax+1);
              end
              temp = sort(([msk(1:nfirms-encfirm);zeros(encfirm,1)]),1);
              etable2(i) = encode2(temp(end:-1:1))-1;
            end
            disp 'etable2 done';
          end
        end    % load or rebuild  dtable etable1 and mask %
      else       % nfirms > maxfirms_qencode  --> build mask only %
        two_n = 2^(nfirms-1);
        clear dtable;
        oneton = 1:1:nfirms;
        % Build a mask of all binary numbers from 0 to two_n - 1 %
        mask = zeros(nfirms-1,two_n);
        for i = 1:1:two_n                 % was:  i=0;  while i < two_n; ... i= i+1; end %
        msk = []; % Build one mask %
        j=2; k=i-1;
        while j <= two_n
          if mod(k , j) == 0
        msk = [0;msk];
          else k = k - j/2; msk = [1;msk];
          end
          j = j*2;
        end
        mask(:,i) = msk(1:nfirms-1);
        end
        disp 'Mask done.  Not building etable1 for qencode since nfirms > maxnfirms_qencode';
        % Make a table for quick decoding %
        dtable = zeros(nfirms,wmax);
        for i = 1:1:wmax,   dtable(:,i) = decode(i);   end
      end

%       InitializeZeros:
%       InitializeZeros = 1;
%       while InitializeZeros == 1
%           InitializeZeros = 0;
       
      % Update values, or define starting values.  RLG moved next block from update() to here %
      oldx     = zeros(wmax,nfirms); 
      oldexit  =  ones(wmax,nfirms);    % oldexit init as 1's since 1 is the newexit for inactive firms %
      oldvalue = zeros(wmax,nfirms); 
      oldentry = zeros(wmax,1);      
      v_entry  = zeros(wmax,1);

      if size(RLG_x,1)==size(oldx,1)
        fprintf(file, 'initializing with RLG_* values from previous model...\n');
        oldx     =  RLG_x     ;
        oldexit  =  RLG_exit  ;
        oldvalue =  RLG_value ;
        oldentry =  RLG_entry ;
        v_entry  =  RLG_v_entry  ;
      end

      foc          = zeros(wmax,2*nfirms+1);   % these 4 globals simplify multithreading with foc %
      rflag        = zeros(wmax,1);
      x_exit_entry = zeros(wmax,2*nfirms+1);
      foc_multi_eq = zeros(wmax,1);            % # of iterations for which contract_w() finds a 2nd equilibrium, in the chk_foc_multi_eq random restarts to check for multi_eq %

      if delta==0 && RLG_wstar<=0          % Be sure simulations start with a firm at kmax %
        temp = 1:wmax;
        do_w = temp(dtable(1,:)'<kmax);	% then will ALWAYS have a firm at kmax, so skip states with 1st firm below kmax % 
        oldx(do_w,:)     = oldx(do_w,:)-1;
        oldexit(do_w,:)  = oldexit(do_w,:)-1;
        oldvalue(do_w,:) = oldvalue(do_w,:)-1;
        oldentry(do_w,:) = oldentry(do_w,:)-1;
        do_w = temp(dtable(1,:)'==kmax);
        fprintf(file, 'Only solving for values and policies over states with firm 1 at the frontier, since delta = 0\n');
      else
        do_w = 1:1:wmax;
      end
      wmax2 = size(do_w,1);
      %% THE LINE BELOW has been modified to reflect random order
      do_w = do_w(round(rand(wmax2,1)*(wmax2-1)+1));
%       do_w = sort( [do_w,rand(wmax2,1)], 2);   % randomize the order states are to be functionessed  %
      do_w = do_w(:,1);                       % so threads have similarly time-consuming states %

      newx     = oldx;      % initialize new here, after set -1 for ignored states %
      newexit  = oldexit;
      newvalue = oldvalue;
      isentry  = oldentry;

      % update(); %         % skipping update() which uses values for nfirms-1 equilibrium as starting values %

      maxiter = 400/(1-RLG_damp);

%       StartIterating:
      StartIterating = 1;
      while StartIterating == 1
          StartIterating = 0;

    %   format /m1 /re;
      fprintf(file, 'Contraction...over %s  states of %s in dtable with RLG_foc = %s  min frontier firm innov = %f  maxiter= %s  tol= %f  foc_tol= %s ...\n',padr(wmax2,1,0),padr(wmax,1,0),padr(RLG_foc,1,0),aeff(size(aeff,1))*RLG_minx1/(1+aeff(size(aeff,1))*RLG_minx1),padr(maxiter,1,0),tol,foc_tol);
      if chk_foc_multi_eq,  disp 'In each iteration and at each state, restarting search for FOC solution at ';  disp(padr(chk_foc_multi_eq,1,0)); disp ' random policies to investigate presence of multiple equilibria.';  end
    %   format /m1 /rd 8,4;
      norm = tol + 99;  iRLG=1;   normx = norm;  normexit= norm;  normentry = norm;
      t = clock;
      t = round(sum(t(end-2:end).*[60*60,60,1]));
      avgnorm = norm;  prevnorm = norm+1;  prevavgnorm = norm+1;  t0 = t;   % hsec= hundredths of seconds since midnight %
      norm100 = norm+9;

      while (norm > tol) && (iRLG < maxiter || (norm-prevnorm<.01 && iRLG<400) || (avgnorm-prevavgnorm<1e-7 && iRLG<400));  %  && (avgnorm > 0.0001*tol);   <--  REMOVED by RLG the avgnorm criteria %    
        prevnorm = norm;

        if mod(iRLG,100)==0;  norm100= (norm+prevnorm)/2;     % attempt to detect non-converg to avoid wasted effort.  avg avoids getting peak norm of a cycle. %
        else
          if iRLG>100 && iRLG-100*floor(iRLG/100)>31 && norm>norm100;           % if iRLG > 100 and higher norm than 31 ago, likely not converging %
        fprintf(file,'aborting contraction since current norm is worse than prev 100th iteration, suggesting convergence issues, ...');
        break;
          end
        end

        %
        if max(max(oldexit))>1
          xxx = indexcat(max(oldexit')>1,1);
          printd(size(xxx,1), ' states have oldexit>1  at iteration ', iRLG+1);
          disp 'Element no./Probability Exit (wmax(nfirms-1)): '; disp([dtable(:,xxx)',oldexit(xxx,:)]); disp ' ';
        end

        xxx = sum( (oldexit(:,1:nfirms-1) > oldexit(:,2:nfirms) + .01)' &   (oldvalue(:,1:nfirms-1)-profit(:,1:nfirms-1) > oldvalue(:,2:nfirms)-profit(:,2:nfirms) + .01)' )>0 ;
        if sum(xxx)>0             % -profit() in prev line since exit is based on continuation values (ie, profits are earned regardless of entry) %
          aaaaa = indexcat(xxx,1);
    %       format /m1 /rd 6,2;
          printd( 'HEY: exit goes down when value goes up as w increases for ', size(aaaaa,1), ' states:   w ~ value ~ exit ~ x');
          disp([dtable(:,aaaaa)',-9*ones(size(aaaaa,1),1),oldvalue(aaaaa,:),-9*ones(length(aaaaa),1),oldexit(aaaaa,:),-9*ones(length(aaaaa),1),oldx(aaaaa,:)]);
        end
        %

        if BatchJob== -1 && iRLG>200,  RLG_damp = .99;  maxiter = 200000; end
        save workspace
        contract();

        abs_error = abs(oldvalue - newvalue);
        norm      = max(max(abs_error));
        avgnorm   = mean(mean(abs_error));
        normx     = max(max(oldx-newx));
        normexit  = max(max(oldexit-newexit));
        normentry = max(oldentry-isentry);
       % 'Sup norm: ' padr(norm,8,4) '; Mean norm: ' padr(avgnorm,8,4) '    \r'; %

        iRLG = iRLG+1;
        'TEST',iRLG,nfirms
        if mod(iRLG ,10)==0 || nfirms>9 
          fprintf(file,'iteration: %s avg~sec %s  norm~avgnorm: %s ~ %s  norm x~exit~entry: %s   %s   %s   maxV %s   ',padr(iRLG,6,0),padr((hsec-t0)/100/(1+9*(nfirms<=9)),6,1),padr(norm,8,5),padr(avgnorm,11,8),padr(normx,10,6),padr(normexit,10,6),padr(normentry,10,6),padr(max(newvalue(:,1)),7,2)) ; % print /flush;
          t0 = hsec;
        end

        % if iRLG<4;  oldvalue~newvalue;   end %

        if 0 && iRLG > 100
    %       format /m1 /rd 8,2;  
          xxx = ones(length(isentry),1)-1.11;
          disp 'dtable ~ values (old then new then diff)';  disp([dtable',oldvalue,newvalue,xxx,newvalue-oldvalue]);
          disp 'dtable ~ x (old then new then diff)'; disp([dtable',oldx,newx,xxx,newx-oldx]);
          disp 'dtable ~ exit ~ entry (old then new then diff)'; disp([dtable',oldexit,newexit,xxx,newexit-oldexit,xxx,oldentry,isentry,xxx,isentry-oldentry]);
        end

        if nfirms<7 && norm > .001 && iRLG>100;

%           screen off;

          normind1 = maxindc(max(abs_error'));    % oldvalue is wmax by nfirms.  so max(abs_error') is wmax vector and normind1 is state with highest error %
          normind2 = maxindc(max(abs_error));     % max(abs_error) is firm vector of highest error across states, so normind2 is firm with highest error %

          if BatchJob==0   && 0  ,  normind1 = qencode([7;7;5]);  normind2 = 1;  end

          normcode = dtable(:,normind1)';
    %       format /rd 7,4;

          aaaaa = abs( newx-oldx );           xind1 = maxindc(max(aaaaa'));   xcode = dtable(:,xind1)';
          bbbbb = abs( newexit-oldexit );     eind1 = maxindc(max(bbbbb'));   ecode = dtable(:,eind1)';
          %
          if sum(sum(aaaaa))==1
        RLG1 = maxindc(max(aaaaa'));
        RLG2 = maxindc(max(aaaaa));
          else
        RLG1 = 1; RLG2 = 0;
          end
          %
          % 'length~cols oldvalue:'  length(oldvalue)~cols(oldvalue) 'length~cols oldx:'  length(oldx)~cols(oldx); %

          fprintf(file,'Norm~Avgnorm: %s   %s   ',padr(norm,8,5),padr(avgnorm,11,8));
          fprintf(file,'Max Norm: firm %s at state %8.4f',padr(normind2,1,0),normcode);
          fprintf(file,' Old~New value: %8.4f ~ %8.4f   Old~New x: %8.4f ~ %8.4f   Old~New exit: %8.4f',oldvalue(normind1,:),newvalue(normind1,:),oldx(normind1,:),newx(normind1,:),oldexit(normind1,:),newexit(normind1,:));
          fprintf(file,'   Entry: %8.4f %8.4f   Old~New x at %8.4f is %8.4f~%8.4f',isentry(normind1), v_entry(normind1),xcode,oldx(xind1,:),newx(xind1,:));

          fprintf(file,'   Old~New exit at %8.4f is %8.4~%8.4',ecode,oldexit(eind1,:),newexit(eind1,:));
          % ' sole exit change at firm ' padr(RLG2,1,0) ' state ' dtable(:,RLG1)'  %  ;
          fprintf(file,'\n');
%           screen on;
        end

        if RLG_damp || (  0 &&  iRLG > 100 && norm>prevnorm) 
          % i = RLG_damp/max(iRLG|(iRLG/100));  %              % can decrease dampening as iRLG increases %
          i = max([RLG_damp;(.5*(norm>prevnorm))*rand(1,1)]);    % random weight on oldx, with max at RLG_damp or .5 if norm>prevnorm %
          if norm>prevnorm && RLG_damp==0,  disp 'norm > prevnorm, dampening weight on old = '; disp(i); disp ' at iter= '; disp(iRLG);  end
          oldx     =     oldx*i + (1-i)*newx;
          oldexit  =  oldexit*i + (1-i)*newexit;
          oldvalue = oldvalue*i + (1-i)*newvalue;
          oldentry = oldentry*i + (1-i)*isentry;
        else
          oldx = newx;  oldvalue = newvalue;  oldexit = newexit;   oldentry = isentry;
        end

    %     format /rd 9,6;   
%         output off;   output on;

        %
%         i = seqa( length(oldvalue)-9,1,10);
%         'iter = '  
%         iRLG  
%         '  qencode_index  value   x   exit '
    %     format /m1 /rds 8,3;  i~oldvalue[i,.]~oldx[i,.]~oldexit[i,.];
        %

        % if RLG_foc==0 AND iRLG > 50; disp ' ';  'Experimental hardcoded break of non-foc effort to switch to foc ';  disp ' ';    break;  end %

      end  % while norm > tol %
      end % StartIterating


    %   format /rd 7,3;
      d2 = date;
      ethsec = (datenum(d2) - datenum(d1))*24*60*60;
      fprintf(file,'%7.3f firm(s) scenario completed; execution time = ',compact(nfirms));
      fprintf(file,'%s minutes, #iterations= %s',padr(ethsec/6000,8,2), padr(iRLG,6,0) ); 
    %   format /m1 /re 6,1; 
      fprintf(file,'tol= %6.1f  norm=%6.1f  prevnorm~avgnorm=%6.1f %6.1f  norm x~exit~entry= %6.1f %6.1f %6.1f \n',tol,norm,prevnorm,prevavgnorm,normx,normexit,normentry);

      % isentry(99999999); %    % deliberate abort %


      if max( foc_multi_eq ) > 0
        temp = 1:wmax;
        d2 = temp(foc_multi_eq > 0);          % states with multiple eq %
        if length(d2)>30
          d2 = sort(d2,1);
          d2 = d2(end:-1:1);
          d2 = d2(1:30);
        end
        kk = 1.111 + zeros(length(d2),1);
    %     format /m1 /rd 4,0;  
        fprintf(file,'Multiple FOC solutions were obtained most frequently at the following states:   # multieq  ~  state  ~  entry  ~  investment  ~  exit');
    %     format /m1 /rd 8,4;  
        fprintf(file,'%8.4f   %8.4f   %8.4f   %8.4f   %8.4f   %8.4f   %8.4f   %8.4f   %8.4f   ',foc_multi_eq(d2),kk,dtable(d2,:),kk,isentry(d2),kk,newx(d2,:),kk,newexit(d2,:));
      end

      if norm>tol,   fprintf(file,' ');     % if BatchJob and RLG_minx1>0, do NOT try FOC -- unlikely to help since already tried FOC with minx1=0 and non-FOC with minx1>0 %
        if RLG_foc==0 && pp(RLG_FOC)==0 && (BatchJob==0 || RLG_minx1==0); 
          RLG_foc = 1;
          maxiter = 200;  fprintf(file,'non-FOC failed to converge... retrying from current with FOC method...');   StartIterating = 1; break;%goto StartIterating; 
        end
        if  0 &&  BatchJob>0 && RLG_minx1==0                  % since this overwrites the config pp() only do this when in Batch mode %
          fprintf(file,'non-FOC and possibly FOC failed to converge ... retrying non-FOC from current after setting pp(RLG_MIN_INNOV1)= .01...   ');
          pp(RLG_MIN_INNOV1)= .01;
          RLG_minx1 = pp(RLG_MIN_INNOV1)/((1-pp(RLG_MIN_INNOV1))*aeff(rlnfirms+1));  % ax = (1+ax) m  --> ax(1-m) = m -->  x = m/a(1-m) % 
          RLG_foc = 0;  maxiter = 200;
%           goto StartIterating;
          StartIterating = 1; break;
        end

        if RLG_foc >0,    'FOC method failed to converge...';      end
        if RLG_foc==0,    'non-FOC method failed to converge...';  end
      end

      % Now find if there is any investment, entry at the highest level %

      probw = reshape(aeff(dtable'+1),wmax,nfirms).*newx;     probw = probw./(1+probw);

      d2 = ' ';
      for kk = nfirms:1:9,  d2 =  ['0 ' , d2];  end   
      d2 = ['  ' , d2];   % build string of 0 for filler columns in output below %
      fprintf(file,'output below corresponds to 10 firm case since matlab code plotting policy functions across models uses column counts of 10 firm industry   \n');
      fprintf(file,' V and policies at select states:     w   ~  V  ~  entry  ~  innovation  ~  exit   \n');
      w = ones(nfirms,1);
    %   format /m1 /rd 3,0;  
      fprintf(file,'%3.0f   ',w'); 
    %   format /m1 /rd 8,4; 
      printd(d2, newvalue(qencode(w),:), d2,  isentry(qencode(w)));
      printd(probw(qencode(w),:),  d2, newexit(qencode(w),:), d2);
      w
      qencode(w)
      newvalue
      newvalue(qencode(w),:)
      isentry(qencode(w))
      fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2, newvalue(qencode(w),:), d2,  isentry(qencode(w)));
      fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),  d2, newexit(qencode(w),:), d2);
      fprintf(file,' ');
      w = kmax*ones(nfirms,1);      
    %   format /m1 /rd 3,0;  
      fprintf(file,'%3.0f   ',w');  
      fprintf(file,'\n');
    %   format /m1 /rd 8,4; 
      w
      qencode(w)
      newvalue
      newvalue(qencode(w),:)
      isentry(qencode(w))
      fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
      fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
      for kk = nfirms:-1:2
        w(kk)=0;    % format /m1 /rd 3,0;  
        fprintf(file,'%3.0f   ',w'); 
    %    format /m1 /rd 8,4; 
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w))); 
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
      end
      w = entry_k*ones(nfirms,1);  wretp = kmax;    % format /m1 /rd 3,0;  
      fprintf(file,'%3.0f   ',w'); 
    %    format /m1 /rd 8,4; 
      fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
      fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
      for kk = nfirms:-1:2
        w(kk)=0;   
    %     format /m1 /rd 3,0;  
        fprintf(file,'%3.0f   ',w'); 
    %     format /m1 /rd 8,4; 
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
      end
      w = zeros(nfirms,1);  wretp = kmax;
      for kk = kmax:-1:1
        w(2)=kk;   
    %     format /m1 /rd 3,0;  
        fprintf(file,'%3.0f   ',w'); 
    %     format /m1 /rd 8,4; 
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
      end
      w = ones(nfirms,1);  wretp = kmax;  w(nfirms) = 0;
      for kk = kmax:-1:1
        w(2)=kk;   % format /m1 /rd 3,0;  
        fprintf(file,'%3.0f   ',w');  % format /m1 /rd 8,4; 
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
      end
      w = ones(nfirms,1);  wretp = kmax;
      for kk = kmax:-1:1;
        w(2)=kk;   
    %     format /m1 /rd 3,0;  
        fprintf(file,'%3.0f   ',w');  
    %     format /m1 /rd 8,4; 
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
        fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
      end
      if nfirms>2
        w = zeros(nfirms,1);  wretp = kmax;  w(2) = kmax;
        for kk = kmax:-1:1
          w(3)=kk;   
    %       format /m1 /rd 3,0;  
          fprintf(file,'%3.0f   ',w');  
    %       format /m1 /rd 8,4; 
          fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
          fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
        end
      end
      disp ' ';
      if nfirms>2
        w = zeros(nfirms,1);  wretp = kmax;  w(2) = kmax-1;
        for kk = kmax-1:-1:1
          w(3)=kk;    
    %       format /m1 /rd 3,0;  
          fprintf(file,'%3.0f   ',w');  
    %       format /m1 /rd 8,4; 
          fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
          fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);
        end
      end
      disp ' ';
      w = nfirms:-1:1;  w = min([w';kmax*ones(1,nfirms)]);
      w(1)=kmax;  
    %   format /m1 /rd 3,0;  
      fprintf(file,'%3.0f   ',w'); 
    %   format /m1 /rd 8,4; 
      fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',d2,newvalue(qencode(w),:),d2,isentry(qencode(w)));
      fprintf(file,'%8.4f %8.4f %8.4f %8.4f   ',probw(qencode(w),:),d2,newexit(qencode(w),:),d2);

      w= sum(newexit'>0)==nfirms;
      if sum(w)>1
        temp = 1:wmax;
        iRLG = temp(w);
        fprintf(file,'All firms exit at the following states: w ~ V ~ x ~ entry ~ exit');
    %     format /m1 /rds 5,2; 
        fprintf(file,'%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f   \n',dtable(:,iRLG)',newvalue(iRLG,:),-9.99*ones(length(iRLG),1),newx(iRLG,:),-9.99*ones(length(iRLG),1),isentry(iRLG),-9.99*ones(length(iRLG),1),newexit(iRLG,:));
      end

      if RLG_wstar>0
        w = zeros(nfirms,1);  wretp = kmax-1;
        kk = newx(qencode(w):(qencode([kmax;zeros(nfirms-1,1)])-1),1); % firm 1 investment over all states with firm 1 at kmax-1 %
        iRLG = max(kk);
        if iRLG > 0
          fprintf(file,'Warning: Positive investment of  %f  recorded at 2nd highest efficiency level, yielding success prob = %f',iRLG,aeff(kmax)*iRLG/(1+aeff(kmax)*iRLG));
          fprintf(file,'Please consider increasing the maximum efficiency level (kmax) since using PM wstar to bound statespace.');
        end
      else
        w = zeros(nfirms,1);  wretp = kmax;
        kk = newx(qencode(w):wmax,1); % firm 1 investment over all states with firm 1 at frontier %
        iRLG = max(kk);
        fprinf(file,'Max investment of  %f  recorded at highest efficiency level, yielding success prob = %f (which is okay since not using wstar to bound statespace)', iRLG,  aeff(kmax+1)*iRLG/(1+aeff(kmax+1)*iRLG));
        iRLG = kk-RLG_minx1<1e-15;
        if sum(iRLG)>0
          temp = qencode(w):1:length(kk);
          iRLG = temp(iRLG);
          fprintf(file,'Minimum investment of %f by frontier firm at %f of %f states with leader at frontier:  w ~ V ~ x ~ entry ~ exit', RLG_minx1, length(iRLG), length(kk));
    %       format /m1 /rds 6,2;  
          fprintf(file,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n',dtable(:,iRLG)',newvalue(iRLG,:),-9.99*ones(length(iRLG),1),newx(iRLG,:),-9.99*ones(length(iRLG),1),isentry(iRLG),-9.99*ones(length(iRLG),1),newexit(iRLG,:));
        end
      end
      disp ' ';
      w = sum( (dtable(:,do_w)==1 & newexit(do_w,:)'<1-1e-8));      % states with firm not exiting for sure at lowest rung %
      if sum(w)>0
        kk = do_w(w>0);
        fprintf(file,'Firms on lowest rung w/ Pr(exit)<1 at %s states.  They get free bump if RLG_w0_exit=0, else exit w/ avg scrap when fall off ladder.  Hardcoded RLG_w0_exit is %s', padr(length(kk),1,0), padr(RLG_w0_exit,1,0));
        fprintf(file,'First few of such states:  value  ~  invest  ~  entry  ~  exit');
        kk = kk(1:min([10;length(kk)]));  w = -9.99*ones(length(kk),1);
    %     format /m1 /rds 6,2;  
        fprintf(file,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f	\n',dtable(:,kk)',newvalue(kk,:),w,newx(kk,:),w,isentry(kk),w,newexit(kk,:));
      end
      w = sum( (dtable(:,do_w)==1 & newexit(do_w,:)'<1-1e-8 & newx(do_w,:)'>1e-8));      % indicator for states with firm at lowest rung with Pr(exit)<1 and investing x>0 %
      if sum(w)>0
        kk = do_w(w>0);
        fprintf(file,'Firms on lowest rung w/  invest> 0 at %s states.  They get free bump if RLG_w0_exit=0, else exit w/ avg scrap when fall off ladder.  Hardcoded RLG_w0_exit is %s', padr(length(kk),1,0), padr(RLG_w0_exit,1,0));
        fprintf(file,'First few of such states:  value  ~  invest  ~  entry  ~  exit');
        kk = kk(1:min([10;length(kk)]));  w = -9.99*ones(length(kk),1);
    %     format /m1 /rds 6,2;  
        fprintf(file,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f	\n',dtable(:,kk)',newvalue(kk,:),w,newx(kk,:),w,isentry(kk),w,newexit(kk,:));
      end
      disp ' ';

      % Store data in file for inspection %
%       screen off;
      fprintf(file,'Element no./ Value Function (wmax x nfirms): %f %f   \n', dtable',newvalue);
      fprintf(file,'Element no./Investment (wmax x nfirms): %f %f   \n', dtable',newx);
      fprintf(file,'Element no./Probabilities of p rising (wmax x nfirms) ~ Investment ~ Efficiency: %f %f %f %f	\n',dtable,probw,newx,reshape(aeff(dtable'+1),wmax,nfirms));
      if phiH ~= 0
        fprintf(file, 'Element no./Probability Exit (wmax(nfirms-1)): %f %f   \n', dtable',newexit);
      end
      fprintf(file,'Element no./Probability ~ Value of entry (wmax(nfirms-1)): %f %f %f   \n', dtable',isentry,v_entry);
      temp = [dtable',isentry,v_entry];
      fprintf(file,'Element no./Probability ~ Value of entry (wmax(nfirms-1)): only states with room for entrant %f   \n', temp(dtable(nfirms,:)'==0));
%       screen on;
%       output off;
      fclose(file);
      % Store data in file for comparative statics program to read %
      printd( 'Generating data for comparative statics --> ', filenam2);
    %   format /m1 /rds 16,12;
%       output file = ^filenam2 reset;
      file = fopen(strcat(filenam2,'.mat'),'wt');
%       screen off;
      fprintf(file,sprintf('%16.12f %16.12f %16.12 ',newvalue, newx, probw));
      if phiH ~= 0
        fprintf(file,sprintf('%16.12f ',newexit));
      end
      fprintf(file,sprintf('%16.12f %16.12f ',isentry,v_entry));
%       output off;
      fclose(file);
%       screen on;
%       output file = ^filenam1 on;
      file = fopen(strcat(filenam1,'.mat'),'wt');
      nfirms = nfirms+1;

      % XP = (w./(1+w));  output file=XX ;  output on;  print  (kmax*ones(length(newx),1))~dtable'~newx~XP~XP./max(XP);  output off; %

    %   format /m1 /rd 5,2;
      if delta==0 && RLG_wstar<=0
        iRLG = sum(newx')-RLG_minx1*sum(dtable==kmax)<1e-15 & isentry==0 & sum(newexit')==0;
        % if sum(iRLG)>0;  disp ' '; 'absorbing states with MIN INVESTMENT, NO ENTRY, NO EXIT: state ~ value ~ exit ~ entry';  selif( dtable'~newvalue~newexit~isentry, iRLG ); disp ' ';  %
        if sum(iRLG)>0,  fprintf(file,sprintf(' absorbing states with MIN INVESTMENT, NO ENTRY, NO EXIT:%5.2f ',  dtable( iRLG' )));  
        end
      end

    end   % while nfirms <= rlnfirms.   nfirms started at stfirm which itself is set to  stfirm = pp(START_FIRMS);  in init.h  %

    if norm <= tol             % save in global variables to use as initial values for next model %
      RLG_x        = oldx    ;
      RLG_exit     = oldexit ;
      RLG_value    = oldvalue;
      RLG_entry    = oldentry;
      RLG_v_entry  = v_entry ;
    end

%     output off;
    fclose(file);
%     screen on;
    pp(EQL_DONE)= (norm <= tol);

    % print 'Press any key to return to the main menu'; wait; %

    % save ^configfile=pp;
    save(configfile,pp);

    if BatchJob==0,  pmgshell;  end


end % eql_ma_foc

function retp = get_foc(x_exit_ent, w, doo)
  % obtain foc for state w at candidate policy x_exit_ent  %
  % x_exit_ent   are policies being changed to find root   %
  % doo          are indices of policies being changed     %
  % using parameters, not globals, since multithreaded     %
  % requires each thread to have its own copy.             %
  %                                                        %
  % modified optimize() and chkentry() to accept candidate policies %
  % as inputs instead of using global oldx oldexit isentry %
  load workspace
  x_exit_entry_w = ( [oldx(w,:),oldexit(w,:),oldentry(w) ])';
	
  if length(doo)>0,   x_exit_entry_w(doo) = x_exit_ent;   end   % possible doo = {} %
  locw  = dtable(:,w)';
  for i = 2:1:nfirms
    if locw(i)==locw(i-1);                                    % impose symmetry  %
      x_exit_entry_w(i)        = x_exit_entry_w(i-1);          % same  x   policy %
      x_exit_entry_w(i+nfirms) = x_exit_entry_w(i-1+nfirms);   % same exit policy %
    end
  end
    
  %
%   if locw(1)==locw(2)
%     format /m1 /re 10,1;  
%     printd( 'before chkentry(): x ' ,x_exit_entry_w(1:nfirms)' , '  exit ', x_exit_entry_w(nfirms+1:2*nfirms)' , '  entry ', x_exit_entry_w(2*nfirms+1)) ;
%     format /m1 /rd 1,0 ;   
%     printd(' locw ', locw , ' doo: ' ,doo');
%   end
  %
  chkentry_w(x_exit_entry_w, w);                                             % implicitly returns v_entry(w) and isentry(w) %

  [newx(w,:), newvalue(w,:), newexit(w,:)]= optimize_w(x_exit_entry_w,w);    % THIS is the global updating of policies since the last call for each state will yield foc = 0 and hence be equilibrium for the state this iteration %
  
  % format /m1 /re 10,1; print 'x_exit_entry_w'  x_exit_entry_w' '  newx ' newx(w,:) '  exit ' newexit(w,:) '  entry ' isentry(w);  %
  foc_w = x_exit_entry_w-[newx(w,:),newexit(w,:),isentry(w)]';

  if length(doo)==0
      retp  = 1e-99 ;
  else
      retp = foc_w(doo);
  end
end
  

function contract_w(w1)
  % This function created so that one line of code needs to be called    %
  % for each state in contract(), which facilitates clean multithreading. %
  % See documentation for contract() for more details %
  % global:  x_exit_entry  foc  rflag  (among others) %
       % w2, w3, x9, y9 are only used in checking code usually turned off %
  load workspace
  locw  = dtable(:,w1);
  
  if RLG_foc<=0    % do not use Newton.  Call get_foc() just once to get best-response as in standard PM %
    
    foc(w1,:) = get_foc(x_exit_entry(w1,:)',w1,oneto2nfirms1)';                               %  get_foc() updates GLOBAL x_exit_entry and GLOBAL newx newexit isentry newval %
    
  else             % use Newton to solve foc at this state %
    temp = 1:1:nfirms;
    dosymm = temp( locw~=([-1;locw(1:nfirms-1)]) & locw>0 );            % dosymm is list of policies to change after imposing symmetry (i.e., firms with same w take same actions) %
    dosymm = [dosymm;(nfirms+dosymm);2*nfirms+1];                                               % newton search is restricted to params in dosymm to impose symmetry %
   
    do1 = dosymm( x_exit_entry(w1,dosymm)'~=0 & x_exit_entry(w1,dosymm)'~=1 );    % most 0,1 policies persist, so first do newton() holding them fixed.  Will check their optimality later %

    if ismiss(do1),  do1 = 1;    end                                                       % always search at least over the leader's x %
  
    % format /m1 /rd 2,0;  'locw ' locw'  '  dosymm ' dosymm' '  do1 ' do1'; %
    
    if length(do1)>0
      x_exit_entry(w1,do1) = newton3( x_exit_entry(w1,do1)', @get_foc, w1, do1, foc_tol)';    % updates GLOBAL x_exit_entry and GLOBAL newx newexit isentry newval %
    end
      
    % call with optimal value to ensure global newx newexit isentry newvalue properly updated with correct values, and get foc for ALL params %
    
    foc(w1,:) = get_foc(x_exit_entry(w1,:)',w1,oneto2nfirms1)';    % foc = 0 where locw = 0.  get_foc() imposes symmetry, so could use dosymm instead of oneto2nfirms1 %
    
    rflag(w1) = max(abs(foc(w1,:)'));
    
    % dpolicy = ( oldx(w1,:)~oldexit(w1,:)~oldentry(w1) )' - ( newx(w1,:)~newexit(w1,:)~isentry(w1) )'; %   % change in policies %
    
    %  % need to uncomment x9 y9 w2 w3 in local declarations %
    if 0 && iRLG==33 && w1==30         % newton failed, so get best-response functions (holding entry, exit fixed) to inspect %
      x9 = zeros(200,1);                  % 1's BR x to 2's x %
      y9 = zeros(200,1);                  % 2's BR x to 1's x %
      for w2 = 1:200
      x_exit_entry(w1,2) = w2/100;   do1 = 1;  w3 = newton3( oldexit(w1,do1), @get_foc, w1, do1, 1e-4);  x9(w2) = w3;
      x_exit_entry(w1,1) = w2/100;   do1 = 2;  w3 = newton3( oldexit(w1,do1), @get_foc, w1, do1, 1e-4);  y9(w2) = w3;
      end
%       format /m1 /rd 8,4; 
      printd [x9,y9];
      printd x9(999);  % aborts %
    end
    %
    
    if length(do1)<length(dosymm) && rflag(w1,:)>foc_tol;              % redo newton3() with do1 = dosymm if do1 did not include all of dosymm AND  % 
      do1 = dosymm;                                                 % some foc > foc_tol (either do1 did not converge OR some fixed policies need to be in do1) %
      x_exit_entry(w1,do1)  = newton3( x_exit_entry(w1,do1)', @get_foc, w1, do1, foc_tol)';
      foc(w1,:) = get_foc(x_exit_entry(w1,:)',w1,oneto2nfirms1)';
      rflag(w1) = max(abs(foc(w1,:)'));
    end
    
    if rflag(w1)>1e-2 ;
%       format /m1 /rd 4,0;  
      printd( 'foc: iter~s ', [iRLG,w1], ' w1 ', dtable(:,w1)'); 
%       format /m1 /rd 7,3;  
      printd( [newx(w1,:),-1.111,newexit(w1,:),-1.111,isentry(w1),log(max(abs(foc(w1,:)')))]);
%       format /m1 /rd 3,0;  
      printd(do1') ;
      % format /m1 /rd 4,0;   print 'x_exit_entry*  ' iRLG~w1 ' w ' dtable(:,w1)'; format /m1 /rd 7,3; print x_exit_entry_X(1:nfirms)'~-9~x_exit_entry_X(1+nfirms:2*nfirms)'~-9~x_exit_entry_X(2*nfirms+1)~log(max(abs(foc)));  format /m1 /rd 3,0; print do1'; %
    end

    if exist(chk_foc_multi_eq)         
      t0 = 0;  do1 = dosymm;  x_exit_entry_val = [x_exit_entry(w1,:),newvalue(w1,:)];
      while t0 < chk_foc_multi_eq  
	t0 = t0+1;
	x_exit_entry(w1,do1)  = newton3( LB(do1)+rand(length(do1),1).*(min([UB(do1)';10*ones(1,length(do1))])-LB(do1)), @get_foc, w1, do1, foc_tol)';    % random starting policies %
	foc_w = get_foc(x_exit_entry(w1,:)',w1,oneto2nfirms1)';
	if max(abs(foc_w))<foc_tol && max( abs( x_exit_entry(w1,do1) - x_exit_entry_val(do1) ))>.001;
	  foc_multi_eq(w1) = 1+foc_multi_eq(w1);
	  break;  % multi-eq found %
	end
      end
      x_exit_entry(w1,:) = x_exit_entry_val(1:2*nfirms+1);               % restore initial equilibrium %
      newvalue(w1,:)     = x_exit_entry_val(2*nfirms+2:3*nfirms+1);
      newx(w1,:)         = x_exit_entry_val(1:nfirms);
      newexit(w1,:)      = x_exit_entry_val(1+nfirms:2*nfirms);
      isentry(w1,:)      = x_exit_entry_val(2*nfirms+1);
    end
    
  end    % RLG_foc %
end


function contract
  % This function does one iterative step on investment and the value fn %
  % Implicit parameters are  oldx, oldexit, oldvalue (passed in)          %
  %                     and  newx, newexit, newvalue (returned)           %
 
  % First: check for which values of w_s would a firm want to enter %
  %        Entry decision is made at very beginning of each period  %
  %        BEFORE exit is implemented (done in optimize),           %
  %        BUT I think isentry should account for expected exit.    %
  %        Is this an error in this PM code (C also) ?              %
  %        chkentry() calls calcval() which integrates over invest  %
  %        outcomes, but does not implement exit, which is done in  %
  %        optimize() prior to calling calcval() when getting newx  %
  %        Of course, when I modify code to allow for random scrap  %
  %        calcval() will integrate over exit as well as invest.    %
  
  % OLD call: {newx(w,:), newvalue(w,:), newexit(w,:)}= optimize(w) %
  % Now the global new values (newx, newvalue, newexit, isentry)    %
  % are updated in  get_foc()  which is called by  optimize_w()     %
  
  % chkentry(); %     % Turn back on, and copy isentry to oldentry to restore PM updating of entry before investment/exit %

  % t0 = hsec;  'starting contract()'; %
  load workspace
  foc          = zeros(wmax,2*nfirms+1);     % these 3 are global to facilitate multithreaded code %
  rflag        = zeros(wmax,1);
  x_exit_entry = [oldx,oldexit,oldentry] ;

  Nthreads = 8;     % comment/uncomment out ThreadBegin lines below to match value of Nthreads %
  
  if nfirms<3 || Nthreads==1
      
    for w1 = 1:wmax2, contract_w( do_w(w1) );  end
  
  else      
    
    wmaxT = floor(wmax2/Nthreads);    % recall wmax2 = length(do_w) %

    parfor w =1+ (Nthreads-1)*wmaxT:wmax2,  contract_w( do_w(w) );  end     % the residual thread %

    parfor w1  =1:wmaxT,  contract_w( do_w( w1           ) );  end  
    parfor w2  =1:wmaxT,  contract_w( do_w( w2+wmaxT     ) );  end  
    parfor w3  =1:wmaxT,  contract_w( do_w( w3+wmaxT*2   ) );  end  
    parfor w4  =1:wmaxT,  contract_w( do_w( w4+wmaxT*3   ) );  end  
    parfor w5  =1:wmaxT,  contract_w( do_w( w5+wmaxT*4   ) );  end  
    parfor w6  =1:wmaxT,  contract_w( do_w( w6+wmaxT*5   ) );  end  
    parfor w7  =1:wmaxT,  contract_w( do_w( w7+wmaxT*6   ) );  end  
%
    parfor w8  =1:wmaxT,  contract_w( do_w( w8+wmaxT*7   ) );  end  
    parfor w9  =1:wmaxT,  contract_w( do_w( w9+wmaxT*8   ) );  end  
    parfor w10 =1:wmaxT,  contract_w( do_w( w10+wmaxT*9  ) );  end  
    parfor w11 =1:wmaxT,  contract_w( do_w( w11+wmaxT*10 ) );  end  
    parfor w12 =1:wmaxT,  contract_w( do_w( w12+wmaxT*11 ) );  end  
    parfor w13 =1:wmaxT,  contract_w( do_w( w13+wmaxT*12 ) );  end  
%
    ThreadJoin;
  end  % use threads %
  
  % ' seconds in contract() = ' (hsec-t0)/100; %

  if RLG_foc>0  
      w1= max(rflag);  
      if w1>1e-5 
%           format /m1 /rd 6,2;  
          printd('mean~max log10(foc) ', log([mean(rflag(do_w)),w1])); 
      end
  end

end     % Implicit returned parameters: newx, newvalue, newexit %


%
% function update
% % This function takes the solved newx, newvalue matrix for the nfirms - 1 problem %
% % and puts them into the nfirms matrices oldx, oldexit, oldvalue, for use as starting values %
% 
%   printd('This function has not been updated for use of FOC method.  Out of bounds index being used to abort...',  nfirms(2));
%   
%   if nfirms == 1
%     for i =1:wmax,
%       oldvalue(i,:) = 1 + 0.1*i;
%     end
%   else
%     for w =1:wmax
%       tuple = dtable(:,w);  % qdecode(w); %
%       nfirms = nfirms - 1;
%       n = encode2(tuple(1:nfirms));
%       oldx(w,1:nfirms)     = newx(n,1:nfirms);
%       oldexit(w,1:nfirms)  = newexit(n,1:nfirms);
%       oldvalue(w,1:nfirms) = newvalue(n,1:nfirms);
%       nfirms = nfirms + 1;
%       tuple(nfirms-1) = tuple(nfirms);
%       tuple(nfirms) = 0;
%       oldvalue(w,nfirms) = oldvalue(encode2(tuple),nfirms-1);
%       oldexit(w,nfirms)  = oldexit(encode2(tuple),nfirms-1);
%       oldx(w,nfirms)     = oldx(encode2(tuple),nfirms-1);
%       % format /m1 /rds 3,0; 'w ' dtable(:,w)';  format /m1 /rds 7,3; ' v ' oldvalue(w,:) ' x ' oldx(w,:) ' exit ' oldexit(w,:); %
%     end
%   end
%   % Implicit returned value: oldx, oldexit, oldvalue %
% end
%

function [retp1,retp2,retp3] = optimize_w(x_exit_entry,w)
% This function calculates optimal investment, and value fn., for a %
% given industry structure w. Thus, a vector nfirms long of each is returned. %
% Implicit parameters are oldx, oldvalue, isentry %
  load workspace

  locw  = dtable(:,w);  % qdecode(w); %     % profit(w,j) is used --> product market competition is BEFORE entry, exit, investment outcomes %
  locwx = locw;
  oval  = oldvalue(w,:)';
  ox    = x_exit_entry(1:nfirms);            % non-foc: oldx(w,:)' %
  oexit = x_exit_entry(nfirms+1:2*nfirms);   % non-foc: oldexit(w,:)' %
  nval  = zeros(nfirms,1);
  nx    = zeros(nfirms,1);
  nexit = zeros(nfirms,1);
  
  % dop = 0; %
  
  % if nfirms==3;  if locw(1)==7 and locw(2)==2 and locw(3)==2;  dop = 1;  end  end %
  
  if phiH==0    % FIXED scrap --> use old PM code.  Could change to use the exit mask in calcval() %
    % Find out which firms want to exit %
    disp ' HEY: probably should not run fixed scrap code anymore';
    i = (min(oval) == phi)*(minindc(oval)-1) + (min(oval) > phi)*nfirms;
    
    % Replace efficiency levels of exitors with zero     %
    if i < nfirms,    locwx(i+1:nfirms) = zeros(nfirms-i,1);   end
  end         % else integrating over random exit in calcval %

  % entry probability based on call to chkentry() just before calling optimize() %
  % if phiH==0 then isentry based on pre-exit, otherwise based on integrating over random exit in calcval %

  entered = x_exit_entry(2*nfirms+1);    % non-foc: isentry(qencode(rev(sort(locwx,1)))); %

  locwe  = locwx;               % locwe doesn't ensure room for entrant since entered = 0 if not --> locwe ignored %
  locwK = locwx;                % locw when entrant starts at kmax %
  locwe(nfirms) = entry_k;      % entrant will get moved up to its proper position by sort in calcval %
  locwk(nfirms) = kmax;
  eexit = oexit;
  eexit(nfirms) = 0;            % ensures entrant does not exit in integration in calcval() %

  % when nfirms > 2, only do leap if TWO open spots so can still check frequency of full industry for need to increase nfirms %
  doleap = (RLG_leap>0 & nfirms>1);   
  if doleap && nfirms>2,   doleap = locw(nfirms-1)==0 ;   end
  
  % Now calculate the optimal policies for this industry structure, %
  % given that entry and exit are as specified.   %
  for j =1:nfirms
    if RLG_w0_exit && locw(j) == 0              %  w=0  assumed to exit since w=0 is signal of space for entrant %
      v1 = phi;  
      if phiH ~= 0,  v1 = (phiH+phi)/2;   end     %  if locw(j) is 0 then so are all higher j %
      nval(j:nfirms) = v1*ones(nfirms-j+1,1);     %  can shut this off if entry costs, scrap chosen s.t. no entry/exit %
      nx(j:nfirms)   = zeros(nfirms-j+1,1);
      nexit(j:nfirms)= ones(nfirms-j+1,1);
      break;
    end

    v1=0; v2=0;
    if entered < 1
      % First: Calculate v, without entry %
      [v1, v2] = calcval(j,locwx,ox,oexit,locw(j));
    
      % if ddebug == 9;  'XX start';  calcvalXX(j,locwx,ox,oexit,locw(j)); ' XX done';  print /flush;  ddebug = 0;  end %
    
    end

    if entered > 0
      % A firm wants to enter with positive probability %
      % if nfirms>1;  dop=1;  end %
      [tempv1, tempv2] = calcval(j,locwe,ox,eexit,locw(j));
    
      % if ddebug == 9;  'XX start';  calcvalXX(j,locwx,ox,oexit,locw(j)); ' XX done';  print /flush;  ddebug = 0;  end   %
    
      if doleap 
	[tempk1, tempk2] = calcval(j,locwk,ox,eexit,locw(j));
	tempv1 = (1-RLG_leap)*tempv1 + RLG_leap*tempk1;
	tempv2 = (1-RLG_leap)*tempv2 + RLG_leap*tempk2;
      end
      v1 = entered*tempv1 + (1-entered)*v1;
      v2 = entered*tempv2 + (1-entered)*v2;
      % if nfirms>1; format /m1 /rds 3,0; 'j ' j ' locwe ' locwe';  format /m1 /rd 12,8; 'entered ' entered ' tempv1~v2 ' tempv1 tempv2 ' v1~v2 ' v1 v2 ;  end %
    end
    
    % Calculate values for firm, given that it is not leaving %

    if v1 <= v2
        r = 1.0;    % Avoid division by zeros %
    else
        r = 1.0/(Beta*aeff(locw(j)+1)*(v1-v2));
    end
    
    % r now contains the value r = (1 - p)^2. => p = 1 - sqrt(r)),  %
    % where p is the optimal prob. of having k rise, cond. on world %
    r = min([max([r;0.0000000000001]);1]);
    p = 1.0 - sqrt(r);
    nx(j) = p/(aeff(locw(j)+1) - aeff(locw(j)+1) * p);
  
    if BatchJob == 0 && locw(j)==kmax && locw(min(nfirms|2))>0,  nx(j) = max([RLG_minx1;nx(j)]);  end   % force min investment by frontier firms, unless monopolist %
  
    % if nx(j)< .018;  nx(j) = .018;  end %         % RLG imposed min x to avoid non convergence %

    % Now calculate the value from staying in %
    % Ask: given this optimal investment level, will there be exit? %

    % PM timing:  nval(j) = profit(w,j) - nx(j) + Beta*(v1*p + v2*(1-p));  %
    %             then check nval(j) <= phi to determine exit              %
    % RG timing:  profit(w,j) earned whether EXIT or not                   %
    %             so check  Beta*(v1*p + v2*(1-p)) - nx(j) <= phi for exit %
    nval(j) = Beta*(v1*p + v2*(1-p)) - nx(j);       % add  profit(w,j)  after exit decision %

    % 'w~j~nval(j)~profit(w,j)~nx(j)~p~v1~v2  ';     w~j~nval(j)~profit(w,j)~nx(j)~p~v1~v2; %

    if phiH==0                                     % fixed scrap %
      if nval(j) <= phi
        nexit(j:nfirms)= ones(nfirms-j+1,1);        % all weakly lower firms exit too      %
        nval(j:nfirms) = ones(nfirms-j+1,1) * phi;  % ox already 0 for lower firms         %
        nx(j) = 0  ;                                % fixed scrap --> next firms also exit %
        break;
      end
    else                                           % random scrap %
      if nval(j) < phiH                            % Pr(exit) > 0 %
        if nval(j) <= phi                          % Pr(exit) = 1 %
          nexit(j) = 1;
          nval(j)  = (phi+phiH)/2;                  % random scrap --> next j may not exit %
          nx(j)    = 0;
        else                                       % assumes scrap ~ U(phi,phiH) %
          nexit(j) = (phiH-nval(j))/(phiH-phi);     % nval(j) uses  E(scrap|exit) %
          nval(j)  = (1-nexit(j))*nval(j) + nexit(j)*(phiH+nval(j))/2;
        end 
      end
    end
    
    nval(j) = nval(j) + profit(w,j);   % RLG: profits earned even when exit at END of period %
    
    if 0          % using this leads to identical firms taking different actions: 1st of identical pair 'moves first' %
      oexit(j) = oexit(j)*RLG_damp + (1-RLG_damp)*nexit(j);    % nval(j)>= phiH --> nexit(j) remains 0, nval(j) as-is  %
      ox(j)    =    ox(j)*RLG_damp + (1-RLG_damp)*nx(j);       % ox, oexit updates optional? since oldx = newx later ? %
    end
    
    if phiH==0
      locwx(j) = (nval(j) > phi)*locw(j);      % implements new exit policy when fixed scrap %
      locwe(j) = locwx(j);
    end

  if 0
%       format /m1 /rds 3,0; 
      printd('j ', j ,' w ', locw');
%       format /m1 /rd 12,8; 
      printd(' v1~v2 ', v1, v2, ' nval ', nval', ' nx ', nx' ,' nexit ', nexit');
  end
    
    if j>1 && 0 
      if locw(j)==locw(j-1) && (abs(nval(j)-nval(j-1))>1e-8 || abs(nx(j)-nx(j-1))>1e-8 || abs(nexit(j)-nexit(j-1))>1e-8)
% 	format /m1 /rds 3,0; 
    printd('j ', j, ' w ', locw');
%     format /m1 /rds 7,3;
    printd(' v ', nval', ' x ', nx', ' newexit ', nexit');
%     format /m1 /re 10,1;
    printd(nval(j)-nval(j-1), nx(j)-nx(j-1) ,nexit(j)-nexit(j-1));
      end
    end
    
  end
  retp1 = nx';
  retp2 = nval';
  retp3 = nexit';
end


function chkentry_w(x_exit_entry,w)
% This function calculates for which value of other people's omegas, would %
% a firm want to enter, given that the market has room for another firm.    %
% Implicit parameters are oldx, oldvalue (passed in) and isentry (returned) %
   % val = Value from entering %

  % With FIXED scrap (i.e., orig PM), exit is implemented before looking up    %
  % the isentry() element, so do NOT implement exit here when getting isentry  %
  % With RANDOM scrap, exit is integrated over in calcval()                    %
  load workspace
  % for w (1,wmax,1); %     % non-foc approach called chkentry() for ALL states prior to loop over w in which call optimize(w) %
    locw = dtable(:,w);     % qdecode(w); %    % calcval integrates over exit, investment outcomes later this period %
    if locw(nfirms) == 0   % room for entrant %
      [vgarbage,v1] = calcval(nfirms,locw,x_exit_entry(1:nfirms),x_exit_entry(nfirms+1:2*nfirms),entry_k);   % non-foc:  calcval(nfirms,locw,oldx(w,:)',oldexit(w,:)',entry_k); %
      
    % if ddebug == 9;  'XX start in chkentry_w';  calcvalXX(nfirms,locw,x_exit_entry(1:nfirms),x_exit_entry(nfirms+1:2*nfirms),entry_k);  ' XX done in chkentry_w';   ddebug = 0; print /flush;  end %
    
      val = Beta * v1;

      % when nfirms > 2, only do leap if TWO open spots so can still check frequency of full industry for need to increase nfirms %
      doleap = RLG_leap>0 & nfirms>1;   
      if doleap && nfirms>2,   doleap = locw(nfirms-1)==0 ;   end
    
      if doleap        % enter at frontier with prob RLG_leap %
	[vgarbage,v1] = calcval(nfirms,locw,x_exit_entry(1:nfirms),x_exit_entry(nfirms+1:2*nfirms),kmax);
	val = (1-RLG_leap)*val + RLG_leap * Beta * v1;
      end
      % format /m1 /rd 7,3; 'w~val ' locw' val; %   % This val does not subtract entry costs -->  It's the post entry cont. value %
      v_entry(w) =  val;
      isentry(w) = max([0,min([1,((val - x_entryl) / (x_entryh - x_entryl))])]);
    end
  % end
%   isentry = minc((isentry~ones(wmax,1))');
%   isentry = maxc((isentry~zeros(wmax,1))');
  %
end


function retp = newton3(p,objfunk,w,doo,tol)
% This function performs a simple Newton-Raphson search to find the root of the function objfunk.
%   LB(doo) <= p <= UB(doo) is enforced
%   The policy parameters in p are ordered  x ~ exit ~ entry.
%   doo  has the indices from x_exit_entry being tweaked to find a root.
%   tol is param in case want to start with loose tol and tighten it as near converged value function
%

  LBw = LB;   LBw(1:nfirms) = RLG_minx1*(dtable(:,w)==kmax);     % frontier firms' min investment %
  
  np = length(p);    minstep = 1;                   % stepsize is random uniform between (minstep,1) %
  trymax= 200*np;  deriv   = zeros(np,np);
  iter = 0;        maxiter = trymax*20;
  
  while iter < maxiter
    iter=iter+1;

    x = objfunk(p,w,doo);      % Calculate function at p %

    if max(abs(x))<tol,  break;  end

    for i =1:np            % Calculate derivative matrix %
      dp = p;
      if UB(doo(i))-p(i) < 0.0001
          dp(i) = p(i)-.0001;   % decrease p to get deriv instead since too near upper bound %
      else
          dp(i) = p(i)+.0001;
      end
      deriv(:,i) = (objfunk( dp, w,doo) - x)/(dp(i)-p(i));
    end

    % if iter>2500;  format /m1 /rds 6,2;  print 'w ' w '  doo ' doo' '  iter ' iter '  p ' p' '  foc ' x';  end %

    trap 1;   % print error message?
    pnew = inv(deriv);  
    trap 0;   % print error message?

    if scalerr(pnew) || (mod(iter,trymax)==0)
%       format /m1 /rd 2,0;
      if scalerr(pnew)                         % reset iter to avoid quick restart due to iter%trymax check %
	printd('restarting newton() at rand() since singular deriv at iter ', iter, ' with w= ', dtable(:,w)', ' doo = ', doo');
%     format /m1 /re 8,1; 
    printd(' p= ', p', ' foc= ', x', ' deriv= ',deriv);   % vecr(deriv)'; %
	iter = trymax*(1+floor(iter/trymax));
      else
	printd('restarting newton() at rand() since not converging at iter ', iter, ' with w= ', dtable(:,w)', ' doo = ', doo');
%     format /m1 /re 8,1;
    printd(' p= ', p', ' foc= ', x');
      end
      pnew = LBw(doo)+rand(np,1).*(min([UB(doo)';10*ones(1,np)])-LBw(doo));
      disp ' pnew= ';disp( pnew');                          % rand between LB and UB, with cap at 10 above LB      %
      minstep = max([.3;minstep-.1]);            % singular deriv less likely if use smaller stepsize?  %
    else
      if mod(iter,trymax)<5
          pnew = p - pnew * ((.1+mod(iter,trymax)/10)*rand(1,1)) *x;     % start slow %
      else
          pnew = p - pnew *   (minstep+(1-minstep)*rand(1,1)) *x;
      end
    end
    pnew = min([UB(doo)';pnew']);
    pnew = max([LBw(doo)';pnew']);
    p = pnew;
  end
  if max(abs(x))>tol
%       format /m1 /rds 2,0;
      printd('Newton failed with w = ', dtable(:,w)', ' doo = ', doo');
%       format /m1 /re 8,1; 
      printd('  returning p= ' ,p' ,'  foc= ', x') ;
  end
  retp = p;
end


function [retp1,retp2] = calcvalXX(place,w,x,isex,k)
% This function calculates val = EEEV(:,.,.,:)p(.)p(.)p(.), where E %
% represents sums, and this is the calculation of the 4-firm problem %
% Vars: place = place of own omega, for calculating value function (v) %
%       w = the vector of omegas; already decoded %
%       x = the vector of investments (nfirms of them) %
%    isex = the vector of exit probabilities (nfirms of them) %
% Implicit parameter: oldvalue %
% For efficiency reasons, it outputs the following vector:  %
% { calcval(k_v+1,w,x), calcval(k_v,w,x) }  %

  z1 = zeros(nfirms,1);
  z2 = kmax*ones(nfirms,1);

  % ddebug = 9; %
  
  unboundedGG = 0;      % last AND checks whether 1st firm alone at frontier %
  if RLG_wstar==0 && w(place)==kmax && ( nfirms==1 )   % USED to have AND (nfirms==1 OR (w(1)-w(min(2|nfirms)))>0) %
    unboundedGG = 1;                                      % but now realize unbounded GG (RLG_wstar=0) does NOT work in oligopoly %
  end                                                  % so use bounded GG (RLG_wstar= -1) to get value and policy functions %
  % unboundedGG = 1 --> adjust valA since firm at kmax %  % but  unbounded GG in welf_ma.g to get higher dynamic CS %

  % Expand mask to allow for non-inclusion of the  place_th  firm (since not integrating over its outcomes)  %
  if nfirms > 1
    if place == 1
        locmask = [zeros(1,two_n);mask];
    elseif place == nfirms
        locmask = [mask;zeros(1,two_n)];
    else
        locmask = [mask(1:place-1,:);zeros(1,two_n);mask(place:nfirms-1,:)];
    end
  else                      locmask = zeros(1,1);
  end
  x(place) = 0;
  w(place) = k;
  isex(place) = 0;             % place firm stays put --> no investment or exit %
  justone = zeros(nfirms,1);
  justone(place) = 1;

  w1 = w;  % w1 = rev(sort(w1,1)); %     % HEY:  w1 is only used in printf statements  %

  i = aeff(w+1).*x;    p_up = i./(1+i);    % p_down = 1 - p_up; %

  valA = 0;  valB = 0;           % valA is value if w up, valB if w same %
  probmask_chksum = 0;
  
  % outer loop over    exit    outcomes %
  % inner loop over investment outcomes %
  % two_n = 2^(nfirms-1)  NOT  2^nfirms %

  if max(isex)>1
%       format /m1 /rds 3,0;
      printd(' w ', w'); 
%       format /m1 /rd 7,3; 
      printd(' isex ', isex');
  end

  % if ddebug==9;  format /m1 /rd 2,0; 'locmask (nfirms by 2^(nfirms-1))' locmask; 'G: place= ' place ' w= ' w'; format /m1 /rds 16,14;  print ' isex= ' isex' ' p_up= ' p_up';   end %

  for iEXIT = 1:two_n       % break issued at endfor if phiH==0 --> fixed scrap --> exit handled in optimize()  %
    EXprobmask = prod(2 .* locmask(:,iEXIT) .* isex + 1 - locmask(:,iEXIT) - isex);
    
  % if ddebug==9;  format /m1 /rd 2,0;  'G: iEXIT= ' iEXIT ' EXITmask= ' locmask(:,iEXIT)'  ' EXprobmask = '; format /m1 /rds 16,14; print EXprobmask '  p_up= ' p_up';   end %
  
    if 0
%       format /m1 /rd 6,0;
      disp ' ';
      disp 'two_n = ';disp( two_n);
      printd('size locmask: ', length(locmask), length(locmask'));
      printd('size    isex: ', length(isex), length(isex'));
      printd('locmask',  locmask);
    end
  
    % if nfirms>1;  if w(1)==12 AND w(2)==1 AND place==2;  if iEXIT==1; disp ' ';  end format /m1 /rds 3,0; place~9~w'~9~locmask(:,iEXIT)'; format /m1 /rds 8,4; isex'~EXprobmask;  end  end %

    % if dop;  format /m1 /rd  3,0;  'iEXIT ' iEXIT ' w ' w' ' place~k ' place k ' locmask ' locmask(:,iEXIT)' ' EXprobmask>0' (EXprobmask>0);  end  %
    
    if EXprobmask>0 || phiH==0  % else skip inner loop since zero probability %
      
      for i = 1:1:two_n
        % probmask = prod(mask(:,i) .* p_up + (1 - mask(:,i)) .* p_down); %
	probmask = prod(2 .* locmask(:,i) .* p_up + 1 - locmask(:,i) - p_up);
	
	% if ddebug==9;  format /m1 /rd 2,0;  'G: iEXIT= ' iEXIT ' i= ' i ' EXITmask= ' locmask(:,iEXIT)'  ' locmask= ' locmask(:,i)' ' EX~probmask = '; format /m1 /rds 12,10; print EXprobmask~probmask;   end %
	
	if probmask>0  % nonindented if-then %
	  
	  
	d = w+locmask(:,i);

	if phiH ~= 0                        % phiH >0 --> random scrap --> integrate over exit here in calcval %
	  d = d.*(1-locmask(:,iEXIT));    % phiH==0 --> exit handled in optimize before calling calcval      %
	  probmask = probmask*EXprobmask; % d = 0 for exiting firms %
	end
	
	probmask_chksum = probmask_chksum + probmask;    % should be 1 when done %
	
	temp = sort([d,justone],1);   % sorts via column 1 (i.e., d) so can look up in qencode %
    temp = temp(end:-1:1);
	if 0
% 	  format /m1 /rd 3,0;
	  printd('d ~ justone ~ w ~ locmask(:,i) ~ temp(2 cols) ~ isex   place= ', place, [d,justone,w,locmask(:,i),temp,isex]);
	  temp(99);
	end
	
	d = temp(:,1);                    % 2nd col of temp has 1 in row of firm j, obtained by pl1 = maxindc() %
	e = d - 1;
	% Check for evaluation of value fn. at -1 %
	e = max(([e,z1])');

	if RLG_out==0
	  if e(1)<kmax
          e = e + (e>0)*(kmax-e(1));	
      end   % all frontier firms must hav exited --> move remaining firms up s.t. highest at kmax %
	  if d(1)<kmax
          d = d + (d>0)*(kmax-d(1));	
      end   % all frontier firms must hav exited --> move remaining firms up s.t. highest at kmax %
	end
	
	if RLG_no_force_exit             % Never bump firms off lowest rung of ladder %
	  % '1 ';  'e= ' e';  'd= ' d'; 'e+d' e'+(d'.==1); %
	  e = e+(d==1);                  % This can move the focal firm into a tie with more firms, % 
	end                            % but his pl1 based on temp(:,2) still valid since tied firms have same values %

	pl1 = maxindc(temp(:,2));   % sum(d(1:place).>=k) + sum(d(place:nfirms).>k);% 
	
        % if ddebug==9;  format /m1 /rd 8,4;   'G: e = ' edisp ' '; B pl1 qencode(e) v '  pl1 qencode(e) oldvalue(qencode(e),pl1);  format /m1 /rds 16,14;  '  chksum= ' probmask_chksum;  end %
      
	if RLG_wstar <=0 && sum(d==kmax+1) > 0  % RLG: at least one firm beyond 'frontier' so use delta = 1.0 --> only use  e  %
	  valB = valB + (                                            oldvalue(qencode(e),pl1) )*probmask;
	else
	  d = min(([d,z2])');  % really only needed by PM statespace (ie, if pp_(WSTAR)>0) but harmless otherwise since previous if already checked d==kmax+1) %
	  if delta>0;
	    valB = valB + ( (1-delta)*oldvalue(qencode(d),pl1) + delta*oldvalue(qencode(e),pl1) )*probmask;
	  else
	    % if ddebug==9;  format /m1 /rd 8,4;   'G: d = ' ddisp ' '; B pl1 qencode(d) v '  pl1 qencode(d) oldvalue(qencode(d),pl1);  end %
	    valB = valB +             oldvalue(qencode(d),pl1) * probmask;
	  end
	end
	
	% format /m1 /rds 3,0;  'B: place= ' place   ' w1= ' w1'  ' d= ' d'  ' e= ' e'  ' just1= ' justone' '  pl1= ' pl1  '  k= ' k;  %
	
	
	d = w+locmask(:,i)+justone;       % +justone is the successful innovation %

	if phiH ~= 0                     % phiH >0 --> random scrap --> integrate over exit here in calcval %
	  d = d.*(1-locmask(:,iEXIT));    % phiH==0 --> exit handled in optimize before calling calcval      %
	end
	
	temp = sort([d,justone],1);
    temp = temp(end:-1:1);
	d = temp(:,1);
	e = d - 1;
	% Check for evaluation of value fn. at -1 %
	e = max(([e,z1])');
	
	if RLG_out==0
	  if e(1)<kmax
          e = e + e(e>0)*(kmax-e(1));	
      end   % all frontier firms must have exited --> move remaining firms up s.t. highest at kmax %
	  if d(1)<kmax
          d = d + d(d>0)*(kmax-d(1));
      end   % all frontier firms must have exited --> move remaining firms up s.t. highest at kmax %
	end
	
	if RLG_no_force_exit             % Never bump firms off lowest rung of ladder %
	  % '2 ';  'e= ' e';  'd= ' d'; 'e+d' e'+(d'.==1); %
	  e = e+(d==1);                  % This can move the focal firm into a tie with more firms, % 
	end                            % but his pl1 based on temp(:,2) still valid since tied firms have same values %

	pl1 = maxindc(temp(:,2));   % sum(e(1:place).>=k) + sum(e(place:nfirms).>k);%

	if RLG_wstar <= 0 && sum(d==kmax+1) > 0  % RLG: at least one firm beyond 'frontier' so use delta = 1.0 --> only use  e  %

	  % see   eql_ma.approximation_notes %

	  %valA = valA + (   % REMOVED: (1-delta)*() + delta*  %     oldvalue(qencode(e),pl1) )*probmask;
      valA = valA + (  oldvalue(qencode(e),pl1) )*probmask;

	  % if ddebug==9;  format /m1 /rd 8,4;   'G: e = ' edisp ' '; A pl1 qencode(e) v '  pl1 qencode(e) oldvalue(qencode(e),pl1)  '  here';  end  %

	
	  if unboundedGG    % only adjust valA (not valB) since approximation needed only when the 'place' firm advances beyond kmax.   nfirms==1 required for unboundedGG=1 %

	    % Unbounded GG approximation occurs here.  Bounded GG (RLG_wstar== -1) and Unbounded GG (RLG_wstar==0) with the nfirms==1 restriction %
	    
	    % valA = valA + ( ( profit(encode2(e+justone),1) - profit(qencode(e),1) ) *(1-delta) /(1-Beta) * 1  )*probmask; %    % <-- NAILS approximation when 1 firm for any delta % 
	    valA   = valA + ( ( profit(qencode(e+justone),1) - profit(qencode(e),1) ) *(1-delta) /(1-Beta) * 1  )*probmask;     % <-- NAILS approximation when 1 firm for any delta %

	    if 0 && nfirms > 1 && e(nfirms)==24
% 	      format /m1 /rds 3,0; 
          printd('place= ', place,   ' w1= ' ,w1',  ' d= ', d',  ' e= ' ,e' , ' just1= ', justone', ' e-1= ', max(([e-1,z1])')' , ' pl1= ', pl1) ;
% 	      format /m1 /rds 6,3;  
          printd(' probmask= ', probmask ,  ' v(e)= ' , oldvalue(qencode(e),pl1) ,  '  v(e-1)= ' ,oldvalue(qencode(max(([e-1,z1])')),pl1));
	      printd(' v diff=', oldvalue(qencode(e),pl1)-oldvalue(qencode(max(([e-1,z1])')),pl1));
	      printd(' profdiff~ /1-Beta= ',  ( profit(encode2(e+justone),1) - profit(qencode(e),1) ),  ( profit(encode2(e+justone),1) - profit(qencode(e),1) )/(1-Beta)*(1-delta));
	    end
	    
	    % valA = valA + ( oldvalue(qencode(e),pl1) - oldvalue(qencode(max((e-1~z1)')),pl1)  )*probmask; %          % <-- Nails approximation when 1 firm AND delta = 0 %
	    
	    % valA = valA + ( profit(encode2(e+justone),1) - profit(qencode(e),1)  + oldvalue(qencode(e),pl1)-oldvalue(qencode(max((e-1~z1)')),pl1)  )*probmask; %
	    
	    % valA = valA + ( profit(encode2(e+justone),1) - profit(qencode(e),1) ) * probmask; %
	    
	    % valA = valA + ( ( pp(MKT_SIZE)*pp(RLG_WSCALE)*(1-1/(1+exp(kmax-1))) ) *(1-delta) /(1-Beta) )*probmask; %
	    
	  end
	else
	  d = min(([d,z2])');  % really only needed by PM statespace (ie, if pp_(WSTAR)>0) but harmless otherwise since previous if already checked d==kmax+1) %
	  if delta>0
	    valA = valA + ( (1-delta)*oldvalue(qencode(d),pl1) + delta*oldvalue(qencode(e),pl1) )*probmask;
	  else
	    valA = valA +             oldvalue(qencode(d),pl1) * probmask;
	  
	  % if ddebug==9;  format /m1 /rd 8,4;  'G: d = ' d'  ' A pl1 qencode(d) v '  pl1 qencode(d) oldvalue(qencode(d),pl1);  end %

	  end
	end
	
	end  % if probmask>0 %
	%
	if dop
% 	  format /m1 /rd  3,0;  
      printd('iExit~i ', iExit, i, ' w ', w', ' place~k ', place ,k, ' locmask ', locmask(:,i)' ,' d ', d', '  e ', e', ' pl1 ' ,pl1, ' qen(d~e) ' , qencode(d), qencode(e));  
% 	  format /m1 /rd 12,8; 
      printd(' valA~B ', valA, valB,  ' probmask ' ,probmask ,' oldvalue(qencode(d~e),pl1) ' ,oldvalue(qencode(d),pl1) ,oldvalue(qencode(e),pl1)) ;
	end
	%
	
	end  % loop over investment outcomes %
	
      end     % check whether EXprobmask==0 in which case inner loop skipped %
    
      if phiH==0
          break;
      end   % exit handled by optimize() before calling calcval --> ignore outer loop %

    end    % loop over    exit    outcomes %

    if abs(1-probmask_chksum)>1e-10  
      printd(' ',  'Gauss HEY: probmask_chksum not 1: ~probmask~EXprobmask  ', probmask_chksum, probmask, EXprobmask, '  locmask next block:');
%       format /m1 /rds 3,0; 
      printd(locmask,  ' ',    'place ~ k ~ 9 ~ w ', [place,k,9,w']); 
%       format /m1 /rds 8,4;
      printd(' isex ', isex',  ' p_up ', p_up');
    end
  
  % if ddebug==9;  format /m1 /re 6,1;  'Gauss code  1-probmask_chksum= ' 1-probmask_chksum;   format /m1 /rds 8,4;  'valA~B ' valA valB;  end %
  
  retp1 = valA;
  retp2 = valB;
end


%  function (2) = calcvalXX(place,w,x,isex,k);  retp(1,1);  end %
 

function [retp1,retp2] = calcval(place,w,x,isex,k)                                              %  C version of calcval()   Hardcode the desired version by renaming the undesired version  %
% This function calculates val = EEEV(:,.,.,:)p(.)p(.)p(.), where E %
% represents sums, and this is the calculation of the 4-firm problem %
% Vars: place = place of own omega, for calculating value function (v) %
%       w = the vector of omegas; already decoded %
%       x = the vector of investments (nfirms of them) %
%    isex = the vector of exit probabilities (nfirms of them) %
% Implicit parameter: oldvalue %
% For efficiency reasons, it outputs the following vector:  %
% { calcval(k_v+1,w,x), calcval(k_v,w,x) }  %
  load workspace
  % Expand mask to allow for non-inclusion of the  place_th  firm (since not integrating over its outcomes)  %
  if nfirms > 1
    if place == 1;
        locmask = [zeros(1,two_n);mask];
    elseif place == nfirms
        locmask = [mask;zeros(1,two_n)];
    else
        locmask = [mask(1:place-1,:);zeros(1,two_n);mask(place:nfirms-1,:)];
    end
  else
      locmask = zeros(1,1);
  end
  x(place) = 0;
  w(place) = k;
  isex(place) = 0;             % place firm stays put --> no investment or exit %

  p_up = aeff(w+1).*x;    
  p_up = p_up./(1+p_up);      % p_down = 1 - p_up; %

  
  % ddebug = -1;  while ddebug ~= 0 AND ddebug ~= 9; %
  
  valA = 0;  valB = 0;         % valA is value if w up, valB if w same %
  probmask_chksum = 0;

%   dllcall call in dynamic library the function should be renamed to
%   calcval_mat
% dllcall calcval( valA, valB, probmask_chksum, place, w, x, isex, locmask, kmax, nfirms, p_up, two_n, phiH, RLG_no_force_exit, RLG_wstar, oldvalue, etable1, multfac1, binomv, delta, RLG_out);
%   calcval_mat( valA, valB, probmask_chksum, place, w, x, isex, locmask, kmax, nfirms, p_up, two_n, phiH, RLG_no_force_exit, RLG_wstar, oldvalue, etable1, multfac1, binomv, delta, RLG_out);
  
  if abs(1-probmask_chksum)>1e-10 
%     format /m1 /re 6,1; 
    printd( ' ',  'Ccode HEY: probmask_chksum not 1: 1-chksum= ', 1-probmask_chksum);
%     format /m1 /rds 3,0; 
    printd( ' ',    'place ~ k ~ 9 ~ w ' ,[place,k,9,w']);
%     format /m1 /rds 8,4; 
    printd(' isex ', isex',  ' p_up ', p_up, ' ', 'valA~B',  valA, valB);
    %
    if ddebug==1
        ddebug = 9; disp 'setting ddebug=9';     % will now call gauss version at same state for comparison %
    else
        ddebug = 1;  disp 'setting ddebug=1';     % will repeat C call with debug turned on %
    end
  else
    ddebug = 0;
   %

  end
  
  if 0
%       format /m1 /re 6,1;
      printd(' C code:  1-probmask_chksum= ', 1-probmask_chksum);
%       format /m1 /rds 12,8;  
      printd('valA~B ', valA, valB);  
  end

% end %
  
  retp1 = valA;
  retp2 = valB;
end


function retp = encode2(ntuple)    % already have a  function encode() in welf_ma.g %
% This function takes a weakly descending n-tuple (n = nfirms)       %
% with min. elt. 0, max. elt. kmax, and encodes it into an integer    %
load workspace
retp = 1+sum(binomv((ntuple+nfirms_oneton)*binomcols+ntuple+1));   % RLG: binomv is vecr(binom) and index based on gauss being row-major %
%
%   local code,digit,i;
%   code = 1;               % Coding is from 1 to wmax %
%   for i (1,nfirms,1);
%     digit = ntuple(i);
%     code = code + binom(digit+nfirms+1-i,digit+1);
%   end
%   retp(code);
%
end


function retp = qencode(ntuple)
% This function does a quick encode of any n-tuple given in weakly descending order. %
% Encoding uses a table lookup. Each column of the table consists of an n-tuple;      %
% the ith column is the ith n-tuple to be decoded.  The table is stored in etable.    %

% RLG: binomv is vecr(binom) and index based on gauss being row-major                 %
% <10% slower than old qencode, but does not need etable1 (massive for high nfirms)   %
  load workspace

  if nfirms < maxfirms_qencode 
    retp = etable1(sum(ntuple.*multfac1)+1);
  else
    retp = 1+sum(binomv((ntuple+nfirms_oneton)*binomcols+ntuple+1));    % about 10% slower than etable1 %
  end
  %
%   if nfirms <= encfirm;
%     retp = etable1(sum(ntuple.*multfac1)+1);
%   else    
%     retp = etable1(sum(ntuple.*multfac1)+1) + etable2(sum(ntuple.*multfac2)+1);
%   end
  %
end

%
% function retp = qdecode(code)
% % This function does a quick decode of a previously encoded number into   %
% % a weakly descending n-tuple. Decoding is done using a table lookup. Each %
% % column of the table consists of an n-tuple; the ith column is the ith    %
% % n-tuple to be decoded. The table is stored in the variable 'dtable'.     %
%   retp = dtable(:,code);    % Why a function call to only access a matrix?? %
% end
%

function retp = decode(code)
% This function takes a previously encoded number, and decodes it into %
% a weakly descending n-tuple (n = nfirms)                              %
% ONLY called to construct dtable (table of w's)                        %
  load workspace
  code = code-1;
  ntuple = zeros(nfirms,1);
  for i  = 1:nfirms
    digit = 0;
    while binom(digit+nfirms-i+2,digit+2) <= code
      digit=digit+1;
    end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i+1,digit+1);
  end
  retp = ntuple;
end

