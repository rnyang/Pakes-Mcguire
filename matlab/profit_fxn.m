function profit_fxn
    %% profit.g
    % LINES TO CHECK: 353
    %%

    % This is a shell for different types of profit functions %

    % new; %   % new; must be first line in program, else it terminates the program.  Hence, my idea to retain interactive PM code by running  new;  only when BatchJob==0 does not w||k %
    load('workspace');
    if BatchJob==0
        load pmg.mat;
    end
    save('workspace');
    init;
    load('workspace');

    % format /rd 12,8;

    nfmax=pp(MAX_FIRMS); % max # of active firms %
    kkmax=pp(KMAX);      % max efficiency level attainable %
    it=pp(IND_TYPE);     % investment type (quality/mc/capacity) %
    et=pp(EQL_TYPE);     % equilibrium type (Nash/monopoly/social planner) %

    dop = 0; % debugging print flag %

    % small = 1e-8;  <-- moved to init.h %

    if RLG_SH_CAP < .5
        disp '\n   SHARE CAP must be between .5 && 1 ... aborting by setting #firms = -9   \n';   
        nfmax=-9;   
    end
    RLG_fix_sh1 = 0;  fixp1 = 0;    % globals related to excluding lead firm price from newton() when imposing share cap %

    kmax=iif(et==COMPETITION, kkmax, kkmax+1);
    % ? - might change this %
    if it==CAPACITY && et==PLANNER
       kmax=kmax+1;
    end

    if RLG_wstar==0 && et==COMPETITION  
        kmax = kmax+1;  
    end  % so can access payoff when frontier improves beyond kmax %


    % clrscr();

    disp('***Computing profit function***');
    % Set up binomial coefficients for decoding/encoding of n-tuples %
    binom = eye(nfmax+kmax+2);
    binom = [zeros(nfmax+kmax+2,1),binom]; %% CHECK
    i=2;
    while i <= nfmax+kmax+2
       binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
       i=i+1;
    end

    nfirms=1;  maxp = 9999;
    while nfirms <= nfmax
        printd('Number of firms = ',compact(nfirms));
        %  Number of descending n-tuples %
        descn = binom(nfirms+kmax+1,kmax+2);
        printd('Total industry structures to compute: ', compact(descn));

        nagents = iif(et==COMPETITION, nfirms, 1);
        profit = zeros(descn,nagents);

        % The following variables are used for comparative statics %
        agprof   = zeros(descn,nfirms);  % Aggregate profits per industry structure %
        csurplus = zeros(descn,1);       % Consumer surplus %
        share    = zeros(descn,nfirms);  % Market share of each firm %
        price    = zeros(descn,nfirms);  % price of each firm %
        pmcmarg  = zeros(descn,1);       % Price/mc margin, average by sales %
        concent  = zeros(descn,1);       % One-firm concentration ratio %
        save('workspace');
        % now, call appropriate profit function %
        if it==QUALITY
           mc=pp(MC);
           M=pp(MKT_SIZE);
           wstar=pp(WSTAR);
           clear w;clear egw;clear egwp;clear p;clear profstar;
           sigma=zeros(nfirms,1);
           if 1 && nfirms==1 && RLG_y==0
               w = [kmax;1]*RLG_wscale + RLG_wshift;  % max differentiated duop %
               p = [mc+.6;mc+.5];
               p = newton(p,@cfunk); 
               maxp = pp(RLG_MAXP)*p(1);
    %            format /m1 /rds 8,4;
%                sprintf('max differentiated duopoly prices: %8.4f max monopoly price is %8.4f * p_leader = %8.4f',p',pp(RLG_MAXP),maxp);
               printd('max differentiated duopoly prices: ',p',' max monopoly price is ',pp(RLG_MAXP),' * p_leader = ',maxp)
           end
           p=(mc+.5)*ones(nfirms,1);
           if et==COMPETITION
              save('workspace');
              cqprofit(nfirms, descn);
           elseif et==MONOPOLY
              save('workspace');
              mqprofit(nfirms,descn);
           elseif et==PLANNER
              save('workspace');
              sqprofit(nfirms,descn);
           end
        elseif it==COST
           D=pp(INTERCEPT);
           f=pp(FIXED_COST);
           ggamma=pp(GAMMA);
           clear quan; clear profstar;clear w;clear theta;clear pstar;
           if et==COMPETITION
              save('workspace');
              ccprofit(nfirms, descn);
           elseif et==MONOPOLY
              save('workspace');
              mcprofit(nfirms,descn);
           elseif et==PLANNER
              save('workspace');
              scprofit(nfirms,descn);
           end
        elseif it==CAPACITY
            D=pp(INTERCEPT);
            mc=pp(MC);
            tau=pp(TAU);
            clear quan;clear profstar;clear w;clear pstar;
            if et==COMPETITION
                save('workspace');
                cpprofit(nfirms, descn);
            elseif et==MONOPOLY
                save('workspace');
                mpprofit(nfirms,descn);
            elseif et==PLANNER
                save('workspace');
                spprofit(nfirms,descn);
            end
        end

        % write output %
        if nfirms==nfmax || BatchJob==0
            [prefix,etype,itype]=acronym(et,it);
            file1 = [prefix,'pr.',compact(nfirms),'f'];
            file2 = [prefix,'cons.',compact(nfirms),'f'];
            printd('Generating output --> ',file1,file2);
    %         format /m1 /rds 16,12;
    %         screen off;
    %         output file = ^file1 reset;
            file = fopen(file1,'wt');
            fprintf(file,sprintf('%16.12f ',profit));
    %         output off;
            fclose(file);
    %         output file = ^file2 reset;
            file = fopen(file2,'wt');
            fprintf(file,'%16.12f ',agprof);
            fprintf(file,'%16.12f ',csurplus);
            fprintf(file,'%16.12f ',share);
            fprintf(file,'%16.12f ',pmcmarg);
            fprintf(file,'%16.12f ',concent);
            fprintf(file,'%16.12f ',price);
    %         screen on;
        end

        if RLG_w0==1
            Xprofit   = profit   ;  % nfirm-1 case provides values for lowest w states when RLG_w0 = 1 %
            Xagprof   = agprof   ;
            Xcsurplus = csurplus ;
            Xshare    = share    ;
            Xprice    = price    ;
            Xpmcmarg  = pmcmarg  ;
            Xconcent  = concent  ;
        end

        nfirms = nfirms + 1;
    end
    pp(PROFIT_DONE)=1;

    % save ^configfile=pp;
    save(configfile,'pp');

    if BatchJob==0
    %     run pmgshell.g;
        pmgshell
    end

% print "Done; press any key to return to main menu";  wait(); %
end % profit

function progress(i)       % report progress %
    if mod(i,1000) == 0
        sprintf('Computed: %d\n',compact(i));  % changed from \r by RLG %
    end
end



%*************************************************************************
%   Following subroutines compute static profit function for different
%   industry structure. See See section 4 of "Implementing the Pakes-McGuire
%   Algorithm for Computing Markov Perfect Equilibria in Gauss" (by Pakes,
%   Gowrisankaran, and McGuire, forthcoming), for a description of the
%   computation method of the static profit functions used.
% 
%   The following is the original code by Gautam Gowrisankaran with
%   a few modifications.
%*************************************************************************%

%*************************************************************************
%    following three functions compute static profit function for "Bertrand"
%    model (differentiated products, investment in quality) for different
%    types of industry conduct (competition, multi-plant monopolist,
%    multi-plant social planner).
%*************************************************************************%

function cqprofit(nfirms, descn)
% competition %
    % matlab check (single plant monopoly): a= 2; M=5; mc=5; y = 10;  w = ( 6 ); g= .001; p1= 8:g:10-g; n = length(p1);  u = exp(w + a*(log(y-p1) - log(y)));  s1 = u./(1+u);   prof = M*( (p1-mc).*s1 );  (blah i) = max(prof);  i=i(1); myprint(( prof(i), p1(i), s1(i)));  plot(p1,prof); %
    load('workspace');
	i = 1;
	while i <= descn
       	progress(i);
	   	w = qdecode2(i,nfirms+1);  % col vector %
		if RLG_w0==1 && length(w)== RLG_wshift
		  price(i,nfirms) = 0;                    % p was initialized as mc+.5 before calling cqprofit() %
		  if nfirms > 1
		    k = encode3(round((w(1:nfirms-1)-RLG_wshift)/RLG_wscale));      % use prices & outcomes from correponding nfirm-1 case %
		    % format /m1 /rd 3,0;  print "w " w' " k " k " w_k " qdecode2(k,nfirms); %
		    
		    profit(i,1:nfirms-1) = Xprofit(k,:);  % nfirm-1 case provides values for lowest w states when RLG_w0 = 1 %
		    agprof(i,1:nfirms-1) = Xagprof(k,:);
		    share(i,1:nfirms-1)  = Xshare(k,:);
		    price(i,1:nfirms-1)  = Xprice(k,:);
		    pmcmarg(i)           = Xpmcmarg(k);
		    csurplus(i)          = Xcsurplus(k);
		    concent(i)           = Xconcent(k);
		  end
		  
		else  % NO firms in the "out" position (i.e., w==RLG_wshift) if RLG_w0 = 1 %
		  
		  % reset starting value since w jumped back to last firm being at 0 so prev w is poor starting point %
		  if w(nfirms)==RLG_wshift    % should never be here if RLG_w0 = 1 --> above code will instead apply %
		    if RLG_y>0
                p = min( [(mc+1+.2*w') ;(RLG_y-.001*RLG_alpha)*ones(1,nfirms)] ) ;  %   "p " p;  "w " w; %
            else
                p(nfirms) = max([mc+.2;mc+1+.2*w(nfirms)/RLG_alpha]);  
                for k = nfirms-1:-1:1  
                    p(k) = p(k+1)+(w(k)-w(k+1))/RLG_alpha;  
                end
		    end
		  end

		  %
		  p(nfirms) = max([mc+.2;mc+1+.2*w(nfirms)/RLG_alpha]);  
		  for k = nfirms-1:-1:1  
              p(k) = p(k+1)+(w(k)-w(k+1))/RLG_alpha;  
          end
		  sprintf('before newton() w %d  p %d',w', p');
		  %
		  
		  % dop=0;  if nfirms==2; if w(1)==27 && w(2)==0;  dop=1;  end  end %
		  
		  if nfirms==1
		      if RLG_y>0
		          if RLG_out==0
                      p = RLG_y - small;                                                      % ignore share cap if RLG_out = 1 && nfirms = 1    %
                  else
                      p = newton( mc+.5, @cfunk);
                      if sigma(1)>RLG_sh_cap
                          p = RLG_y*(1-exp( ln( RLG_sh_cap / ((1-RLG_sh_cap)*eg(w(1))) ) /   RLG_alpha ));    % simplification of expression for fixp1 in cfunk() %
                      end
                  end
              else
                if RLG_out==0  
                    p = maxp;                                                               % ign||e share cap if RLG_out = 1 && nfirms = 1    %
                else
                    p = newton( mc+.5, @cfunk);
                    if sigma(1)>RLG_sh_cap
                        p = ln( RLG_sh_cap / ((1-RLG_sh_cap)*eg(w(1))) ) / (-RLG_alpha);                    % simplification of expression for fixp1 in cfunk() %
                    end
                end
              end
		    k= cfunk(p);
		  else
		    p = newton(p,@cfunk);

		    if sigma(1)>RLG_sh_cap                 % RLG_sh_cap >= .5 so only need to check lead firm  %
		      fixp1 = 0;                            % global, set in cfunk to yield RLG_sh_cap for lead %
		      RLG_fix_sh1 = 1;                      % flag for cfunk to know that sh_cap being imposed  %
		      p = newton( p(2:size(p,1)), @cfunk);
		      p = [fixp1;p];
		      RLG_fix_sh1 = 0;  fixp1 = 0;
		    end
		  end
		  
		  for j = 1:nfirms-1   % force firms with same w to have EXACT same price by using avg. %
		    k = (w==w(j));       % Avoids slight differences in values, policies for identical w. %
		    if sum(k)>1
              indexcat = 1:length(k);
		      k = indexcat(k==1);  % indices of firms with same w as j's %
		      p(k) = mean(p(k))*ones(size(k,1),1);
		    end
          end
		    
		  k=cfunk(p);          % call cfunk() with adjusted p to get new global sigma, used below %
		  
		  profstar = M*p.*sigma - M*mc*sigma - pp(FIXED_COST).*(w~=RLG_wshift);
		  
		  %
% 		  "nfirms     = "  nfirms;
% 		  "profit(i,.)= " profit(i,.);
% 		  "profstar'  = " profstar';
% 		  "p          = " p';
		  %
		  profit(i,:) = profstar';
		  agprof(i,:) = profstar';
		  egw = eg(w);
		  %   ###########################################   Consumer Surplus issues   #################################################     %
		  % RLG_out is just an approximation -- could equivalently have RLG_wshift really high, but might get exp() overflow issues         %
		  % TRUE absence of outside good --> utility outside good = -infinity -->  ...                                                      %
		  % --> ANY period with monopolist yields -infinity for that period && hence discounted EU is -infinity if EVER get monopolist.    %
		  % CS is infinity, or arbitrarily high, if no outside good or if approximating no outside good w/ arbitrarily high RLG_wshift.     %
		  %                                                                                                                                 %
		  % One "fix" --> redefine utility w/out the industry as 0 if RLG_y = 0 or RLG_alpha * log(RLG_y) if RLG_y > 0.  But CS<0 possible  %
		  % Another   --> redefine utility w/out the industry as utility when 1 firm at lowest w.                                           %
		  %                                                                                                                                 %
		  % Or, in welf_ma.g focus on growth in CS, perhaps as  mean CS_finalperiod / CS_initperiod.  Good for delta=0, but delta>0 ?       %
		  % This ratio works well with linear, but less so with log(y-p) since really want to measure Compensating Y needed to get disc.EU  %
		  % With log(y-p), ratio of CS_fullrun / (CS_initperiod/(1-Beta)) might be appropriate.                                             %
		  %                                                                                                                                 %
		  % With log(y-p), period CS quickly explodes as  huge y* needed for log(y*) to match  +RLG_wscale * RLG_forced (ie. grid shifts)   %
		  %                                                                                                                                 %
		  % csurplus below is  M * log(sum(exp( utility ))) /RLG_alpha  including the outside good, which therefore depends on RLG_wshift.  %
		  % If RLG_out = 1 utility for out.good, excluding log(y) if present, hardcoded as RLG_wshift-30 to avoid -infinity if monopolist.  %
		  %                                                                                                                                 %
		  % If RLG_y > 0 then in welf_ma.g will compute the period CS as  M * (exp( csurplus * RLG_alpha )/ RLG_alpha - RLG_y), since       %
		  % y* solving   alpha*log( y* ) =  log( sum( exp( utility )))   is income per person needed this period to match EU with industry. %
		  %                                                                                                                                 %
		  % If RLG_y > 0 cannot compute discounted CS by simply dividing by RLG_alpha, since marginal utility of $ not constant.  Options:  %
		  %   1. discounted sum of period CS computed as above:    M * ( exp( csurplus * RLG_alpha )/ RLG_alpha - RLG_y )                   %
		  %   2. M*(y*-y)/(1-Beta) where  alpha*log(y*)/(1-Beta) = disc.sum  Beta^t * csurplus * RLG_alpha                                  %
		  %                            --> y* = exp( (1-Beta)/alpha * sum( Beta^t * (csurplus + RLG_forced*RLG_wscale) * RLG_alpha / M )    %
		  %   1 is closer to st&&ard/linear case of discounting CS in each period, but 2 is the true compensating variation.               %
		  %                                                                                                                                 %
		  % If RLG_wstar < 0 (ie, bounded GG)  or  RLG_wstar > 0 (ie, PM)  welf_ma.g  does NOT add anything to EU = csurplus*RLG_alpha/M    %
		  % If RLG_wstar = 0 (ie, unbounded GG) then welf_ma.g adds  RLG_forced grid shifts * RLG_wscale *Beta^t to cumsum of EU            %
		  %                    where EU is csurplus*RLG_alpha.  If RLG_y==0, multiply this added amount by M/RLG_alpha to get CS addition   %
		  
		  % Arbitrary shift inside w at least 30 if not RLG_out.  Exact is +infinity --> CS = infinity %
		  if RLG_y  
              csurplus(i) = M * ln( exp(ln(RLG_y)) + sum( exp( (1-RLG_out)*max([0;10-RLG_wshift]) + ln(egw) + RLG_alpha*ln(RLG_y-p)))) /RLG_alpha;
          else
              csurplus(i) = M * ln( 1              + sum( exp( (1-RLG_out)*max([0;10-RLG_wshift]) + ln(egw) - RLG_alpha*p)) ) /RLG_alpha;
		  end 
		  share(i,:) = sigma';
		  price(i,:) = p';
		  pmcmarg(i) = sum(p.*sigma) / mc / sum(sigma);
		  concent(i) = max(sigma)/sum(sigma);
		  
		end    %  RLG_w0 == 1 %

		% print added by RLG 8/2007 to inspect more stuff %
		if nfirms<3   % || (nfirms==nfmax && w(1)==kkmax*RLG_wscale+RLG_wshift) %  ;      %  && w(1)==2.1; %
% 		  format /m1 /rd 10,3; 
%           sprintf(' w %10.3d prices %10.3f mktshares %10.3f sumshares %10.3f profit %10.3f csurplus %10.3f',w',price(i,:),share(i,:),sum(share(i,:)'),profit(i,:),csurplus(i));
          if RLG_y>0
              sprintf(' real CS %f', M*(exp(csurplus(i)*RLG_alpha/M)/RLG_alpha-RLG_y));  
          else
              sprintf(' ');
          end
		end
	      
		i = i+1;
	end

	if   0  &&  RLG_y==0 && RLG_out==0
	  % print "shifting csurplus up by min surplus with an inside good present"; %
	  csurplus = csurplus + csurplus(2);  % 2 is the 1st state with 1 firm in lowest quality position %
	end
end % cqprofit %


function mqprofit(nfirms, descn)
% monopoly %
    load('workspace');
	if RLG_alpha ~= 1  
        sprintf(' '); 
        sprintf(' HEY: need to add RLG_alpha to monopoly (cartel) case in profit.g ... aborting ');
        sprintf(' ');  
        nfirms = 0;  
    end
	  
	i = 1;
	while i <= descn
       progress(i);
	   w = qdecode2(i,nfirms+1);

 	   % As zero here represents being out, move everything down by one,  except for zero, as there is no negative efficiency level %
           w = max([(w-RLG_wscale),(RLG_wshift*ones(nfirms,1))]');

	   % matlab code to check solution:  M=5; mc=5; y = 10;  w = ( 6 4 ); g= .001; p1= 9:g:10-g; p2= p1; n = length(p1);  pp = ( repmat(p1',n,1) reshape( repmat(p2,n,1), n*n,1));  u = exp(repmat(w,n*n,1) + log(y-pp) - log(y));  s1 = u(:,1)./(1+sum(u,2));  s2 = u(:,2)./(1+sum(u,2));  prof = M*( (pp(:,1)-mc).*s1 + (pp(:,2)-mc).*s2 );  (blah i) = max(prof);  i=i(1); myprint(( prof(i), pp(i,:), s1(i), s2(i))) %
	   
	   % reset starting value since w jumped back to last firm being at 0 so prev w is po|| starting point %
	   if nfirms>1 && w(nfirms)==RLG_wshift
	     if RLG_y>0
             p = min( ([mc+.5*(w+1),(RLG_y-.001)*ones(nfirms,1)])' ) ;   % "p " p;  "w " w; %
         else
             p = (mc+1+.2*w);
	     end
	   end             %  " p before newton() " p';  %
	   
	   p = newton(p,@mfunk);
	   
           %
           if RLG_wscale == 3
               w = max([(w-3),(-7*ones(nfirms,1))]');
           else
               w = max([(w-1),(  zeros(nfirms,1))]');
           end
           %

	   profstar = m*p.*sigma - m*mc*sigma;
	   profit(i) = sum(profstar-pp(FIXED_COST.*(w~=RLG_wshift)));  % FC per active plant %
	   agprof(i,:) = profstar';
	   sprintf(' '); 
       'HEY: need to check whether fixed costs are correctly being deducted only for active plants in mqprofit() in profit.g';	   

	   egw = eg(w);   
	   if RLG_y>0
           csurplus(i) = M*ln(RLG_out + sum(exp(w+ln(RLG_y-p)-ln(RLG_y))));
       else
           csurplus(i) = M*ln(RLG_out + sum(exp(ln(egw)-p)));
	   end
	   share(i,:) = sigma';
           price(i,:) = p';
	   pmcmarg(i) = sum(p.*sigma) / mc / sum(sigma);
	   concent(i) = max(sigma)/sum(sigma);

           if nfirms<3
               pritd(' w ', w', ' prices ', p',  ' mktshares ', sigma',  ' sumshares ', sum(sigma),  ' profit ', profstar', ' CS ', csurplus(i)); 
               
           end

	   i = i+1;
	end
end % mqprofit %


function sqprofit(nfirms, descn)
% social planner %
    load('workspace');
	i = 1;
	while i <= descn
       progress(i);
 	   w = qdecode2(i,nfirms+1);

 	   % As zero here represents being out, move everything down by one,  except for zero, as there is no negative efficiency level %
           w = max(([(w-RLG_wscale),(RLG_wshift*ones(nfirms,1))])');

           %
           if RLG_wscale == 3
               w = max([(w-3),(-7*ones(nfirms,1))]');
           else
               w = max([(w-1),(  zeros(nfirms,1))]');
           end
           %

	   p = mc*ones(nfirms,1);
	   egw = eg(w);   
	   if RLG_y>0
           egwp = exp(w+ln(RLG_y-p)-ln(RLG_y));
       else
           egwp = exp(ln(egw)-p);
	   end
	   n = egwp;
	   sigma = n./(1 + sum(n));
	   profit(i) = M*ln(RLG_out + sum(egwp))-pp(FIXED_COST)*sum(w~=RLG_wshift);  % Consumer surplus - plant fixed costs %
	   sprintf(' '); 
       sprintf('HEY: need to check whether fixed costs are correctly being deducted only for active plants in sqprofit() in profit.g');
	   share(i,:) = sigma';
           price(i,:) = p';
	   pmcmarg(i) = sum(p.*sigma) / mc / sum(sigma);
	   concent(i) = max(sigma)/sum(sigma);

           if nfirms<3 
               printd(' w ',w',' prices ',p',' mktshares ',sigma','sumshares ',sum(sigma),' CS ',profit(i));  
           end

	   i = i+1;
	end
    csurplus=profit;
end % sqprofit %


function retp = cfunk(p)
% used for quality competition profit function, net of fixed costs %

	if RLG_fix_sh1  % fixing firm 1 share at RLG_sh_cap %
	  p = [1;p];       % input p excludes firm 1 price     %
	end           % 1 is filler until p(1) yielding share cap is solved for %
	
	if RLG_y>0
        n = eg(w).*exp(RLG_alpha*(ln(RLG_y-p)-ln(RLG_y)));   % eg(w) = exp(w) if RLG_wstar=0 %
    else
        n = eg(w).*exp(-RLG_alpha*p);
	end

	if RLG_fix_sh1
	  if RLG_y>0 
          fixp1 = RLG_y*(1-exp( ln( (RLG_out+sum(n(2:size(n,1))))*RLG_sh_cap / ((1-RLG_sh_cap)*eg(w(1))) ) /   RLG_alpha ));
      else
          fixp1 =               ln( (RLG_out+sum(n(2:size(n,1))))*RLG_sh_cap / ((1-RLG_sh_cap)*eg(w(1))) ) / (-RLG_alpha);      % fixp1 is global %
	  end
	  p(1) = fixp1;
	  n(1) = eg(w(1)).*exp(-RLG_alpha*p(1));
	end
	
	sigma = n./(RLG_out + sum(n));
	
	% if RLG_fix_sh1;   "shcap " RLG_sh_cap "  sigma " sigma' "  p " p' "  w " w';   end %
	
	% if nfirms==2 && w(1)==10;  print "w " w'  " prices " p'  " mktshares " sigma'  " sumshares " sum(sigma)  " exp() " n' " cfunk " (-(p-mc).*(1-sigma)./(RLG_y-p) + 1)' " profit " (M*p.*sigma-M*mc*sigma)'; 
% 	end %
	if RLG_y>0
	  foc = -RLG_alpha*(p-mc).*(1-sigma)./(RLG_y-p) + 1;    % -RLG_alpha%(RLG_y-p) is du/dp.    dsigma/dp = sigma*(1-sigma)*du/dp  where u = utility %
	else
	  % foc = -(p-mc).*sigma.*(1-sigma) + sigma;  %         % Caplin-Nalebuff FOC vect|| = 0 at solution.  divide by sigma to get rid of one of them  %
	  foc = -RLG_alpha*(p-mc).*(1-sigma) + 1;               % the -RLG_alpha* is  du/dp (ie, utility slope) %
	end

% 	%     EXPERIMENTAL CODE
% 	if RLG_sh_cap>0
% 	  i = sigma<RLG_sh_cap;
% 	  foc = foc.*i + M*RLG_sh_cap.*(sigma-RLG_sh_cap).*(1-i);                   % dprofit/dp = max quantity = M*sharecap if share exceeds sharecap.  .*(sigma-RLG_sh_cap) ensures smooth transition to foc=0 %
% 	  if max(i)==1, printd('i = ',  i', '  w=', w', '  p=', p', '  sigma= ', sigma', '  foc=', foc'); end
% 	end
	%s
	
	retp = foc(RLG_fix_sh1+1:size(foc,1));    % omit 1st value if RLG_fix_sh1 = 1 %
end

function retp = mfunk(p)
% used for quality-investment, monopoly model %
%   % Calculate the profit derivative wrt price for the multi-plant monopolist.
%   The derivative is as in (18) in PGM paper. However, we can simplify this
%   expression as follows:
%   first, divide through by sigma(n).
%   left with: FOCn = -(pn-mc)(1-sigman) + 1 + sum(k<>n)(sigmak(pk-mc))
%   second, define A = sum( sigmak(pk-mc) ).
%   left with: FOCn = -(pn-mc)(1-sigman) + 1 + A - sigman(pn-mc).
%     ==> FOCn = -(pn-mc) + 1 + A
%     ==> FOC = -(p-mc) + 1 + A.
% 
%    when RLG_y>0, the simplification is .
%    FOCn = -(pn-mc)(1-sigman)/(RLG_y-p) + 1 + sum(k<>n)(sigmak(pk-mc))/(RLG_y-p)
%    multiply by RLG_y-p && then use same A as above when RLG_y=0 
%    
    
  % Might need to fix this %

  if RLG_y>0  
    n = exp(w+ln(RLG_y-p)-ln(RLG_y));
  else
    if RLG_wstar>0
        n = eg(w).*exp(-p);
    else
        n = exp(w-p);
    end
  end
  
  sigma = n./(1.0 + sum(n));        % global sigma is used as well --> implicit returned variable %
  
  A = sum(sigma.*(p-mc));
  
  % "p= " p'; %
  
  if RLG_y>0
      retp = -(p-mc) + (RLG_y-p) + A; 
  else
      retp = -(p-mc) + 1 + A;
  end
  % retp((p-mc).*(-1 + sum(sigma)) + 1); %
end

function retp = eg(w)          
  % Calculates e^g(w)   used for quality competition profit function %
  
  % checking RLG_wstar is probably redundant now since checking RLG_wstar before calling eg() %
  % to avoid overflows with large w && large p:  exp(w).*exp(-p) less accurate than exp(w-p) %
  
  wret = exp(w);
  if RLG_wstar>0
    for i = 1:size(w,1)
      if w(i) > wstar
          wret(i) = exp(wstar)*(2.0-exp(-(w(i)-wstar)));               % w is shifted, scaled by _RLG_WSCALE, RLG_WSHIFT %
      else
          break;   % w sorted High to Low, so no need to check other i %
      end
    end
  end
  retp = wret;
end


%*************************************************************************
%    following three functions compute static profit function for "Cournot"
%    model (homogenous products, investment in marginal cost) for different
%    types of industry conduct (competition, multi-plant monopolist,
%    multi-plant social planner).
%*************************************************************************%

function ccprofit(nfirms,descn)
% competition %

	i = 1;
	while i <= descn
      progress(i);
	  w = cdecode(i,nfirms+1);
	  theta = ggamma * exp(-w);  % marginal cost %
	  % quan = solveeq(theta); % % inlined %
      % Solve for equilibrium with n firms; reduce n until all firms
      % want to produce quantity > 0 %
      n=nfirms;
      p = (D + sum(theta(1:n)))/(n+1);
      while ~(p - theta(n) >= 0) && ~(n==1) % RECHECK used to be until (p - theta(n) >= 0) || (n==1)
         n=n-1;
         p = (D + sum(theta(1:n)))/(n+1);
      end
      q = zeros(nfirms,1);
      if (p - theta(n)) > 0
         q(1:n) = p - theta(1:n);
      end
      quan=q;

	  pstar = D - sum(quan);   % Equilibrium price %
	  profstar = (pstar>theta).*(pstar-theta).*quan - f; % Equilibrium profits %
	  profit(i,:) = profstar';
	  csurplus(i) = 0.5*sum(quan)*sum(quan);
	  agprof(i,:) = profstar';
	  share(i,:) = quan';
	  if sum(quan) > 0
	    pmcmarg(i) = pstar / (sum(theta.*quan)) * sum(quan);
	    concent(i) = max(quan)/sum(quan);
	  else
        pmcmarg(i) = 1;
	  end
	  i = i+1;
	end
end % ccprofit %

function mcprofit(nfirms,descn)
% monopolist %

	i = 1;
	while i <= descn
      progress(i);
  	% The monopolist will always choose to produce everything from the lowest
% 	  priced firm, so it acts like a 1-plant firm for the static profits. %
% 	  w = cdecode(i,nfirms+1);
% 	  numin = sum(w .> -4); % No. of firms in, for fixed-fee computation %
% 	  agprof(i,:) = -f * ((w .> -4)');
% 	  w = max([(w-1),(-4*ones(nfirms,1))]');
	  % As zero here represents being out, move everything down by one, except for
	  % zero, as there is no negative efficiency level %

	  theta = ggamma * exp(-w(1));  % marginal cost %
	  pstar = 0.5*(D + theta);   % One-plant monopolist price %
	  quan = (pstar>theta)*(pstar-theta);
	  profstar = quan*(pstar-theta) - f*numin; % Monopolist profits %
	  profit(i) = profstar;
	  agprof(i,1) = agprof(i,1) + profstar;
	  csurplus(i) = 0.5*quan*quan;
	  pmcmarg(i) = pstar / theta;
	  concent(i) = 1;
	  if nfirms > 1
	    quan = [quan;zeros(nfirms-1,1)];
	  end
	  share(i,:) = quan';
      i = i+1;
	end
end % mcprofit %

function scprofit(nfirms, descn)

	i = 1;
	while i <= descn
      progress(i);
	  % The social planner will always choose to produce everything from the lowest
	  % priced firm, so it acts like a 1-plant firm for the static profits. %
	  w = cdecode(i,nfirms+1);
	  agprof(i,:) = -f*((w>-4)');
	  numin = sum(w > -4); % No. of firms in, for fixed-fee computation %
	  w = max([(w-1),(-4*ones(nfirms,1))']);
	  % As zero here represents being out, move everything down by one, except for
	  % zero, as there is no negative efficiency level %

	  theta = ggamma * exp(-w(1));  % marginal cost %
	  pstar = theta;  % Set price = mc, for social planner solution %
	  quan = (D>theta)*(D-pstar);
	  profstar = 0.5*(D-theta)*quan; % Consumer surplus %
	  profit(i) = profstar - f*numin;  % Producer surplus %
	  csurplus(i) = profstar;
	  pmcmarg(i) = 1;
	  concent(i) = 1;
	  if nfirms > 1
	    quan = [quan;zeros(nfirms-1,1)];
	  end
	  share(i,:) = quan';
	  i = i+1;
	end
end % scprofit %



%*************************************************************************
%    following three functions compute static profit function for "Capacity"
%    model (homogenous products, investment in capacity) for different
%    types of industry conduct (competition, multi-plant monopolist,
%    multi-plant social planner).
%*************************************************************************%


function cpprofit(nfirms,descn)
% competition %
	i = 1;
	while i <= descn
      progress(i);
	  w = pdecode(i,nfirms+1);
	  % quan = solveeq(w*tau); % % inlined %
      cap=w*tau;
      % Solve for equilibrium with n firms; reduce n until all firms
      % want to produce quantity > 0 %
       n=nfirms;
       sub = 0;
       q = D / (n + 1) * ones(n,1);
       while q(n) > cap(n)  % Capacity constraint violated %
          q(n) = cap(n);
          sub = sub + q(n);
          % Find reduced dem&& curve %
          n=n-1;
          if n == 0
              break; 
          end
          q(1:n) = (D - sub) / (n + 1) * ones(n,1);
      end
      quan=q;
	  pstar = D - sum(quan);   % Equilibrium price %
	  profstar = pstar.*quan; % Equilibrium profits %
	  profit(i,:) = profstar';
	  csurplus(i) = 0.5*sum(quan)*sum(quan);
	  agprof(i,:) = profstar';
	  share(i,:) = quan';
	  if sum(quan) > 0
	    pmcmarg(i) = (pstar+mc)/mc;
	    concent(i) = max(quan)/sum(quan);
	  else
        pmcmarg(i) = 1;
	  end
	  i = i+1;
	end
end % cpprofit %


function mpprofit(nfirms, descn)
% monopoly %

	% Monopolist produces quantity D/2, && charges price on dem&& curve, unless
	% total capacity is less than D/2, in which case quantity = capacity %
	i = 1;
	while i <= descn
      progress(i);
	  w = pdecode(i,nfirms+1);
	  w = max([(w-1),zeros(nfirms,1)']);
	  % As zero here represents being out, move everything down by one, except for
	  % zero, as there is no negative efficiency level %
	  % cap = w*tau;
	  % quan = solveeq(cap); % % inlined %

  	  % Solve for total production, && then allocate the production among the
      % firms putting the maximum possible production on the first firms. %
      q = zeros(nfirms,1);
      totprod = min([sum(cap);0.5*D]);
      j=1;
      while (totprod > 0) && (j <= nfirms)
         q(j) = min([totprod;cap(j)]);
         totprod = totprod - q(j);
         j=j+1;
      end
      quan=q;
      pstar = D-sum(quan);   % One-plant monopolist price %
	  profit(i) = sum(quan)*pstar;       % Monopolist profits %
	  agprof(i,:) = pstar*quan';
	  pmcmarg(i) = (pstar + mc)/mc;
	  csurplus(i) = 0.5*sum(quan)*sum(quan);
	  if sum(quan) > 0
	    concent(i) = quan(1)/sum(quan);
	    share(i,:) = quan';
	  end
	  i=i+1;
	end
end % mpprofit %


function spprofit(nfirms,descn)
% social planner%

	% Social planner produces quantity D, && charges price on dem&& curve, unless
	% total capacity is less than D, in which case quantity = capacity %
	i = 1;
	while i <= descn
      progress(i);
	  w = pdecode(i,nfirms+1);
	  w = max([(w-1),zeros(nfirms,1)']);
	  % As zero here represents being out, move everything down by one, except for
	  % zero, as there is no negative efficiency level %
	  cap = w*tau;
	  % quan = solveeq(cap); % % inlined %
      % Solve for total production, && then allocate the production among the
      % firms putting the maximum possible production on the first firms. %
      q = zeros(nfirms,1);
      totprod = min([sum(cap);D]);
      j=1;
      while (totprod > 0) && (j <= nfirms)
         q(j) = min([totprod,cap(j)]);
         totprod = totprod - q(j);
         j=j+1;
      end
      quan=q;
	  pstar = D-sum(quan);   % One-plant monopolist price %
	  agprof(i,:) = pstar*quan';
	  csurplus(i) = 0.5*sum(quan)*sum(quan);
	  profit(i) = sum(agprof(i,:)')+csurplus(i); % Actually, total surplus %
	  pmcmarg(i) = (pstar + mc)/mc;
	  if sum(quan) > 0
	    concent(i) = quan(1)/sum(quan);
	    share(i,:) = quan';
	  end
	  i=i+1;
	end
end % spprofit %



% decoding states for various models %

function retp = puttuple(ntuple,place,val)
% This function puts val in position 'place' of ntuple, && then re||ders %
% ntuple to make sure that it is in descending ||der %
  ntuple(place) = val;
  retp = sort(ntuple,1);
  retp = retp(end:-1:1);
end

function retp = qdecode2(code,nfirms)   % RLG changed to qdecode2() since welf_ma.g has a different qdecode() && batch_ma.g merges them %
% Bertr&& %
% This function takes a previously encoded number, && decodes it into %
% a weakly descending n-tuple (n = nfirms - 1)                          %
  load('workspace');
  code = code-1;
  ntuple = zeros(nfirms-1,1);
  i = 1;
  while i <= nfirms - 1
    digit = 0;
    while binom(digit+nfirms-i+1,digit+2) <= code
      digit=digit+1;
    end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i,digit+1);
    i = i+1;
  end

  % Now convert to format of starting at -7, && jumping by 3's %
  % if RLG_wscale == 3;    ntuple = (ntuple.*3  - 7);  end %
  ntuple = (ntuple.*RLG_wscale + RLG_wshift);

  retp = ntuple;
end

function retp = cdecode(code,nfirms)
% Cournot %
% This function takes a previously encoded number, && decodes it into %
% a weakly descending n-tuple (n = nfirms - 1)                          %

  code = code-1;
  ntuple = zeros(nfirms-1,1);
  i = 1;
  while i <= nfirms - 1
    digit = 0;
    while binom(digit+nfirms-i+1,digit+2) <= code
      digit=digit+1;
    end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i,digit+1);
    i = i+1;
  end

  % Now convert to format of starting at -4, && jumping by 1's %
  ntuple = ntuple-4;
  retp = ntuple;
end


function retp = pdecode(code,nfirms)
% Capacity %
% This function takes a previously encoded number, && decodes it into %
% a weakly descending n-tuple (n = nfirms - 1)                          %

  code = code-1;
  ntuple = zeros(nfirms-1,1);
  i = 1;
  while i <= nfirms - 1
    digit = 0;
    while binom(digit+nfirms-i+1,digit+2) <= code
      digit=digit+1;
    end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i,digit+1);
    i = i+1;
  end
  retp = ntuple;
end

function retp = newton( p, objfunk )
% Newton-Raphson search to find the root of the function objfunk. %
% LB <= p <= UB is enforced %

  epsilon = .0001;  tol = 1e-6;   np = size(p,1);   minstep = 1;   % stepsize is r&&om uniform between (minstep,1) %
  trymax  = 200*np;                                              % trymax = iterations before restart search at rand() %

  if pp(IND_TYPE)==QUALITY
    LB = zeros(np,1) + pp(MC);
    if pp(RLG_y) >0
        UB = zeros(np,1) + pp(RLG_y) - small;    % global small defined in init.h %
    else
        UB = zeros(np,1) + maxp;                  % global maxp set as pp(RLG_MAXP)*price of maximal leader in duopoly %
    end
  else
    disp 'HEY: in Newton() you might need to add UB && LB definitions for Cournot versions of static profit game'
    UB = zeros(np,1)+99999;
    LB = zeros(np,1);
  end
  
  maxiter = trymax*20;   iter = 0;  deriv = zeros(np,np);
  while (iter < maxiter)   % PM used:  ( (n||m > tol) || (max(abs(x))>tol/100) ) &&... but why care about p-newp when searching for foc? %
    iter = iter+1;
    x    = objfunk(p);    % Calculate function at p     %
  
    if max(abs(x))<tol
        break;  
    end     % FOC satisfied %
  
    for i = 1:np       % Calculate derivative matrix %
      dp = p;
      if UB(i)-p(i) < 0.0001
          dp(i) = p(i)-epsilon;    % decrease p to get deriv instead since too near upper bound %
      else
          dp(i) = p(i)+epsilon;
      end
      deriv(:,i) = (objfunk( dp ) - x)/(dp(i)-p(i));
    end

%     trap 1;    
    pnew = inv(deriv); 
%     trap 0;

    if (mod(iter,trymax)<4 || mod(iter,trymax)>(trymax-4)) && iter>10 
%         format /m1 /rd 7,2; 
        sprintf('Newton at iter %d with w= %d p= %d foc= ',iter,w',p');  
%         format /m1 /re 8,1; 
        sprintf('%d deriv= %d sigma= %d',x',vecr(deriv)',sigma');    
    end
    
    if scalerr(pnew) || mod(iter,trymax)==0
%       format /m1 /rd 7,2;
      if scalerr(pnew)
	iter = trymax*(1+floor(iter/trymax));   % reset iter to avoid quick restart due to iter%trymax check %
	sprintf('restarting newton() at rand() since singular deriv at iter %d with w= %d p= %d foc= ',iter,w',p');  
%     format /m1 /re 8,1; 
    sprintf('%d deriv= %d sigma= %d',x',vecr(deriv)',sigma');
      else              
	sprintf('restarting newton() at rand() since not converging at iter %d with w= %d p= %d foc= ',iter,w',p');  
%     format /m1 /re 8,1; 
    disp(x');
      end
      pnew = sort(rand(np,1),1);
      pnew = LB+pnew(end:-1:1);
      minstep = max([.1;minstep-.1]);          % non-converg less likely if use smaller stepsize? %
    else
      if mod(iter,trymax)<5
          pnew = pnew * ((.1+mod(iter,trymax)/10)*rand(1,1)) *x;     % start slow %
      else
          pnew = pnew *   (minstep+(1-minstep)*rand(1,1)) *x;     % these pnews are the change in prices %
      end
      pnew = p - pnew.*(10/max([10;abs(pnew)]));                                              % never change price by m||e than 10 in one step %

      %   THIS CODE WAS EXPERIMENTED WITH for SOLVING EQ w/ SHARE CAP
      if RLG_sh_cap
          pnew = p - pnew.*(10/max([10;abs(pnew)])).*( (RLG_sh_cap<sigma) + (RLG_sh_cap-sigma));     % take small steps if almost hitting RLG_sh_cap from below %
      else          
      end
      %
      
      if size(p,1)>1 && sigma(1)<.001 && p(1)>mc+10
          pnew(1) = max([pnew(1)-1;pnew(2)+.1]);                  
      end  % avoid min with high p, near zero share %
      if size(p,1)>1 && iter%trymax<5 && RLG_y>0
          pnew(1) = min([pnew(1);(RLG_y-.5+mod(iter,trymax)/10)]);   
      end
    end
    pnew = min([UB';pnew']);
    pnew = max([LB';pnew']);
    p = pnew;
  end
  if max(abs(x))>tol
%       format /m1 /rds 8,5; 
      sprintf('newton() failed: foc = %d  returning p= %d',x',p');  
  end
  retp = p;
end



function retp = encode3(ntuple)    % already have a  function encode() in welf_ma.g  &&  encode2(ntuple) in eql_ma.g %
% This function takes a weakly descending n-tuple (n = nfirms)       %
% with min. elt. 0, max. elt. kmax, && encodes it into an integer    %
  code = 1;                % Coding is from 1 to wmax %
  for i = 1:nfirms-1  % hardcoding nfirms-1 instead of nfirms since need index from n-1 industry.  compare to encode2() in eql_ma.g %
    digit = ntuple(i);
    code = code + binom(digit+nfirms-i,digit+1);
  end
  retp = code;
end
