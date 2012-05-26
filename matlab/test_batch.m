
%% test_batch.g
% LINES TO CHECK:
%%

%
%   created by RLG 8/2007 to run   profit.g  then  eql_ma.g  then  welf_ma.g
%   in batch mode (ie, no interaction from keyboard).
%   for loops enable processing of many different model specifications with by one gauss job.
%   
%   profit.g  eql_ma_foc.g  welf_ma.g  are merely #included into this file.
%   BatchJob==0 within these files detects whether they are being called interactively by pmgshell.g
%   or by a batch file like this one.  A few things are specific to batch/non-batch.
%

% new;
load pmg.mat;

BatchJob0 = 100;  BatchJob = BatchJob0;
% outwidth 256;

file_ma = ['batch_ma_welf_',num2str(BatchJob)];
% output file = ^file_ma reset;
file = fopen(strcat(file_ma,'.mat'),'wt');
% screen off;    output on;
fprintf(file,'  Converge Full_Indust invesT0  innov   nfirms    monop   atKMAX   run_id  EQL_TYPE  MktSize   ALPHA   WSCALE   WSHIFT    WSTAR  OUT_GOOD   RLG_y  MAXP_FACT  DELTA  MAX_FIRMS   KMAX  ENTRY w  Pr_LEAP  ENTRY:LOW  HIGH  SCRAP:Low   HIGH  INV_MULT Spillover  wstart  discEU_T log( discCS_T  CS    CS2 )    PS      innov Front_innT innov1   innov2   invest   invest1  invest0  invesT0  nfirms  #invest0  >wstar    entry   share0   share1   share2  share1T  profit1  profit2  tie_1_2    w1-w2  w1T-w2T  1st exit  markup  markup1  markup2 LoRung_stay LoRung x>0  #firms=1 #firms=2 #firms=3 #firms=4   ...   ');                           
% screen on;     output off;
fclose(file);
adj_nfirms = 1;       % if 1 then adj_nfirms as needed to avoid full industry %
MaxN = 9;             % max number of firms in any run %

% for index must be an integer   for XX (start,stop,increment)    %

for Xwstar = -1:1:-1   %  wstar = -1 for bounded GG, 0 for unbounded.  bounded vs. unbounded GG ONLY affects CS, which we currently ignore %
    for Xdelta = 0:1:0     % 0,0,1  unless comparing PM to GG.  Then 0,9,1 %
        for Xy = 15:1:15   % RLG_y  0,0,1 if linear  else  15,15,1  if log %
            for Xa = 2:-1:2    % aeff = .25 .5  1   investment efficiency      %
                for Xspill = 6:1:6    % RLG_INV =  Xspill/10    %
                    for Xentry = 16:-1:16  % 8 or 0 to 0 by -1   shift of entry costs      %
                        for Xleap = 0:1:0    % Prob leap = RLG_leap = Xleap/100              %

                            % Put the "compstat" parameter  LASt  in this list of for() loops %
                            % so can break out of it when need more than 9 firms              %

                            if BatchJob>BatchJob0
                                resultmat = [resultmat;zeros(1,rows(resultmat'))];  
                            end    % separates inner-loop results with row of zeros %

                            atMaxFirms = 0;  rlnfirms = 4;    % reset _MAX_FIRMS lower when restart inner most loop (either scale, entry, or spill for which 1st value yields fewest active firms) %

                            for Xscale = 12:-1:12   % 16,4,-1  or 8,8,-1  scaling of q, p coeff.  alpha = Xscale/10 %

                                configfile = ['_config_',num2str(BatchJob)];

                                pp = default;

                                pp(MAX_FIRMS)  = rlnfirms + ( atMaxFirms > .0001 );   %  atMaxFirms and rlnfirms from previous run, or init value at top %

                                % pp(MAX_FIRMS)  = 2;  %                       % fixed DUOP %

                                pp(RLG_LEAP)   = Xleap/100;

                                pp(RLG_y)      = Xy;                           % 0 for linear model, >0 for log(y-p) %

                                pp(INV_MULT)   =  2/(2^Xa) ;                   % investment efficiency  .25  .5  1   %

                                pp(RLG_INV)    =  Xspill/10 ;                  % proportional spillover in init.h    %

                                pp(WSTAR)      = Xwstar;                       % -1 for bounded GG, 0 for unbounded GG, >0 (say 12) for PM %
                                pp(START_FIRMS)= pp(MAX_FIRMS);
                                pp(RLG_OUTGOOD)= 0 ;                         % outside good in linear model, not in log %
                                pp(DELTA)      = Xdelta/10;

                                pp(ENTRY_LOW)   = Xentry+6;                  % +6 since +6 is default upperbound of scrap %
                                pp(SCRAP_VAL)   = pp(ENTRY_LOW)-1;          % entry is outer loop, so okay to vary it with entry costs %
                                pp(RLG_SCRAP_HI)= pp(ENTRY_LOW)  ;          % set exit just below entry %

                                pp(ENTRY_HIGH) = 2+pp(ENTRY_LOW);

                                pp(RLG_ALPHA) = Xscale/10;                   % increasing both is like decreasing var(logit error) --> higher competition %

                                pp(RLG_WSCALE)= pp(RLG_ALPHA) /2 ;

                                pp(MAX_FIRMS) = min( MaxN|pp(MAX_FIRMS));

                                % if _RLG_OUTGOOD = 1 then _RLG_WSCALE < 1 --> more competition with outside good --> need higher KMAX %
                                % One basis for KMAX is such that at p=MC (or higher, little difference) want outside share near 0                 %
                                % example:  a = .2   y = 10   p = 5    w= 4.6/a - a*(ln(y-p)-ln(y)) = 23.1  --> unfortunately too high for 5 firms %
                                % OR set pp[_RLG_OUTGOOD] = 0 %

                                pp(KMAX) = 7 + 0*(Xy==0) ;

                                pp(ENTRY_AT)  =  pp(KMAX) - 3 ;



                                pp(RLG_FOC) = 0;
                                RLG_foc = pp(RLG_FOC);  % may change in eql_ma_foc.g %

    %                             jobstart: %% DOES NOT EXIST
                                jobstart = 1;
                                while jobstart == 1
                                    jobstart = 0;

                                    pp(START_FIRMS) = pp(MAX_FIRMS);

%                                     format /m1 /rds 8,2;
                                    sprintf('  ');  
                                    disp 'BJ= ';disp(BatchJob);disp 'WSTAR= ';disp(pp(WSTAR));disp '  MAX_FIRMS= ';disp(pp(MAX_FIRMS));disp '  DELTA= ';disp(pp(DELTA));disp '  EQL_TYPE= ';disp(pp(EQL_TYPE));disp '  KMAX= ';disp(pp(KMAX));disp '  a~adj= ';disp(pp(INV_MULT));disp(pp(RLG_INV));

                                    disp '  ENTRY_AT= ';disp(pp(ENTRY_AT));disp '  OUTGOOD= ';disp(pp(RLG_OUTGOOD));disp '  ALPHA= ';disp(pp(RLG_ALPHA)); 

                                    RLG_x      = 0;       % initialize variable used to store previous model's results %
                                    RLG_exit   = 0;       % to enable tracing out equilibria which requires starting computation at previous equilibrium %
                                    RLG_value  = 0;
                                    RLG_entry  = 0;       % set these to zero to restart at zeros for each model %
                                    RLG_v_entry= 0;
                                    if exist(configfile,'file') == 2
                                        save(configfile,'pp','-append');
                                    else
                                        save(configfile,'pp');
                                    end

                                    disp 'Entering profit.g...';
                                    save('workspace');
                                    profit_fxn;
                                    load('workspace');
                                    
                                    disp 'Entering eql_ma.g...';
                                    save('workspace');
                                    eql_ma_foc;
                                    load('workspace');

                                    disp 'Entering welf_ma.g...';
                                    save('workspace');
                                    welf_ma;
                                    load('workspace');


                                    % OUTPUT FILE COMPILING RUNS %
                                    %
                                    if 0                  % old column labeled  w_noinvest %
                                      d1 = dtable(1,:)';   % find lowest frontier w for which frontier invests zero at ALL states %
                                      w1 = d1(wmax);
                                      x1 = x(:,1);         % firm 1's investment at each state.  x is the loaded policy function %
                                      y1 = sum( x1(d1==w1) );
                                      if y1>0
                                        w1 = -max(x1(d1==w1) );   % frontier invests at his highest state %
                                      else
                                        while y1==0 && w1>1
                                          w1 = w1-1;
                                          y1 = sum(x1(d1==w1) );
                                        end
                                        w1 = w1+1;
                                      end
                                    end
                                    %
                                    pp=load(configfile); 

%                                     output off;
%                                     output file=^file_ma ; 
                                    file = fopen(strcat(file_ma,'.mat'),'at');
%                                     output on;
%                                     format /m1 /rd 8,0;  
                                    fprintf(file,sprintf('%8.0f',seqa(1,1,65)'));  % column numbers %
%                                     format /m1 /rd 8,3;
                                    fprintf(file,'Converge Full_Indust invesT0  innov   nfirms    monop   atKMAX   run_id  EQL_TYPE  MktSize   ALPHA   WSCALE   WSHIFT    WSTAR  OUT_GOOD   RLG_y  MAXP_FACT  DELTA  MAX_FIRMS   KMAX  ENTRY w  Pr_LEAP  ENTRY:LOW  HIGH  SCRAP:Low   HIGH  INV_MULT Spillover  wstart  discEU_T log( discCS_T  CS    CS2 )    PS      innov Front_innT innov1   innov2   invest   invest1  invest0  invesT0  nfirms  #invest0  >wstar    entry   share0   share1   share2  share1T  profit1  profit2  tie_1_2    w1-w2  w1T-w2T  1st exit  markup  markup1  markup2 LoRung_stay LoRung x>0  #firms=1 #firms=2 #firms=3 #firms=4   ...');

                                    t = [(pp(EQL_DONE)*(1+RLG_foc+pp(RLG_MIN_INNOV1))),mean([tot_firms(:,pp(MAX_FIRMS)),totxT0,totinnov,totnfirms,tot_firms(:,1),hitkmax])',BatchJob/1000,pp(EQL_TYPE),pp(MKT_SIZE),pp(RLG_ALPHA),pp(RLG_WSCALE),pp(RLG_WSHIFT),pp(WSTAR),pp(RLG_OUTGOOD),pp(RLG_y),pp(RLG_MAXP),pp(DELTA),pp(MAX_FIRMS),pp(KMAX),pp(ENTRY_AT),pp(RLG_LEAP),pp(ENTRY_LOW),pp(ENTRY_HIGH),pp(SCRAP_VAL),pp(RLG_SCRAP_HI),pp(INV_MULT),pp(RLG_INV),sum(v(encode(wstart),:)'),mean(discEU_T),log(mean([discCS_T,consurp,consurp2])'),mean([prodsurp,totinnov,totinnovT,totinnov1,totinnov2,totinvest,totinves1,totinves0,totxT0,totnfirms,totinve0n,totwstar,totentry,totshare0,totshare1,totshare2,totshar1T,totprof1,totprof2,totwsame,totw2,totw2T,(totlspi<=numtimes),totmarkup,totprice1/pp(MC),totprice2/pp(MC),totLRstay,totLRx,tot_firms])'] ;

                                    atMaxFirms = mean(tot_firms(:,rlnfirms));

                                    if adj_nfirms
                                      if (atMaxFirms>.2 && rlnfirms>=9) || rlnfirms==MaxN 
                                        fprintf(file, 'HEY: not increasing nfirms since already at nfirms=9 and fills up more often than .2 of periods or at nfirms= 11');
                                      else
                                        if atMaxFirms>.001           % <=3 means N fixed as duop.  otherwise increase N if binding constraint %
                                          pp(MAX_FIRMS)  = rlnfirms + 1;
                                          fprintf(file,' redoing this run with higher N firms since hitting max more than .001 often');
                                          fprintf(file,num2str(t)) ;
                                          fclose(file);
%                                           goto jobstart; 
                                          jobstart = 1;
                                        end
                                      end
                                    end

                                    if jobstart == 0
                                        if BatchJob== BatchJob0
                                          resultmat = t;
                                        else
                                          i = rows(t');
                                          j = rows(resultmat');
                                          if i>j
                                              resultmat = [resultmat,(-ones(rows(resultmat),i-j))];  
                                          end
                                          if i<j
                                              t =         [t,(-ones(1,j-i))];
                                          end
                                          resultmat = [resultmat;t];
                                        end

                                        disp(resultmat);

%                                         output off;
                                        fclose(file);

                                        BatchJob = BatchJob+1;  % becomes part of prefix of all filenames so output retained across runs %
                                        if   0 &&   rlnfirms<9 && mean(totxT0)>.9 && pp(RLG_MIN_INNOV1)==0       % only do if rlnfirms <9 since takes too long for potentially low payoff %
                                          sprintf('  ');  disp '##########   Re-do with  _RLG_MIN_INNOV1 = .01  #############';
                                          pp(RLG_MIN_INNOV1) = .01;
                                          pp(EQL_DONE)= 0;
    %                                       goto jobstart; %% DOES NOT 
                                          jobstart = 1;
                                        end
                                        if jobstart == 0
                                            if atMaxFirms>.001 && rlnfirms>2   && 0       % inner-most loop is the comp.stat for which #firms is increasing %
                                              sprintf('  ');  disp 'Breaking inner-loop since too many active firms already --> later runs in this loop will be worthless';  sprintf('  ');
                                              break;
                                            end 
                                        end
                                    end
                                end % while
                            end
                        end
                    end
                end
            end
        end
    end
end




% function retp = default
%     p=zeros(NPARMS,1);
%     p(MAX_FIRMS)= 3;
%     p(START_FIRMS)= 1;
%     p(EQL_TYPE)= COMPETITION;
%     p(IND_TYPE)= QUALITY;
%     p(ENTRY_TYPE)= RAN_ENTRY;
%     p(ENTRY_LOW)= 0.3;
%     p(ENTRY_HIGH)= 0.5;
%     p(ENTRY_SUNK)= 0.2;
%     p(ENTRY_AT)= 5;
%     p(Beta)= 0.925;
%     p(DELTA)= 0.7;
% 
%     p(INV_MULT)= 3;
%     p(INV_COST)= 1;
%     p(MC)= 5;
%     p(MKT_SIZE)= 5;
%     p(KMAX)= 15;
%     p(WSTAR)= 12;           % 12 for Rand94 ? %
%     p(INTERCEPT)= 3;
%     p(FIXED_COST)= 0.0;
%     p(GAMMA)= 1;
%     p(TAU)= 0.1;
% 
%     p(RLG_OUTGOOD)=  1;
%     p(RLG_ALPHA)=    1;
%     p(RLG_WSCALE)=   1;    %  3 for Rand94 ? %
%     p(RLG_WSHIFT)=   0;    % -7 for Rand94 ? %
%     p(RLG_SH_CAP)=   1;    % .55, .65 for Rand94 %
%     p(RLG_INV)=   0;       % rate at which investment efficiency increases in w %
%     p(RLG_y)=  10;         % -alpha*log(RLG_y - price) %
%     p(RLG_LEAP)=  0;       % prob entrant starts at KMAX %
%     p(RLG_MAXP)=  1.5;     % maximimum monopolist price =  RLG_MAXP*max_duop_price in competitive equilibrium when RLG_y =  0 (ie, linear).  default =  9999 %
%     p(RLG_FOC) = 0;        % if 0, use FOC method in  eql_ma_foc.g (where this value is hardcoded)  ONLY if non-FOC fails to converge.   If >0, do not try non-FOC.  If <0, never use FOC %
%     p(RLG_MIN_INNOV1)=0;   % minimum innovation rate by frontier firms.  Only used if RLG_WSTAR <= 0 (i.e., not using wstar of original PM) %
% 
%     p(SCRAP_VAL)   =  5;   % lower bound of scrap %
%     p(RLG_SCRAP_HI)=  6;   % upper bound of scrap value, 0 for fixed scrap %
% 
%     p(PROFIT_DONE)= 0;
%     p(EQL_DONE)= 0;
%     p(ACTIVE_CFG)= 'default';
%     retp = p;
% end % default %
