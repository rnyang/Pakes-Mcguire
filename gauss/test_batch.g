/*
  created by RLG 8/2007 to run   profit.g  then  eql_ma.g  then  welf_ma.g
  in batch mode (ie, no interaction from keyboard).
  for loops enable processing of many different model specifications with by one gauss job.
  
  profit.g  eql_ma_foc.g  welf_ma.g  are merely #included into this file.
  BatchJob==0 within these files detects whether they are being called interactively by pmgshell.g
  or by a batch file like this one.  A few things are specific to batch/non-batch.
*/

new;
#include pmg.h;

BatchJob0 = 100;  BatchJob = BatchJob0;
outwidth 256;

file_ma = "batch_ma_welf_" $+ ftos(BatchJob,"%*.*lf",1,0);
output file = ^file_ma reset;   screen off;    output on;
" "; 
"Converge Full_Indust invesT0  innov   nfirms    monop   atKMAX   run_id  EQL_TYPE  MktSize   ALPHA   WSCALE   WSHIFT    WSTAR  OUT_GOOD   RLG_y  MAXP_FACT  DELTA  MAX_FIRMS   KMAX  ENTRY w  Pr_LEAP  ENTRY:LOW  HIGH  SCRAP:Low   HIGH  INV_MULT Spillover  wstart  discEU_T log( discCS_T  CS    CS2 )    PS      innov Front_innT innov1   innov2   invest   invest1  invest0  invesT0  nfirms  #invest0  >wstar    entry   share0   share1   share2  share1T  profit1  profit2  tie_1_2    w1-w2  w1T-w2T  1st exit  markup  markup1  markup2 LoRung_stay LoRung x>0  #firms=1 #firms=2 #firms=3 #firms=4   ...";
" ";                            screen on;     output off;

adj_nfirms = 1;       @ if 1 then adj_nfirms as needed to avoid full industry @
MaxN = 9;             @ max number of firms in any run @

@ for index must be an integer   for XX (start,stop,increment)    @

for Xwstar (-1,-1,1);   @  wstar = -1 for bounded GG, 0 for unbounded.  bounded vs. unbounded GG ONLY affects CS, which we currently ignore @
for Xdelta (0,0,1);     @ 0,0,1  unless comparing PM to GG.  Then 0,9,1 @
for Xy     (15,15,1);   @ RLG_y  0,0,1 if linear  else  15,15,1  if log @
for Xa     (2,2,-1);    @ aeff = .25 .5  1   investment efficiency      @
for Xspill (6,6, 1);    @ RLG_INV =  Xspill/10    @
for Xentry (16,16,-1);  @ 8 or 0 to 0 by -1   shift of entry costs      @
for Xleap  (0,0, 1);    @ Prob leap = RLG_leap = Xleap/100              @

@ Put the "compstat" parameter  LASt  in this list of for() loops @
@ so can break out of it when need more than 9 firms              @

if BatchJob>BatchJob0;  resultmat = resultmat|zeros(1,rows(resultmat'));  endif;    @ separates inner-loop results with row of zeros @

atMaxFirms = 0;  rlnfirms = 4;    @ reset _MAX_FIRMS lower when restart inner most loop (either scale, entry, or spill for which 1st value yields fewest active firms) @
  
for Xscale (12,12,-1);   @ 16,4,-1  or 8,8,-1  scaling of q, p coeff.  alpha = Xscale/10 @

configfile = "_config_" $+ ftos(BatchJob,"%*.*lf",1,0);

pp = default();

pp[_MAX_FIRMS]  = rlnfirms + ( atMaxFirms > .0001 );   @  atMaxFirms and rlnfirms from previous run, or init value at top @

@ pp[_MAX_FIRMS]  = 2;  @                       @ fixed DUOP @

pp[_RLG_LEAP]   = Xleap/100;

pp[_RLG_y]      = Xy;                           @ 0 for linear model, >0 for log(y-p) @

pp[_INV_MULT]   =  2/(2^Xa) ;                   @ investment efficiency  .25  .5  1   @

pp[_RLG_INV]    =  Xspill/10 ;                  @ proportional spillover in init.h    @

pp[_WSTAR]      = Xwstar;                       @ -1 for bounded GG, 0 for unbounded GG, >0 (say 12) for PM @
pp[_START_FIRMS]= pp[_MAX_FIRMS];
pp[_RLG_OUTGOOD]= 0 ;                         @ outside good in linear model, not in log @
pp[_DELTA]      = Xdelta/10;

pp[_ENTRY_LOW]   = Xentry+6;                  @ +6 since +6 is default upperbound of scrap @
pp[_SCRAP_VAL]   = pp[_ENTRY_LOW]-1;          @ entry is outer loop, so okay to vary it with entry costs @
pp[_RLG_SCRAP_HI]= pp[_ENTRY_LOW]  ;          @ set exit just below entry @

pp[_ENTRY_HIGH] = 2+pp[_ENTRY_LOW];

pp[_RLG_ALPHA] = Xscale/10;                   @ increasing both is like decreasing var(logit error) --> higher competition @

pp[_RLG_WSCALE]= pp[_RLG_ALPHA] /2 ;

pp[_MAX_FIRMS] = minc( MaxN|pp[_MAX_FIRMS]);

@ if _RLG_OUTGOOD = 1 then _RLG_WSCALE < 1 --> more competition with outside good --> need higher KMAX @
@ One basis for KMAX is such that at p=MC (or higher, little difference) want outside share near 0                 @
@ example:  a = .2   y = 10   p = 5    w= 4.6/a - a*(ln(y-p)-ln(y)) = 23.1  --> unfortunately too high for 5 firms @
@ OR set pp[_RLG_OUTGOOD] = 0 @

pp[_KMAX] = 7 + 0*(Xy==0) ;

pp[_ENTRY_AT]  =  pp[_KMAX] - 3 ;



pp[_RLG_FOC] = 0;
RLG_foc = pp[_RLG_FOC];  @ may change in eql_ma_foc.g @

jobstart:


pp[_START_FIRMS] = pp[_MAX_FIRMS];

format /m1 /rds 8,2;
" ";  "BJ= " Batchjob " WSTAR= " pp[_WSTAR] "  MAX_FIRMS= " pp[_MAX_FIRMS] "  DELTA= " pp[_DELTA] "  EQL_TYPE= " pp[_EQL_TYPE]  "  KMAX= " pp[_KMAX]  "  a~adj= " pp[_INV_MULT]  pp[_RLG_INV];;
"  ENTRY_AT= " pp[_ENTRY_AT]  "  OUTGOOD= " pp[_RLG_OUTGOOD]  "  ALPHA= " pp[_RLG_ALPHA]; 

RLG_x      = 0;       @ initialize variable used to store previous model's results @
RLG_exit   = 0;       @ to enable tracing out equilibria which requires starting computation at previous equilibrium @
RLG_value  = 0;
RLG_entry  = 0;       @ set these to zero to restart at zeros for each model @
RLG_v_entry= 0;

save ^configfile=pp;
  
"Entering profit.g...";
#include profit.g ;

"Entering eql_ma.g...";
#include eql_ma_foc.g ;


"Entering welf_ma.g...";
#include welf_ma.g ;


/* OUTPUT FILE COMPILING RUNS */
/*
if 0;                  @ old column labeled  w_noinvest @
  d1 = dtable[1,.]';   @ find lowest frontier w for which frontier invests zero at ALL states @
  w1 = d1[wmax];
  x1 = x[.,1];         @ firm 1's investment at each state.  x is the loaded policy function @
  y1 = sumc( selif(x1,d1.==w1) );
  if y1>0;
    w1 = -maxc(selif(x1,d1.==w1) );   @ frontier invests at his highest state @
  else;
    do while y1==0 AND w1>1;
      w1 = w1-1;
      y1 = sumc( selif(x1,d1.==w1) );
    endo;
    w1 = w1+1;
  endif;
endif;
*/
load pp=^configfile;

output off;
output file=^file_ma ;
output on;
format /m1 /rd 8,0;  print seqa(1,1,65)';  @ column numbers @
format /m1 /rd 8,3;
"Converge Full_Indust invesT0  innov   nfirms    monop   atKMAX   run_id  EQL_TYPE  MktSize   ALPHA   WSCALE   WSHIFT    WSTAR  OUT_GOOD   RLG_y  MAXP_FACT  DELTA  MAX_FIRMS   KMAX  ENTRY w  Pr_LEAP  ENTRY:LOW  HIGH  SCRAP:Low   HIGH  INV_MULT Spillover  wstart  discEU_T log( discCS_T  CS    CS2 )    PS      innov Front_innT innov1   innov2   invest   invest1  invest0  invesT0  nfirms  #invest0  >wstar    entry   share0   share1   share2  share1T  profit1  profit2  tie_1_2    w1-w2  w1T-w2T  1st exit  markup  markup1  markup2 LoRung_stay LoRung x>0  #firms=1 #firms=2 #firms=3 #firms=4   ...";

t = (pp[_EQL_DONE]*(1+RLG_foc+pp[_RLG_MIN_INNOV1]))~meanc(tot_firms[.,pp[_MAX_FIRMS]]~totxT0~totinnov~totnfirms~tot_firms[.,1]~hitkmax)'~BatchJob/1000~pp[_EQL_TYPE]~pp[_MKT_SIZE]~pp[_RLG_ALPHA]~pp[_RLG_WSCALE]~pp[_RLG_WSHIFT]~pp[_WSTAR]~pp[_RLG_OUTGOOD]~pp[_RLG_y]~pp[_RLG_MAXP]~pp[_DELTA]~pp[_MAX_FIRMS]~pp[_KMAX]~pp[_ENTRY_AT]~pp[_RLG_LEAP]~pp[_ENTRY_LOW]~pp[_ENTRY_HIGH]~pp[_SCRAP_VAL]~pp[_RLG_SCRAP_HI]~pp[_INV_MULT]~pp[_RLG_INV]~sumc(v[encode(wstart),.]')~meanc(discEU_T)~log(meanc(discCS_T~consurp~consurp2)')~meanc(prodsurp~totinnov~totinnovT~totinnov1~totinnov2~totinvest~totinves1~totinves0~totxT0~totnfirms~totinve0n~totwstar~totentry~totshare0~totshare1~totshare2~totshar1T~totprof1~totprof2~totwsame~totw2~totw2T~(totlspi.<=numtimes)~totmarkup~totprice1/pp[_MC]~totprice2/pp[_MC]~totLRstay~totLRx~tot_firms)' ;

atMaxFirms = meanc(tot_firms[.,rlnfirms]);

if adj_nfirms;
  if (atMaxFirms>.2 AND rlnfirms>=9) OR rlnfirms==MaxN ;
    "HEY: not increasing nfirms since already at nfirms=9 and fills up more often than .2 of periods or at nfirms= 11";
  else;
    if atMaxFirms>.001 ;          @ <=3 means N fixed as duop.  otherwise increase N if binding constraint @
      pp[_MAX_FIRMS]  = rlnfirms + 1;
      " redoing this run with higher N firms since hitting max more than .001 often";
      print t ;
      output off;
      goto jobstart;
    endif;
  endif;
endif;


if BatchJob== BatchJob0;
  resultmat = t;
else;
  i = rows(t');
  j = rows(resultmat');
  if i>j;  resultmat = resultmat~(-ones(rows(resultmat),i-j));  endif;
  if i<j;          t =         t~(-ones(1,j-i));                endif;
  resultmat = resultmat|t;
endif;
  
print resultmat;

output off;

BatchJob = BatchJob+1;  @ becomes part of prefix of all filenames so output retained across runs @

if   0 AND   rlnfirms<9 AND meanc(totxT0)>.9 AND pp[_RLG_MIN_INNOV1]==0;       @ only do if rlnfirms <9 since takes too long for potentially low payoff @
  " ";  "##########   Re-do with  _RLG_MIN_INNOV1 = .01  #############";
  pp[_RLG_MIN_INNOV1] = .01;
  pp[_EQL_DONE]= 0;
  goto jobstart;
endif;

if atMaxFirms>.001 AND rlnfirms>2   AND 0 ;      @ inner-most loop is the comp.stat for which #firms is increasing @
  " ";  "Breaking inner-loop since too many active firms already --> later runs in this loop will be worthless";  " ";
  break;
endif;

endfor;
endfor;
endfor;
endfor;
endfor;
endfor;
endfor;
endfor;




proc (1) = default();
   local p;

   p=zeros(NPARMS,1);
   p[_MAX_FIRMS]= 3;
   p[_START_FIRMS]= 1;
   p[_EQL_TYPE]= _COMPETITION;
   p[_IND_TYPE]= _QUALITY;
   p[_ENTRY_TYPE]= _RAN_ENTRY;
   p[_ENTRY_LOW]= 0.3;
   p[_ENTRY_HIGH]= 0.5;
   p[_ENTRY_SUNK]= 0.2;
   p[_ENTRY_AT]= 5;
   p[_BETA]= 0.925;
   p[_DELTA]= 0.7;

   p[_INV_MULT]= 3;
   p[_INV_COST]= 1;
   p[_MC]= 5;
   p[_MKT_SIZE]= 5;
   p[_KMAX]= 15;
   p[_WSTAR]= 12;           /* 12 for Rand94 ? */
   p[_INTERCEPT]= 3;
   p[_FIXED_COST]= 0.0;
   p[_GAMMA]= 1;
   p[_TAU]= 0.1;

   p[_RLG_OUTGOOD]=  1;
   p[_RLG_ALPHA]=    1;
   p[_RLG_WSCALE]=   1;    /*  3 for Rand94 ? */
   p[_RLG_WSHIFT]=   0;    /* -7 for Rand94 ? */
   p[_RLG_SH_CAP]=   1;    /* .55, .65 for Rand94 */
   p[_RLG_INV]=   0;       /* rate at which investment efficiency increases in w */
   p[_RLG_y]=  10;         /* -alpha*log(RLG_y - price) */
   p[_RLG_LEAP]=  0;       /* prob entrant starts at KMAX */
   p[_RLG_MAXP]=  1.5;     /* maximimum monopolist price =  RLG_MAXP*max_duop_price in competitive equilibrium when RLG_y =  0 (ie, linear).  default =  9999 */
   p[_RLG_FOC] = 0;        /* if 0, use FOC method in  eql_ma_foc.g (where this value is hardcoded)  ONLY if non-FOC fails to converge.   If >0, do not try non-FOC.  If <0, never use FOC */
   p[_RLG_MIN_INNOV1]=0;   /* minimum innovation rate by frontier firms.  Only used if RLG_WSTAR <= 0 (i.e., not using wstar of original PM) */

   p[_SCRAP_VAL]   =  5;   /* lower bound of scrap */
   p[_RLG_SCRAP_HI]=  6;   /* upper bound of scrap value, 0 for fixed scrap */

   p[_PROFIT_DONE]= 0;
   p[_EQL_DONE]= 0;
   p[_ACTIVE_CFG]= "default";
   retp(p);
endp; /* default */
