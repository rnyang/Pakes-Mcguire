@ This is a slightly modified version of markov.g program that allows
for reading parameters from the configuration file.
Version 18 uses the Allw x nfirms coding method, and allows for starting
at any number of firms.
Written by: Gautam Gowrisankaran
May 16, 1993

  Modified by Ron Goettler April 2011 to solve FOC within each state each iteration.
  Modification consists of adding a function to return foc for  1) entry prob  2) each firm's x  3) each firm's exit prob.
  The foc are obtained simply by subtracting the candidate values from their optimal values as determined within optimize().
  
  To run in iterative BR mode (i.e., standard PM), the code calls  get_foc()  once, but doesn't call the Newton solver.
  One minor difference between this iterative BR and PM is that  isentry  is updated along with other policies, not beforehand as in PM.
  Hence, RLG added  oldentry.  The variable  isentry  is retained, instead of being renamed to  newentry.
  If needed, the old PM timing of updating  isentry  prior to updating investments/exit could be reinstated without too much coding.
  But the foc approach requires solving for entry simultaneously with investment/exit, so cleaner to update entry similarly when not using foc.
  
  When both FOC and non-FOC converge, they appear to converge to the same point.  But sometimes FOC converged when non-FOC doesn't, and vice-versa.
  if RLG_foc = 0, use FOC method in  eql_ma_foc.g (where this value is hardcoded)  ONLY if non-FOC fails to converge.   If >0, do not try non-FOC.  If <0, never use FOC.
@

@ new; @   @ new; must be first line in program, else it terminates the program.  Hence, my idea to retain interactive PM code by running  new;  only when BatchJob==0 does not work @

if BatchJob==0;  
#include pmg.h;
endif;

x9 = 0;  y9=0;  @ globals for storing best-response functions when invesigating non-convergence of nfirm = 2 case @

clrscr();   "**** Computing Dynamic Equilibrium using  eql_ma_foc.g  ****";

#include init.h;

dlibrary ./calcval ;

ddebug = 0;

@ constants not modifiable by user @
tol = 0.001;        @ Tolerance for convergence of value function @

foc_tol = 1e-6;     @ Tolerance for convergence of FOC solution within each state @

rndseed 58921;      @ random restarts in newton3() @

chk_foc_multi_eq = 0;  @ # of random restarts of within state foc search to attempt to find a different solution @

RLG_wstar  = pp[_WSTAR];

RLG_foc = pp[_RLG_FOC];       @ see documentation of RLG_foc at end of preamble @

RLG_damp   = 0;      "weight on old values is RLG_damp  = " RLG_damp;

dop = 0;     @ do print flag for debugging @

if 1;   RLG_w0_exit = 1;      "enforcing exit at w = 0 (as done in original PM code to signal room for entrants) ";
else;   RLG_w0_exit = 0;  "NOT enforcing exit at w = 0 (as done in original PM code to signal room for entrants.  Okay if no entry/exit) ";
endif;

if phiH==0; " ";  " HEY: must revise fixed scrap to address profits earned in exit period --> nval[] for exitors depends on j in optimize !! ... deliberately forcing abort with index out of range";   pp[9999999];  endif;

if pp[_RLG_SCRAP_HI]-pp[_SCRAP_VAL] < .2;  " HEY: solving FOC in each state may be difficult since range of scrap values is less than .2"; " ";  endif;


clear newvalue,newx,newexit,oldvalue,oldx,oldexit,isentry,v_entry;


@ Set up binomial coefficients for decoding/encoding of n-tuples @

if RLG_wstar == 0;  kmax = kmax+1;  endif;      @ need tables 1 beyond kmax @

binom = eye(rlnfirms+kmax+1);
binom = zeros(rlnfirms+kmax+1,1)~binom;
for i (2,rlnfirms+kmax+1,1);
  binom[i,2:i] = binom[i-1,2:i] + binom[i-1,1:i-1];
endfor;
binomv    = vecr(binom);      @ for RLG version of encode() that avoids for loop & @
binomcols = rows(binom');     @ avoids etable1 (huge if nfirms big).  ~10% slower  @
maxfirms_qencode = 9;         @ maximum nfirms for which etable1 constructed for use in qencode() @
etable1 = 0;  multfac1 = 0;   @ initialize here in case nfirms > maxfirms_qencode  @

if RLG_wstar == 0;  kmax = kmax-1;  endif;      @ RESTORE kmax @


@ RLG asks what is value of encfirm < rlnfirms ??  The execution time is same.  And eql_sa.g crashes if encfirm < rlnfirms @
encfirm = rlnfirms;      @ was hardcoded to something like 4,5, or 6.   Max. number of firms to encode in table @

filenam1 = prefix $+ "markov.dp";
print "Execution protocol --> " filenam1;
output file = ^filenam1 reset;

nfirms = stfirm;        @ max. # of firms at initial computation stage @
nfirms = rlnfirms;      @ hardcoded ignoring of option to use equilibrium from industry with fewer firms as starting point @
oneton = seqa(1,1,nfirms);
oneto2nfirms1 = seqa(1,1,2*nfirms+1);

LB = zeros(1+2*nfirms,1);  @ bounds on policies:  x in 1:nfirms, exit in nfirms+1:2*nfirms, entry in 1+2*nfirms @
UB =  ones(1+2*nfirms,1);  UB[1:nfirms] = UB[1:nfirms] + 9999999999;  @ investment has no upper-bound @


@ Next blocks are somewhat messy since PM allowed for possibility of using equilibrium with fewer firms as starting values  @
@ RLG ignores this option and also sets encfirm = nfirms.  Hence, multfac2 is not used and the opening block is skipped.    @
@ Much of code in next few blocks could therefore be removed, but is retained in case desire to use it again someday.       @
@ RLG also added saving and loading of  dtable  etable1  mask, since building them is SLOW (nested loops) when many firms.  @

if   0   AND nfirms > 1;   @ Skip this block when don't want to use model with fewer firms as starting values @
  nfirms = nfirms - 1;
  wmax = binom[nfirms+1+kmax,kmax+2];
  if 0;
    @ Read data: v (value), x (investment), p (probability of state rising), isentry, newexit @
    filename = prefix $+ "markov." $+ padr(nfirms,1,0) $+ "ot";
    load bigread[] = ^filename;
    newvalue = (reshape(bigread[1:wmax*nfirms],wmax,nfirms));
    bigread = bigread[wmax*nfirms+1:rows(bigread)];
    newx = (reshape(bigread[1:wmax*nfirms],wmax,nfirms));
    if phiH /= 0;
      bigread = bigread[wmax*nfirms+1:rows(bigread)];
      newexit = (reshape(bigread[1:wmax*nfirms],wmax,nfirms));
    endif;
    clear bigread;               @ isentry not read in since it's computed first @
  else;
    newx     = zeros(wmax,nfirms);
    newexit  = zeros(wmax,nfirms);
    newvalue = zeros(wmax,nfirms);
  endif;
  oneton = seqa(  1,1,nfirms);

  if nfirms >= encfirm;
    multfac1 = (kmax+1)^(oneton[1:encfirm]-1);
    nfirms = encfirm;
    @ Encode all numbers from 1 to kmax^nfirms @
    etable1 = zeros((kmax+1)^nfirms,1);
    for i (1, rows(etable1), 1);     @ was:  i=0;  do while i < rows(etable1); ... i= i+1; endo; @
      clear msk; @ Build one mask @
      j=kmax+1; k=i-1;
      do while j <= rows(etable1);
        msk = ((k % j)*(kmax+1)/j)|msk;
        k = k - (k % j);
        j = j*(kmax+1);
      endo;
      etable1[i] = encode2(rev(sortc(msk[1:nfirms],1)));
    endfor;
    nfirms = stfirm-1;
    if nfirms > encfirm;
      multfac1 = zeros(nfirms-encfirm,1)|multfac1;
    endif;
  endif;
  nfirms = stfirm;
endif;

do while nfirms <= rlnfirms;
  nfirms_oneton = nfirms-seqa(1,1,nfirms);     @ used in encode3() and qencode() @

  filenam2 = prefix $+ "markov." $+ ftos(nfirms,"%*.*lf",1,0) $+ "ot";
  d1 = date;
  format /rd 9,6;
  @ Number of different combinations of competitors positions faced by a firm @
  wmax = binom[nfirms+1+kmax,kmax+2];
  print "Number of firms: " compact(nfirms);;
  print "; Total number of states: " compact(wmax);
  print "Initializing ...";;

  filename = prefix $+ "pr." $+ ftos(nfirms,"%*.*lf",1,0) $+ "f";
  
  i = wmax;
  if RLG_wstar==0;   i = binom[nfirms+2+kmax,kmax+3];  endif;    @ RLG_wstar==0 extended kmax by one in profit.g @
  load profit[i,nfirms]=^filename;

  if nfirms<maxfirms_qencode;
    @ check if can load  dtable etable1 and mask  since they take LONG time to build if nfirms high @
    filename = "bmarkov.dtable." $+ ftos(nfirms,"%*.*lf",1,0) $+ "f." $+ ftos(kmax,"%*.*lf",1,0) $+ "k";
    if (not filesa(filename) $== "");
      oneton = seqa(1,1,nfirms);          @ these variables also created in the code below that builds  dtable etable1 and mask @
      multfac1 = (kmax+1)^(oneton-1);
      two_n = 2^(nfirms-1);
      
      load dtable[nfirms,wmax]=^filename;
      
      filename = "bmarkov.etable1." $+ ftos(nfirms,"%*.*lf",1,0) $+ "f." $+ ftos(kmax,"%*.*lf",1,0) $+ "k";
      load etable1[(kmax+1)^nfirms,1]=^filename;
      
      filename = "bmarkov.mask." $+ ftos(nfirms,"%*.*lf",1,0) $+ "f." $+ ftos(kmax,"%*.*lf",1,0) $+ "k";
      load mask[nfirms-1,two_n]=^filename;
      
      "Loaded dtable etable1 and mask instead of building them...";
    else;   @ rebuild dtable etable1 and mask @
      
      two_n = 2^(nfirms-1);
      clear dtable;
      if nfirms > 1;
	oneton = seqa(1,1,nfirms);
	@ Build a mask of all binary numbers from 0 to two_n - 1 @
	mask = zeros(nfirms-1,two_n);
	for i (1, two_n, 1);                 @ was:  i=0;  do while i < two_n; ... i= i+1; endo; @
	clear msk; @ Build one mask @
	j=2; k=i-1;
	do while j <= two_n;
	  if k % j == 0;
	    msk = 0|msk;
	  else; k = k - j/2; msk = 1|msk;
	  endif;
	  j = j*2;
	endo;
	mask[.,i] = msk[1:nfirms-1];
	endfor;
	print "Mask done";
      endif;
      
      @ Make a table for quick decoding @
      dtable = zeros(nfirms,wmax);
      for i (1,wmax,1);   dtable[.,i] = decode(i);   endfor;
      
      @ Make a table for quick encoding @
      @ Fill in multfac1, multfac2, for quick encoding @
      if nfirms <= encfirm;
	multfac1 = (kmax+1)^(oneton-1);
	@ Encode all numbers from 1 to kmax^nfirms @
	etable1 = zeros((kmax+1)^nfirms,1);
	for i (1, rows(etable1), 1);     @ was:  i=0;  do while i < rows(etable1); ... i= i+1; endo; @
	  clear msk; @ Build one mask @
	  @ j=kmax+1; @
	  k = i-1;
	  if i%100000==0;  print "etable1 row " i " of " rows(etable1);  endif;
	  @ do while j <= rows(etable1); @       @ loop will execute nfirms times @
	  for t0 (1,nfirms,1);
	    j = (kmax+1)^t0;
	    msk = ((k % j)*(kmax+1)/j)|msk;
	    k = k - (k % j);
	    @ j = j*(kmax+1);  endo; @
	  endfor;
	  etable1[i] = encode2(rev(sortc(msk[1:nfirms],1)));
	endfor;
	print "etable1 done";
	if nfirms>1;     @ save file since creating dtable etable1 and mask are time-consuming @
	  /* write output */
	  "Saving dtable  etable1  and  mask ";
	  format /m1 /rds 1,0;
	  screen off;
	  filename = "bmarkov.dtable." $+ ftos(nfirms,"%*.*lf",1,0) $+ "f." $+ ftos(kmax,"%*.*lf",1,0) $+ "k";
	  output file = ^filename reset;	print dtable;   output off;

	  filename = "bmarkov.etable1." $+ ftos(nfirms,"%*.*lf",1,0) $+ "f." $+ ftos(kmax,"%*.*lf",1,0) $+ "k";
	  output file = ^filename reset;	print etable1;   output off;
	  
	  filename = "bmarkov.mask." $+ ftos(nfirms,"%*.*lf",1,0) $+ "f." $+ ftos(kmax,"%*.*lf",1,0) $+ "k";
	  output file = ^filename reset;	print mask;   output off;
	  
	  output file = ^filenam1;
	  screen on;
	endif;  @ save  dtable etable1 and mask  files @
	
      else;        @ nfirms > encfirm.   NO LONGER EVER HERE since RLG set encfirm = nfirms @
	multfac1 = 0|multfac1;
	multfac2 = (kmax+1)^(oneton[1:nfirms-encfirm]-1)|zeros(encfirm,1);
	/*
	"Multfac1 is " multfac1;
	"Multfac2 is " multfac2;
	*/
	etable2 = zeros((kmax+1)^(nfirms-encfirm),1);
	for i (1, rows(etable2), 1);     @ was:  i=0;  do while i < rows(etable2); ... i= i+1; endo; @
 	  clear msk; @ Build one mask @
	  j=kmax+1; k=i-1;
	  do while j <= rows(etable2);
	    msk = ((k % j)*(kmax+1)/j)|msk;
	    k = k - (k % j);
	    j = j*(kmax+1);
	  endo;
	  etable2[i] = encode2(rev(sortc((msk[1:nfirms-encfirm]|zeros(encfirm,1)),1)))-1;
	endfor;
	print "etable2 done";
      endif;
    endif;    @ load or rebuild  dtable etable1 and mask @
  else;       @ nfirms > maxfirms_qencode  --> build mask only @
    two_n = 2^(nfirms-1);
    clear dtable;
    oneton = seqa(1,1,nfirms);
    @ Build a mask of all binary numbers from 0 to two_n - 1 @
    mask = zeros(nfirms-1,two_n);
    for i (1, two_n, 1);                 @ was:  i=0;  do while i < two_n; ... i= i+1; endo; @
    clear msk; @ Build one mask @
    j=2; k=i-1;
    do while j <= two_n;
      if k % j == 0;
	msk = 0|msk;
      else; k = k - j/2; msk = 1|msk;
      endif;
      j = j*2;
    endo;
    mask[.,i] = msk[1:nfirms-1];
    endfor;
    print "Mask done.  Not building etable1 for qencode since nfirms > maxnfirms_qencode";
    @ Make a table for quick decoding @
    dtable = zeros(nfirms,wmax);
    for i (1,wmax,1);   dtable[.,i] = decode(i);   endfor;
  endif;

  InitializeZeros:
  
  @ Update values, or define starting values.  RLG moved next block from update() to here @
  oldx     = zeros(wmax,nfirms); 
  oldexit  =  ones(wmax,nfirms);    @ oldexit init as 1's since 1 is the newexit for inactive firms @
  oldvalue = zeros(wmax,nfirms); 
  oldentry = zeros(wmax,1);      
  v_entry  = zeros(wmax,1);

  if rows(RLG_x)==rows(oldx);
    "initializing with RLG_* values from previous model...";
    oldx     =  RLG_x     ;
    oldexit  =  RLG_exit  ;
    oldvalue =  RLG_value ;
    oldentry =  RLG_entry ;
    v_entry  =  RLG_v_entry  ;
  endif;
  
  foc          = zeros(wmax,2*nfirms+1);   @ these 4 globals simplify multithreading with foc @
  rflag        = zeros(wmax,1);
  x_exit_entry = zeros(wmax,2*nfirms+1);
  foc_multi_eq = zeros(wmax,1);            @ # of iterations for which contract_w() finds a 2nd equilibrium, in the chk_foc_multi_eq random restarts to check for multi_eq @

  if delta==0 AND RLG_wstar<=0;          @ Be sure simulations start with a firm at kmax @
    do_w = selif(seqa(1,1,wmax), dtable[1,.]'.<kmax);	@ then will ALWAYS have a firm at kmax, so skip states with 1st firm below kmax @ 
    oldx[do_w,.]     = oldx[do_w,.]-1;
    oldexit[do_w,.]  = oldexit[do_w,.]-1;
    oldvalue[do_w,.] = oldvalue[do_w,.]-1;
    oldentry[do_w,.] = oldentry[do_w,.]-1;
    do_w = selif(seqa(1,1,wmax), dtable[1,.]'.==kmax);
    "Only solving for values and policies over states with firm 1 at the frontier, since delta = 0";
  else;
    do_w = seqa(1,1,wmax);
  endif;
  wmax2 = rows(do_w);

  do_w = sortc( do_w~rndu(wmax2,1), 2);   @ randomize the order states are to be processed  @
  do_w = do_w[.,1];                       @ so threads have similarly time-consuming states @

  newx     = oldx;      @ initialize new here, after set -1 for ignored states @
  newexit  = oldexit;
  newvalue = oldvalue;
  isentry  = oldentry;
  
  @ update(); @         @ skipping update() which uses values for nfirms-1 equilibrium as starting values @

  maxiter = 400/(1-RLG_damp);
 
  StartIterating:

  format /m1 /re;
  "Contraction...over " padr(wmax2,1,0) " states of " padr(wmax,1,0) " in dtable with RLG_foc = " padr(RLG_foc,1,0)  "  min frontier firm innov = "  aeff[rows(aeff)]*RLG_minx1/(1+aeff[rows(aeff)]*RLG_minx1) "  maxiter= "  padr(maxiter,1,0)  "  tol= " tol "  foc_tol= " foc_tol  " ...";
  if chk_foc_multi_eq;  "In each iteration and at each state, restarting search for FOC solution at "  padr(chk_foc_multi_eq,1,0) " random policies to investigate presence of multiple equilibria.";  endif;
  format /m1 /rd 8,4;
  norm = tol + 99;  iRLG=1;   normx = norm;  normexit= norm;  normentry = norm;
  avgnorm = norm;  prevnorm = norm+1;  prevavgnorm = norm+1;  t0 = hsec;   @ hsec= hundredths of seconds since midnight @
  norm100 = norm+9;

  do while (norm > tol) AND (iRLG < maxiter OR (norm-prevnorm<.01 AND iRLG<400) OR (avgnorm-prevavgnorm<1e-7 AND iRLG<400));  /*  AND (avgnorm > 0.0001*tol);   <--  REMOVED by RLG the avgnorm criteria */    
    prevnorm = norm;

    if iRLG%100==0;  norm100= (norm+prevnorm)/2;     @ attempt to detect non-converg to avoid wasted effort.  avg avoids getting peak norm of a cycle. @
    else;
      if iRLG>100 AND iRLG-100*floor(iRLG/100)>31 AND norm>norm100;           @ if iRLG > 100 and higher norm than 31 ago, likely not converging @
	"aborting contraction since current norm is worse than prev 100th iteration, suggesting convergence issues, ...";
	break;
      endif;
    endif;
    
    /*
    if maxc(maxc(oldexit))>1;
      xxx = indexcat(maxc(oldexit').>1,1);
      print rows(xxx) " states have oldexit>1  at iteration " iRLG+1;
      "Element no./Probability Exit (wmax[nfirms-1]): "; dtable[.,xxx]'~oldexit[xxx,.]; "";
    endif;
    
    xxx = sumc( (oldexit[.,1:nfirms-1] .> oldexit[.,2:nfirms] + .01)' .AND   (oldvalue[.,1:nfirms-1]-profit[.,1:nfirms-1] .> oldvalue[.,2:nfirms]-profit[.,2:nfirms] + .01)' ).>0 ;
    if sumc(xxx)>0;             @ -profit[] in prev line since exit is based on continuation values (ie, profits are earned regardless of entry) @
      aaaaa = indexcat(xxx,1);
      format /m1 /rd 6,2;
      print "HEY: exit goes down when value goes up as w increases for " rows(aaaaa) " states:   w ~ value ~ exit ~ x";
      dtable[.,aaaaa]'~-9*ones(rows(aaaaa),1)~oldvalue[aaaaa,.]~-9*ones(rows(aaaaa),1)~oldexit[aaaaa,.]~-9*ones(rows(aaaaa),1)~oldx[aaaaa,.];
    endif;
    */

    if BatchJob== -1 AND iRLG>200;  RLG_damp = .99;  maxiter = 200000; endif;

    contract();
    
    abs_error = abs(oldvalue - newvalue);
    norm      = maxc(maxc(abs_error));
    avgnorm   = meanc(meanc(abs_error));
    normx     = maxc(maxc(oldx-newx));
    normexit  = maxc(maxc(oldexit-newexit));
    normentry = maxc(oldentry-isentry);
   /* "Sup norm: " padr(norm,8,4) "; Mean norm: " padr(avgnorm,8,4) "    \r";; */

    iRLG = iRLG+1;
    if (iRLG %10)==0 OR nfirms>9 ;
      "iteration: " padr(iRLG,6,0) " avg~sec " padr((hsec-t0)/100/(1+9*(nfirms<=9)),6,1) "  norm~avgnorm: " padr(norm,8,5) " ~ " padr(avgnorm,11,8) "  norm x~exit~entry: " padr(normx,10,6) padr(normexit,10,6) padr(normentry,10,6) "  maxV " padr(maxc(newvalue[.,1]),7,2) "\n" ;;  print /flush;;
      t0 = hsec;
    endif;

    @ if iRLG<4;  oldvalue~newvalue;   endif; @

    if 0 AND iRLG > 100;
      format /m1 /rd 8,2;  xxx = ones(rows(isentry),1)-1.11;
      "dtable ~ values (old then new then diff)";  print dtable'~oldvalue~newvalue~xxx~newvalue-oldvalue;
      "dtable ~ x (old then new then diff)"; print dtable'~oldx~newx~xxx~newx-oldx;
      "dtable ~ exit ~ entry (old then new then diff)"; print dtable'~oldexit~newexit~xxx~newexit-oldexit~xxx~oldentry~isentry~xxx~isentry-oldentry;
    endif;
    
    if nfirms<7 AND norm > .001 AND iRLG>100;

      screen off;
      
      normind1 = maxindc(maxc(abs_error'));    @ oldvalue is wmax by nfirms.  so maxc(abs_error') is wmax vector and normind1 is state with highest error @
      normind2 = maxindc(maxc(abs_error));     @ maxc(abs_error) is firm vector of highest error across states, so normind2 is firm with highest error @
      
      if BatchJob==0   AND 0  ;  normind1 = qencode(7|7|5);  normind2 = 1;  endif;
    
      normcode = dtable[.,normind1]';
      format /rd 7,4;
      
      aaaaa = abs( newx-oldx );           xind1 = maxindc(maxc(aaaaa'));   xcode = dtable[.,xind1]';
      bbbbb = abs( newexit-oldexit );     eind1 = maxindc(maxc(bbbbb'));   ecode = dtable[.,eind1]';
      /*
      if sumc(sumc(aaaaa))==1;
	RLG1 = maxindc(maxc(aaaaa'));
	RLG2 = maxindc(maxc(aaaaa));
      else;
	RLG1 = 1; RLG2 = 0;
      endif;
      */
      @ "rows~cols oldvalue:"  rows(oldvalue)~cols(oldvalue) "rows~cols oldx:"  rows(oldx)~cols(oldx); @
      
      "Norm~Avgnorm: " padr(norm,8,5) " " padr(avgnorm,11,8) "  "
      "Max Norm: firm " padr(normind2,1,0) " at state " normcode "; Old~New value: "
      oldvalue[normind1,.] " ~ "
      newvalue[normind1,.] "   Old~New x: "
      oldx[normind1,.]     " ~ "
      newx[normind1,.]     "   Old~New exit: "
      oldexit[normind1,.]     " ~ "
      newexit[normind1,.]  "   Entry: "
      isentry[normind1] v_entry[normind1]
      "   Old~New x at "    xcode " is "    oldx[xind1,.] "~"    newx[xind1,.]
      "   Old~New exit at " ecode " is " oldexit[eind1,.] "~" newexit[eind1,.]
      @ " sole exit change at firm " padr(RLG2,1,0) " state " dtable[.,RLG1]'  @  ;
      
      screen on;
    endif;
    
    if RLG_damp OR (  0 AND  iRLG > 100 AND norm>prevnorm) ;
      @ i = RLG_damp/maxc(iRLG|(iRLG/100));  @              @ can decrease dampening as iRLG increases @
      i = maxc(RLG_damp|(.5*(norm>prevnorm))*rndu(1,1));    @ random weight on oldx, with max at RLG_damp or .5 if norm>prevnorm @
      if norm>prevnorm AND RLG_damp==0;  "norm > prevnorm, dampening weight on old = " i " at iter= " iRLG;  endif;
      oldx     =     oldx*i + (1-i)*newx;
      oldexit  =  oldexit*i + (1-i)*newexit;
      oldvalue = oldvalue*i + (1-i)*newvalue;
      oldentry = oldentry*i + (1-i)*isentry;
    else;
      oldx = newx;  oldvalue = newvalue;  oldexit = newexit;   oldentry = isentry;
    endif;
    
    format /rd 9,6;   output off;   output on;

    /*
    i = seqa( rows(oldvalue)-9,1,10);
    "iter = "  iRLG  "  qencode_index  value   x   exit ";
    format /m1 /rds 8,3;  i~oldvalue[i,.]~oldx[i,.]~oldexit[i,.];
    */
    
    @ if RLG_foc==0 AND iRLG > 50; " ";  "Experimental hardcoded break of non-foc effort to switch to foc ";  " ";    break;  endif; @
  
  endo;  @ while norm > tol @

  
  format /rd 7,3;
  d2 = date;
  print compact(nfirms) " firm(s) scenario completed; execution time = ";;
  print padr(ethsec(d1,d2)/6000,8,2) " minutes, #iterations= " padr(iRLG,6,0) ;; format /m1 /re 6,1; " tol= " tol  "  norm=" norm "  prevnorm~avgnorm=" prevnorm prevavgnorm "  norm x~exit~entry= " normx normexit normentry;  " ";

  @ isentry[99999999]; @    @ deliberate abort @
  
    
  if maxc( foc_multi_eq ) > 0;
    d2 = selif( seqa(1,1,wmax), foc_multi_eq .> 0);          @ states with multiple eq @
    if rows(d2)>30;
      d2 = rev(sortc(d2,1));
      d2 = d2[1:30];
    endif;
    kk = 1.111 + zeros(rows(d2),1);
    format /m1 /rd 4,0;  "Multiple FOC solutions were obtained most frequently at the following states:   # multieq  ~  state  ~  entry  ~  investment  ~  exit";
    format /m1 /rd 8,4;  print foc_multi_eq[d2]~kk~dtable[d2,.]~kk~isentry[d2]~kk~newx[d2,.]~kk~newexit[d2,.];
  endif;
  
  if norm>tol;   " ";     @ if BatchJob and RLG_minx1>0, do NOT try FOC -- unlikely to help since already tried FOC with minx1=0 and non-FOC with minx1>0 @
    if RLG_foc==0 AND pp[_RLG_FOC]==0 AND (BatchJob==0 OR RLG_minx1==0); 
      RLG_foc = 1;
      maxiter = 200;  "non-FOC failed to converge... retrying from current with FOC method...";   goto StartIterating;
    endif;
    if  0 AND  BatchJob>0 AND RLG_minx1==0;                  @ since this overwrites the config pp[] only do this when in Batch mode @
      "non-FOC and possibly FOC failed to converge ... retrying non-FOC from current after setting pp[_RLG_MIN_INNOV1]= .01...";
      pp[_RLG_MIN_INNOV1]= .01;
      RLG_minx1 = pp[_RLG_MIN_INNOV1]/((1-pp[_RLG_MIN_INNOV1])*aeff[rlnfirms+1]);  @ ax = (1+ax) m  --> ax(1-m) = m -->  x = m/a(1-m) @ 
      RLG_foc = 0;  maxiter = 200;
      goto StartIterating;
    endif;
    
    if RLG_foc >0;    "FOC method failed to converge...";      endif;
    if RLG_foc==0;    "non-FOC method failed to converge...";  endif;
  endif;

  @ Now find if there is any investment, entry at the highest level @

  probw = reshape(aeff[dtable'+1],wmax,nfirms).*newx;     probw = probw./(1+probw);

  d2 = "   ";  for kk(nfirms,9,1);  d2 =  "0 " $+ d2;  endfor;   d2 = "  " $+ d2;   @ build string of 0 for filler columns in output below @
  " ";
  "output below corresponds to 10 firm case since matlab code plotting policy functions across models uses column counts of 10 firm industry";
  " ";
  " V and policies at select states:     w   ~  V  ~  entry  ~  innovation  ~  exit ";
  w = ones(nfirms,1);
  format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  " ";
  w = kmax*ones(nfirms,1);     format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  @LINE 522@ 
  for kk(nfirms,2,-1);
    w[kk]=0;    format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  endfor;
  " ";
  w = entry_k*ones(nfirms,1);  w[1] = kmax;    format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  for kk(nfirms,2,-1);
    w[kk]=0;   format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  endfor;
  " ";
  w = zeros(nfirms,1);  w[1] = kmax;
  for kk(kmax,1,-1);
    w[2]=kk;   format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  endfor;
  " ";
  w = ones(nfirms,1);  w[1] = kmax;  w[nfirms] = 0;
  for kk(kmax,1,-1);
    w[2]=kk;   format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  endfor;
  " ";
  w = ones(nfirms,1);  w[1] = kmax;
  for kk(kmax,1,-1);
    w[2]=kk;   format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  endfor;
  " ";
  if nfirms>2;
    w = zeros(nfirms,1);  w[1] = kmax;  w[2] = kmax;
    for kk(kmax,1,-1);
      w[3]=kk;   format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
    endfor;
  endif;
  " ";
  if nfirms>2;
    w = zeros(nfirms,1);  w[1] = kmax;  w[2] = kmax-1;
    for kk(kmax-1,1,-1);
      w[3]=kk;    format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
    endfor;
  endif;
  " ";
  w = rev(seqa(1,1,nfirms));  w = minc(w'|kmax*ones(1,nfirms));
  w[1]=kmax;  format /m1 /rd 3,0;  w';; format /m1 /rd 8,4; d2 newvalue[qencode(w),.] d2  isentry[qencode(w)] "   " probw[qencode(w),.]  d2 newexit[qencode(w),.] d2;
  
  " ";
  w= sumc(newexit'.>0).==nfirms;
  if sumc(w)>1;
    iRLG = selif(seqa(1,1,wmax),w);
    "All firms exit at the following states: w ~ V ~ x ~ entry ~ exit";
    format /m1 /rds 5,2;  dtable[.,iRLG]'~newvalue[iRLG,.]~-9.99*ones(rows(iRLG),1)~newx[iRLG,.]~-9.99*ones(rows(iRLG),1)~isentry[iRLG]~-9.99*ones(rows(iRLG),1)~newexit[iRLG,.];
  endif;
  " ";

  if RLG_wstar>0;
    w = zeros(nfirms,1);  w[1] = kmax-1;
    kk = newx[qencode(w):(qencode(kmax|zeros(nfirms-1,1))-1),1]; @ firm 1 investment over all states with firm 1 at kmax-1 @
    iRLG = maxc(kk);
    if iRLG > 0;
      "Warning: Positive investment of  " iRLG "  recorded at 2nd highest efficiency level, yielding success prob = " aeff[kmax]*iRLG/(1+aeff[kmax]*iRLG);
      "Please consider increasing the maximum efficiency level (kmax) since using PM wstar to bound statespace.";
    endif;
  else;
    w = zeros(nfirms,1);  w[1] = kmax;
    kk = newx[qencode(w):wmax,1]; @ firm 1 investment over all states with firm 1 at frontier @
    iRLG = maxc(kk);
    "Max investment of  " iRLG "  recorded at highest efficiency level, yielding success prob = " aeff[kmax+1]*iRLG/(1+aeff[kmax+1]*iRLG) " (which is okay since not using wstar to bound statespace)";
    iRLG = kk-RLG_minx1.<1e-15;
    if sumc(iRLG)>0;
      iRLG = selif(seqa(qencode(w),1,rows(kk)),iRLG);
      "Minimum investment of " RLG_minx1 " by frontier firm at " rows(iRLG) " of " rows(kk) " states with leader at frontier:  w ~ V ~ x ~ entry ~ exit";
      format /m1 /rds 6,2;  dtable[.,iRLG]'~newvalue[iRLG,.]~-9.99*ones(rows(iRLG),1)~newx[iRLG,.]~-9.99*ones(rows(iRLG),1)~isentry[iRLG]~-9.99*ones(rows(iRLG),1)~newexit[iRLG,.];
    endif;
  endif;
  " ";
  w = sumc( (dtable[.,do_w].==1 .AND newexit[do_w,.]'.<1-1e-8));      @ states with firm not exiting for sure at lowest rung @
  if sumc(w)>0;
    kk = selif(do_w,w.>0);
    "Firms on lowest rung w/ Pr(exit)<1 at " padr(rows(kk),1,0) " states.  They get free bump if RLG_w0_exit=0, else exit w/ avg scrap when fall off ladder.  Hardcoded RLG_w0_exit is " padr(RLG_w0_exit,1,0);
    "First few of such states:  value  ~  invest  ~  entry  ~  exit";
    kk = kk[1:minc(10|rows(kk))];  w = -9.99*ones(rows(kk),1);
    format /m1 /rds 6,2;  dtable[.,kk]'~newvalue[kk,.]~w~newx[kk,.]~w~isentry[kk]~w~newexit[kk,.];
  endif;
  " ";
  w = sumc( (dtable[.,do_w].==1 .AND newexit[do_w,.]'.<1-1e-8 .AND newx[do_w,.]'.>1e-8));      @ indicator for states with firm at lowest rung with Pr(exit)<1 and investing x>0 @
  if sumc(w)>0;
    kk = selif(do_w,w.>0);
    "Firms on lowest rung w/  invest> 0 at " padr(rows(kk),1,0) " states.  They get free bump if RLG_w0_exit=0, else exit w/ avg scrap when fall off ladder.  Hardcoded RLG_w0_exit is " padr(RLG_w0_exit,1,0);
    "First few of such states:  value  ~  invest  ~  entry  ~  exit";
    kk = kk[1:minc(10|rows(kk))];  w = -9.99*ones(rows(kk),1);
    format /m1 /rds 6,2;  dtable[.,kk]'~newvalue[kk,.]~w~newx[kk,.]~w~isentry[kk]~w~newexit[kk,.];
  endif;
  " ";
  
  @ Store data in file for inspection @
  screen off;
  "Element no./ Value Function (wmax x nfirms): "; dtable'~newvalue; "";
  "Element no./Investment (wmax x nfirms): "; dtable'~newx; "";
  "Element no./Probabilities of p rising (wmax x nfirms) ~ Investment ~ Efficiency: ";  
  print   dtable'~probw~newx~reshape(aeff[dtable'+1],wmax,nfirms); "";
  if phiH /= 0;
    "Element no./Probability Exit (wmax[nfirms-1]): "; dtable'~newexit; "";
  endif;
  "Element no./Probability ~ Value of entry (wmax[nfirms-1]): "; dtable'~isentry~v_entry; "";
  "Element no./Probability ~ Value of entry (wmax[nfirms-1]): only states with room for entrant"; selif(dtable'~isentry~v_entry, dtable[nfirms,.]'.==0); "";
  screen on;
  output off;
  @ Store data in file for comparative statics program to read @
  print "Generating data for comparative statics --> " filenam2;
  format /m1 /rds 16,12;
  output file = ^filenam2 reset;
  screen off;
  newvalue; ""; 
  newx; "";
  probw; "";
  if phiH /= 0;
    newexit; "";
  endif;
  isentry;
  v_entry;
  output off;
  screen on;
  output file = ^filenam1 on;
  nfirms = nfirms+1;

  @ XP = (w./(1+w));  output file=XX ;  output on;  print  (kmax*ones(rows(newx),1))~dtable'~newx~XP~XP./maxc(XP);  output off; @

  format /m1 /rd 5,2;
  if delta==0 AND RLG_wstar<=0;
    iRLG = sumc(newx')-RLG_minx1*sumc(dtable.==kmax).<1e-15 .AND isentry.==0 .AND sumc(newexit').==0;
    @ if sumc(iRLG)>0;  " "; "absorbing states with MIN INVESTMENT, NO ENTRY, NO EXIT: state ~ value ~ exit ~ entry";  selif( dtable'~newvalue~newexit~isentry, iRLG ); " ";  @
    if sumc(iRLG)>0;  " "; "absorbing states with MIN INVESTMENT, NO ENTRY, NO EXIT:";  selif( dtable', iRLG ); " ";  
    endif;
  endif;
  
endo;   @ while nfirms <= rlnfirms.   nfirms started at stfirm which itself is set to  stfirm = pp[_START_FIRMS];  in init.h  @

if norm <= tol;             @ save in global variables to use as initial values for next model @
  RLG_x        = oldx    ;
  RLG_exit     = oldexit ;
  RLG_value    = oldvalue;
  RLG_entry    = oldentry;
  RLG_v_entry  = v_entry ;
endif;

output off;
screen on;
pp[_EQL_DONE]= (norm <= tol);

@ print "Press any key to return to the main menu"; wait; @

save ^configfile=pp;

if BatchJob==0;  run pmgshell.g;  endif;


proc (1) = get_foc(x_exit_ent, w, doo);
  @ obtain foc for state w at candidate policy x_exit_ent  @
  @ x_exit_ent   are policies being changed to find root   @
  @ doo          are indices of policies being changed     @
  @ using parameters, not globals, since multithreaded     @
  @ requires each thread to have its own copy.             @
  @                                                        @
  @ modified optimize() and chkentry() to accept candidate policies @
  @ as inputs instead of using global oldx oldexit isentry @
  local foc_w, x_exit_entry_w, locw, i;

  x_exit_entry_w = ( oldx[w,.]~oldexit[w,.]~oldentry[w] )';
	
  if rows(doo)>0;   x_exit_entry_w[doo] = x_exit_ent;   endif;   @ possible doo = {} @

  locw  = dtable[.,w]';
  for i (2,nfirms,1);
    if locw[i]==locw[i-1];                                     @ impose symmetry  @
      x_exit_entry_w[i]        = x_exit_entry_w[i-1];          @ same  x   policy @
      x_exit_entry_w[i+nfirms] = x_exit_entry_w[i-1+nfirms];   @ same exit policy @
    endif;
  endfor;
    
  /*
  if locw[1]==locw[2];
    format /m1 /re 10,1;   print "before chkentry(): x " x_exit_entry_w[1:nfirms]'  "  exit " x_exit_entry_w[nfirms+1:2*nfirms]'  "  entry " x_exit_entry_w[2*nfirms+1] ;;
    format /m1 /rd 1,0 ;   print " locw " locw  " doo: " doo';
  endif;
  */
  chkentry_w(x_exit_entry_w, w);                                             @ implicitly returns v_entry[w] and isentry[w] @

  {newx[w,.], newvalue[w,.], newexit[w,.]}= optimize_w(x_exit_entry_w,w);    @ THIS is the global updating of policies since the last call for each state will yield foc = 0 and hence be equilibrium for the state this iteration @
  
  @ format /m1 /re 10,1; print "x_exit_entry_w"  x_exit_entry_w' "  newx " newx[w,.] "  exit " newexit[w,.] "  entry " isentry[w];  @
  
  foc_w = x_exit_entry_w-(newx[w,.]~newexit[w,.]~isentry[w])';

  if rows(doo)==0;  retp( 1e-99 );
  else;             retp(foc_w[doo]);
  endif;
endp;
  

proc (0) = contract_w(w1);
  @ This procedure created so that one line of code needs to be called    @
  @ for each state in contract(), which facilitates clean multithreading. @
  @ See documentation for contract() for more details @
  @ global:  x_exit_entry  foc  rflag  (among others) @
  local t0, do1, locw, dosymm, x_exit_entry_val, foc_w;        @ w2, w3, x9, y9 are only used in checking code usually turned off @

  locw  = dtable[.,w1];
  
  if RLG_foc<=0;    @ do not use Newton.  Call get_foc() just once to get best-response as in standard PM @
    
    foc[w1,.] = get_foc(x_exit_entry[w1,.]',w1,oneto2nfirms1)';                               @  get_foc() updates GLOBAL x_exit_entry and GLOBAL newx newexit isentry newval @
    
  else;             @ use Newton to solve foc at this state @

    dosymm = selif( seqa(1,1,nfirms), locw./=(-1|locw[1:nfirms-1]) .AND locw.>0 );            @ dosymm is list of policies to change after imposing symmetry (i.e., firms with same w take same actions) @
    dosymm = dosymm|(nfirms+dosymm)|2*nfirms+1;                                               @ newton search is restricted to params in dosymm to impose symmetry @
   
    do1 = selif( dosymm, x_exit_entry[w1,dosymm]'./=0 .AND x_exit_entry[w1,dosymm]'./=1 );    @ most 0,1 policies persist, so first do newton() holding them fixed.  Will check their optimality later @

    if ismiss(do1);  do1 = 1;    endif;                                                       @ always search at least over the leader's x @
  
    @ format /m1 /rd 2,0;  "locw " locw'  "  dosymm " dosymm' "  do1 " do1'; @
    
    if rows(do1)>0;
      x_exit_entry[w1,do1] = newton3( x_exit_entry[w1,do1]', &get_foc, w1, do1, foc_tol)';    @ updates GLOBAL x_exit_entry and GLOBAL newx newexit isentry newval @
    endif;
      
    @ call with optimal value to ensure global newx newexit isentry newvalue properly updated with correct values, and get foc for ALL params @
    
    foc[w1,.] = get_foc(x_exit_entry[w1,.]',w1,oneto2nfirms1)';    @ foc = 0 where locw = 0.  get_foc() imposes symmetry, so could use dosymm instead of oneto2nfirms1 @
    
    rflag[w1] = maxc(abs(foc[w1,.]'));
    
    @ dpolicy = ( oldx[w1,.]~oldexit[w1,.]~oldentry[w1] )' - ( newx[w1,.]~newexit[w1,.]~isentry[w1] )'; @   @ change in policies @
    
    /*  @ need to uncomment x9 y9 w2 w3 in local declarations @
    if 0 AND iRLG==33 AND w1==30;         @ newton failed, so get best-response functions (holding entry, exit fixed) to inspect @
      x9 = zeros(200,1);                  @ 1's BR x to 2's x @
      y9 = zeros(200,1);                  @ 2's BR x to 1's x @
      for w2 (1,200,1);
      x_exit_entry[w1,2] = w2/100;   do1 = 1;  w3 = newton3( oldexit[w1,do1], &get_foc, w1, do1, 1e-4);  x9[w2] = w3;
      x_exit_entry[w1,1] = w2/100;   do1 = 2;  w3 = newton3( oldexit[w1,do1], &get_foc, w1, do1, 1e-4);  y9[w2] = w3;
      endfor;
      format /m1 /rd 8,4;  print x9~y9;
      x9[999];  @ aborts @
    endif;
    */
    
    if rows(do1)<rows(dosymm) AND rflag[w1,.]>foc_tol;              @ redo newton3() with do1 = dosymm if do1 did not include all of dosymm AND  @ 
      do1 = dosymm;                                                 @ some foc > foc_tol (either do1 did not converge OR some fixed policies need to be in do1) @
      x_exit_entry[w1,do1]  = newton3( x_exit_entry[w1,do1]', &get_foc, w1, do1, foc_tol)';
      foc[w1,.] = get_foc(x_exit_entry[w1,.]',w1,oneto2nfirms1)';
      rflag[w1] = maxc(abs(foc[w1,.]'));
    endif;
    
    if rflag[w1]>1e-2 ;
      format /m1 /rd 4,0;  print "foc: iter~s " iRLG~w1 " w1 " dtable[.,w1]';; 
      format /m1 /rd 7,3;  print newx[w1,.]~-1.111~newexit[w1,.]~-1.111~isentry[w1]~log(maxc(abs(foc[w1,.]')));;
      format /m1 /rd 3,0;  print do1';
      @ format /m1 /rd 4,0;   print "x_exit_entry*  " iRLG~w1 " w " dtable[.,w1]';; format /m1 /rd 7,3; print x_exit_entry_X[1:nfirms]'~-9~x_exit_entry_X[1+nfirms:2*nfirms]'~-9~x_exit_entry_X[2*nfirms+1]~log(maxc(abs(foc)));;  format /m1 /rd 3,0; print do1'; @
    endif;

    if chk_foc_multi_eq;         
      t0 = 0;  do1 = dosymm;  x_exit_entry_val = x_exit_entry[w1,.]~newvalue[w1,.];
      do while t0 < chk_foc_multi_eq;  
	t0 = t0+1;
	x_exit_entry[w1,do1]  = newton3( LB[do1]+rndu(rows(do1),1).*(minc(UB[do1]'|10*ones(1,rows(do1)))-LB[do1]), &get_foc, w1, do1, foc_tol)';    @ random starting policies @
	foc_w = get_foc(x_exit_entry[w1,.]',w1,oneto2nfirms1)';
	if maxc(abs(foc_w))<foc_tol AND maxc( abs( x_exit_entry[w1,do1] - x_exit_entry_val[do1] ))>.001;
	  foc_multi_eq[w1] = 1+foc_multi_eq[w1];
	  break;  @ multi-eq found @
	endif;
      endo;
      x_exit_entry[w1,.] = x_exit_entry_val[1:2*nfirms+1];               @ restore initial equilibrium @
      newvalue[w1,.]     = x_exit_entry_val[2*nfirms+2:3*nfirms+1];
      newx[w1,.]         = x_exit_entry_val[1:nfirms];
      newexit[w1,.]      = x_exit_entry_val[1+nfirms:2*nfirms];
      isentry[w1,.]      = x_exit_entry_val[2*nfirms+1];
    endif;
    
  endif;    @ RLG_foc @
endp;


proc (0) = contract();
  @ This procedure does one iterative step on investment and the value fn @
  @ Implicit parameters are  oldx, oldexit, oldvalue (passed in)          @
  @                     and  newx, newexit, newvalue (returned)           @
  local t0, w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,  w, wmaxT, Nthreads;
 
  @ First: check for which values of w_s would a firm want to enter @
  @        Entry decision is made at very beginning of each period  @
  @        BEFORE exit is implemented (done in optimize),           @
  @        BUT I think isentry should account for expected exit.    @
  @        Is this an error in this PM code (C also) ?              @
  @        chkentry() calls calcval() which integrates over invest  @
  @        outcomes, but does not implement exit, which is done in  @
  @        optimize() prior to calling calcval() when getting newx  @
  @        Of course, when I modify code to allow for random scrap  @
  @        calcval() will integrate over exit as well as invest.    @
  
  @ OLD call: {newx[w,.], newvalue[w,.], newexit[w,.]}= optimize(w) @
  @ Now the global new values (newx, newvalue, newexit, isentry)    @
  @ are updated in  get_foc()  which is called by  optimize_w()     @
  
  @ chkentry(); @     @ Turn back on, and copy isentry to oldentry to restore PM updating of entry before investment/exit @

  @ t0 = hsec;  "starting contract()"; @

  foc          = zeros(wmax,2*nfirms+1);     @ these 3 are global to facilitate multithreaded code @
  rflag        = zeros(wmax,1);
  x_exit_entry = oldx~oldexit~oldentry ;

  Nthreads = 8;     /* comment/uncomment out ThreadBegin lines below to match value of Nthreads */
  
  if nfirms<3 OR Nthreads==1;
    
    for w1 (1,wmax2,1);   contract_w( do_w[w1] );  endfor;
  
  else;      
    
    wmaxT = floor(wmax2/Nthreads);    @ recall wmax2 = rows(do_w) @

    ThreadBegin;  for w (1+ (Nthreads-1)*wmaxT,wmax2,1);  contract_w( do_w[w] );  endfor;  ThreadEnd;    @ the residual thread @

    ThreadBegin;  for w1  (1,wmaxT,1);  contract_w( do_w[ w1           ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w2  (1,wmaxT,1);  contract_w( do_w[ w2+wmaxT     ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w3  (1,wmaxT,1);  contract_w( do_w[ w3+wmaxT*2   ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w4  (1,wmaxT,1);  contract_w( do_w[ w4+wmaxT*3   ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w5  (1,wmaxT,1);  contract_w( do_w[ w5+wmaxT*4   ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w6  (1,wmaxT,1);  contract_w( do_w[ w6+wmaxT*5   ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w7  (1,wmaxT,1);  contract_w( do_w[ w7+wmaxT*6   ] );  endfor;  ThreadEnd;
/*
    ThreadBegin;  for w8  (1,wmaxT,1);  contract_w( do_w[ w8+wmaxT*7   ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w9  (1,wmaxT,1);  contract_w( do_w[ w9+wmaxT*8   ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w10 (1,wmaxT,1);  contract_w( do_w[ w10+wmaxT*9  ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w11 (1,wmaxT,1);  contract_w( do_w[ w11+wmaxT*10 ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w12 (1,wmaxT,1);  contract_w( do_w[ w12+wmaxT*11 ] );  endfor;  ThreadEnd;
    ThreadBegin;  for w13 (1,wmaxT,1);  contract_w( do_w[ w13+wmaxT*12 ] );  endfor;  ThreadEnd;
*/
    ThreadJoin;
  endif;  @ use threads @
  
  @ " seconds in contract() = " (hsec-t0)/100; @

  if RLG_foc>0;  w1= maxc(rflag);  if w1>1e-5; format /m1 /rd 6,2;  "mean~max log10(foc) " log(meanc(rflag[do_w])~w1);  endif;  endif;

endp;     @ Implicit returned parameters: newx, newvalue, newexit @


/*
proc (0) = update();
@ This procedure takes the solved newx, newvalue matrix for the nfirms - 1 problem @
@ and puts them into the nfirms matrices oldx, oldexit, oldvalue, for use as starting values @
  local w,i,n,tuple;

  "This procedure has not been updated for use of FOC method.  Out of bounds index being used to abort...";  nfirms[2];
  
  if nfirms == 1;
    for i (1,wmax,1);
      oldvalue[i,.] = 1 + 0.1*i;
    endfor;
  else;
    for w (1,wmax,1);
      tuple = dtable[.,w];  @ qdecode(w); @
      nfirms = nfirms - 1;
      n = encode2(tuple[1:nfirms]);
      oldx[w,1:nfirms]     = newx[n,1:nfirms];
      oldexit[w,1:nfirms]  = newexit[n,1:nfirms];
      oldvalue[w,1:nfirms] = newvalue[n,1:nfirms];
      nfirms = nfirms + 1;
      tuple[nfirms-1] = tuple[nfirms];
      tuple[nfirms] = 0;
      oldvalue[w,nfirms] = oldvalue[encode2(tuple),nfirms-1];
      oldexit[w,nfirms]  = oldexit[encode2(tuple),nfirms-1];
      oldx[w,nfirms]     = oldx[encode2(tuple),nfirms-1];
      @ format /m1 /rds 3,0; "w " dtable[.,w]';;  format /m1 /rds 7,3; " v " oldvalue[w,.] " x " oldx[w,.] " exit " oldexit[w,.]; @
    endfor;
  endif;
  @ Implicit returned value: oldx, oldexit, oldvalue @
endp;
*/

proc (3) = optimize_w(x_exit_entry,w);
@ This procedure calculates optimal investment, and value fn., for a @
@ given industry structure w. Thus, a vector nfirms long of each is returned. @
@ Implicit parameters are oldx, oldvalue, isentry @

  local locw,locwx,locwe,  @ Decoded copies of omegas with entry (locwe) and w/o entry (locwx) @
        oval,ox,oexit,     @ Old local values @
        entered,           @ Indicates the probability of an entrant @
        v1,v2,             @ v1: value of investing; v2: value of not investing @
        i,j,p,r,tempv1,tempv2,temp,tempk1,tempk2,
        nval,nx,nexit,     @ Returned values of investment, value fn. @
	eexit,             @ copy of exit policy with 0 in last position for the potential entrant @
	locwk,doleap;      @ locw when entrant starts at kmax @

  locw  = dtable[.,w];  @ qdecode(w); @     @ profit(w,j) is used --> product market competition is BEFORE entry, exit, investment outcomes @
  locwx = locw;
  oval  = oldvalue[w,.]';
  ox    = x_exit_entry[1:nfirms];            @ non-foc: oldx[w,.]' @
  oexit = x_exit_entry[nfirms+1:2*nfirms];   @ non-foc: oldexit[w,.]' @
  nval  = zeros(nfirms,1);
  nx    = zeros(nfirms,1);
  nexit = zeros(nfirms,1);
  
  @ dop = 0; @
  
  @ if nfirms==3;  if locw[1]==7 and locw[2]==2 and locw[3]==2;  dop = 1;  endif;  endif; @
  
  if phiH==0;    @ FIXED scrap --> use old PM code.  Could change to use the exit mask in calcval() @
    @ Find out which firms want to exit @
    " HEY: probably should not run fixed scrap code anymore";
    i = (minc(oval) == phi)*(minindc(oval)-1) + (minc(oval) > phi)*nfirms;
    
    @ Replace efficiency levels of exitors with zero     @
    if i < nfirms;    locwx[i+1:nfirms] = zeros(nfirms-i,1);   endif;
  endif;         @ else integrating over random exit in calcval @

  @ entry probability based on call to chkentry() just before calling optimize() @
  @ if phiH==0 then isentry based on pre-exit, otherwise based on integrating over random exit in calcval @

  entered = x_exit_entry[2*nfirms+1];    @ non-foc: isentry[qencode(rev(sortc(locwx,1)))]; @

  locwe  = locwx;               @ locwe doesn't ensure room for entrant since entered = 0 if not --> locwe ignored @
  locwK = locwx;                @ locw when entrant starts at kmax @
  locwe[nfirms] = entry_k;      @ entrant will get moved up to its proper position by sortc in calcval @
  locwk[nfirms] = kmax;
  eexit = oexit;
  eexit[nfirms] = 0;            @ ensures entrant does not exit in integration in calcval() @

  @ when nfirms > 2, only do leap if TWO open spots so can still check frequency of full industry for need to increase nfirms @
  doleap = RLG_leap>0 AND nfirms>1;   if doleap AND nfirms>2;   doleap = locw[nfirms-1]==0 ;   endif;
  
  @ Now calculate the optimal policies for this industry structure, @
  @ given that entry and exit are as specified.   @
  for j (1,nfirms,1);
    if RLG_w0_exit AND locw[j] == 0;              @  w=0  assumed to exit since w=0 is signal of space for entrant @
      v1 = phi;  
      if phiH /= 0;  v1 = (phiH+phi)/2;   endif;     @  if locw[j] is 0 then so are all higher j @
      nval[j:nfirms] = v1*ones(nfirms-j+1,1);     @  can shut this off if entry costs, scrap chosen s.t. no entry/exit @
      nx[j:nfirms]   = zeros(nfirms-j+1,1);
      nexit[j:nfirms]= ones(nfirms-j+1,1);
      break;
    endif;

    v1=0; v2=0;
    if entered < 1;
      @ First: Calculate v, without entry @
      {v1, v2} = calcval(j,locwx,ox,oexit,locw[j]);
    
      /* if ddebug == 9;  "XX start";  calcvalXX(j,locwx,ox,oexit,locw[j]); " XX done";  print /flush;;  ddebug = 0;  endif; */
    
    endif;

    if entered > 0;
      @ A firm wants to enter with positive probability @
      @ if nfirms>1;  dop=1;  endif; @
      {tempv1, tempv2} = calcval(j,locwe,ox,eexit,locw[j]);
    
      /* if ddebug == 9;  "XX start";  calcvalXX(j,locwx,ox,oexit,locw[j]); " XX done";  print /flush;;  ddebug = 0;  endif;   */
    
      if doleap ;
	{tempk1, tempk2} = calcval(j,locwk,ox,eexit,locw[j]);
	tempv1 = (1-RLG_leap)*tempv1 + RLG_leap*tempk1;
	tempv2 = (1-RLG_leap)*tempv2 + RLG_leap*tempk2;
      endif;
      v1 = entered*tempv1 + (1-entered)*v1;
      v2 = entered*tempv2 + (1-entered)*v2;
      @ if nfirms>1; format /m1 /rds 3,0; "j " j " locwe " locwe';;  format /m1 /rd 12,8; "entered " entered " tempv1~v2 " tempv1 tempv2 " v1~v2 " v1 v2 ;  endif; @
    endif;
    
    @ Calculate values for firm, given that it is not leaving @

    if v1 <= v2;   r = 1.0;    @ Avoid division by zeros @
    else;          r = 1.0/(beta*aeff[locw[j]+1]*(v1-v2));
    endif;
    
    @ r now contains the value r = (1 - p)^2. => p = 1 - sqrt(r)),  @
    @ where p is the optimal prob. of having k rise, cond. on world @
    r = minc(maxc(r|0.0000000000001)|1);
    p = 1.0 - sqrt(r);
    nx[j] = p/(aeff[locw[j]+1] - aeff[locw[j]+1] * p);
  
    if BatchJob == 0 AND locw[j]==kmax AND locw[minc(nfirms|2)]>0;  nx[j] = maxc(RLG_minx1|nx[j]);  endif;   @ force min investment by frontier firms, unless monopolist @
  
    @ if nx[j]< .018;  nx[j] = .018;  endif; @         @ RLG imposed min x to avoid non convergence @

    @ Now calculate the value from staying in @
    @ Ask: given this optimal investment level, will there be exit? @

    @ PM timing:  nval[j] = profit[w,j] - nx[j] + beta*(v1*p + v2*(1-p));  @
    @             then check nval[j] <= phi to determine exit              @
    @ RG timing:  profit[w,j] earned whether EXIT or not                   @
    @             so check  beta*(v1*p + v2*(1-p)) - nx[j] <= phi for exit @
    
    nval[j] = beta*(v1*p + v2*(1-p)) - nx[j];       @ add  profit[w,j]  after exit decision @

    @ "w~j~nval[j]~profit[w,j]~nx[j]~p~v1~v2  ";;     w~j~nval[j]~profit[w,j]~nx[j]~p~v1~v2; @

    if phiH==0;                                     @ fixed scrap @
      if nval[j] <= phi;
	nexit[j:nfirms]= ones(nfirms-j+1,1);        @ all weakly lower firms exit too      @
	nval[j:nfirms] = ones(nfirms-j+1,1) * phi;  @ ox already 0 for lower firms         @
	nx[j] = 0  ;                                @ fixed scrap --> next firms also exit @
	break;
      endif;
    else;                                           @ random scrap @
      if nval[j] < phiH;                            @ Pr(exit) > 0 @
	if nval[j] <= phi;                          @ Pr(exit) = 1 @
	  nexit[j] = 1;
	  nval[j]  = (phi+phiH)/2;                  @ random scrap --> next j may not exit @
	  nx[j]    = 0;
	else;                                       @ assumes scrap ~ U(phi,phiH) @
	  nexit[j] = (phiH-nval[j])/(phiH-phi);     @ nval[j] uses  E(scrap|exit) @
	  nval[j]  = (1-nexit[j])*nval[j] + nexit[j]*(phiH+nval[j])/2;
	endif; 
      endif;
    endif;
    
    nval[j] = nval[j] + profit[w,j];   @ RLG: profits earned even when exit at END of period @
    
    if 0;          @ using this leads to identical firms taking different actions: 1st of identical pair "moves first" @
      oexit[j] = oexit[j]*RLG_damp + (1-RLG_damp)*nexit[j];    @ nval[j]>= phiH --> nexit[j] remains 0, nval[j] as-is  @
      ox[j]    =    ox[j]*RLG_damp + (1-RLG_damp)*nx[j];       @ ox, oexit updates optional? since oldx = newx later ? @
    endif;
    
    if phiH==0;
      locwx[j] = (nval[j] > phi)*locw[j];      @ implements new exit policy when fixed scrap @
      locwe[j] = locwx[j];
    endif;

  if 0; format /m1 /rds 3,0; "j " j " w " locw';;  format /m1 /rd 12,8; " v1~v2 " v1 v2 " nval " nval' " nx " nx' " nexit " nexit';  endif;
    
    if j>1 AND 0 ;
      if locw[j]==locw[j-1] AND (abs(nval[j]-nval[j-1])>1e-8 OR abs(nx[j]-nx[j-1])>1e-8 OR abs(nexit[j]-nexit[j-1])>1e-8);
	format /m1 /rds 3,0; "j " j " w " locw';;  format /m1 /rds 7,3; " v " nval' " x " nx' " newexit " nexit';; format /m1 /re 10,1; nval[j]-nval[j-1] nx[j]-nx[j-1] nexit[j]-nexit[j-1];
      endif;
    endif;
    
  endfor;
  retp(nx',nval',nexit');
endp;


proc (0) = chkentry_w(x_exit_entry,w);
@ This procedure calculates for which value of other people's omegas, would @
@ a firm want to enter, given that the market has room for another firm.    @
@ Implicit parameters are oldx, oldvalue (passed in) and isentry (returned) @
  local @ w,@ locw,v1,vgarbage, doleap, val;                @ val = Value from entering @

  @ With FIXED scrap (i.e., orig PM), exit is implemented before looking up    @
  @ the isentry[] element, so do NOT implement exit here when getting isentry  @
  @ With RANDOM scrap, exit is integrated over in calcval()                    @
  
  @ for w (1,wmax,1); @     @ non-foc approach called chkentry() for ALL states prior to loop over w in which call optimize(w) @
    locw = dtable[.,w];     @ qdecode(w); @    @ calcval integrates over exit, investment outcomes later this period @
    if locw[nfirms] == 0;   @ room for entrant @
      {vgarbage,v1} = calcval(nfirms,locw,x_exit_entry[1:nfirms],x_exit_entry[nfirms+1:2*nfirms],entry_k);   @ non-foc:  calcval(nfirms,locw,oldx[w,.]',oldexit[w,.]',entry_k); @
      
    /* if ddebug == 9;  "XX start in chkentry_w";  calcvalXX(nfirms,locw,x_exit_entry[1:nfirms],x_exit_entry[nfirms+1:2*nfirms],entry_k);  " XX done in chkentry_w";   ddebug = 0; print /flush;;  endif; */
    
      val = beta * v1;

      @ when nfirms > 2, only do leap if TWO open spots so can still check frequency of full industry for need to increase nfirms @
      doleap = RLG_leap>0 AND nfirms>1;   if doleap AND nfirms>2;   doleap = locw[nfirms-1]==0 ;   endif;
    
      if doleap;        @ enter at frontier with prob RLG_leap @
	{vgarbage,v1} = calcval(nfirms,locw,x_exit_entry[1:nfirms],x_exit_entry[nfirms+1:2*nfirms],kmax);
	val = (1-RLG_leap)*val + RLG_leap * beta * v1;
      endif;
      @ format /m1 /rd 7,3; "w~val " locw' val; @   @ This val does not subtract entry costs -->  It's the post entry cont. value @
      v_entry[w] =  val;
      isentry[w] = maxc(0|minc(1|((val - x_entryl) / (x_entryh - x_entryl))));
    endif;
  @ endfor;
  isentry = minc((isentry~ones(wmax,1))');
  isentry = maxc((isentry~zeros(wmax,1))');
  @
endp;


proc newton3(p,&objfunk,w,doo,tol);
@ This procedure performs a simple Newton-Raphson search to find the root of the function objfunk.
  LB[doo] <= p <= UB[doo] is enforced
  The policy parameters in p are ordered  x ~ exit ~ entry.
  doo  has the indices from x_exit_entry being tweaked to find a root.
  tol is param in case want to start with loose tol and tighten it as near converged value function
@
  local objfunk:proc,deriv,pnew,epsilon,x,i,iter,maxiter,np,dp,minstep,trymax,LBw;

  LBw = LB;   LBw[1:nfirms] = RLG_minx1*(dtable[.,w].==kmax);     @ frontier firms' min investment @
  
  np = rows(p);    minstep = 1;                   @ stepsize is random uniform between [minstep,1] @
  trymax= 200*np;  deriv   = zeros(np,np);
  iter = 0;        maxiter = trymax*20;
  
  do while iter < maxiter;
    iter=iter+1;

    x = objfunk(p,w,doo);      @ Calculate function at p @

    if maxc(abs(x))<tol;  break;  endif;

    for i (1,np,1);            @ Calculate derivative matrix @
      dp = p;
      if UB[doo[i]]-p[i] < 0.0001;  dp[i] = p[i]-.0001;   @ decrease p to get deriv instead since too near upper bound @
      else;                         dp[i] = p[i]+.0001;
      endif;
      deriv[.,i] = (objfunk( dp, w,doo) - x)/(dp[i]-p[i]);
    endfor;

    @ if iter>2500;  format /m1 /rds 6,2;  print "w " w "  doo " doo' "  iter " iter "  p " p' "  foc " x';  endif; @

    trap 1;    pnew = inv(deriv);  trap 0;

    if scalerr(pnew) OR (iter%trymax==0);
      format /m1 /rd 2,0;
      if scalerr(pnew);                         @ reset iter to avoid quick restart due to iter%trymax check @
	"restarting newton() at rndu() since singular deriv at iter " iter " with w= " dtable[.,w]' " doo = " doo';; format /m1 /re 8,1; " p= " p' " foc= " x' " deriv= "; deriv;   @ vecr(deriv)'; @
	iter = trymax*(1+floor(iter/trymax));
      else;
	"restarting newton() at rndu() since not converging at iter " iter " with w= " dtable[.,w]' " doo = " doo';;  format /m1 /re 8,1; " p= " p' " foc= " x';;
      endif;
      pnew = LBw[doo]+rndu(np,1).*(minc(UB[doo]'|10*ones(1,np))-LBw[doo]);
      " pnew= " pnew';                          @ rndu between LB and UB, with cap at 10 above LB      @
      minstep = maxc(.3|minstep-.1);            @ singular deriv less likely if use smaller stepsize?  @
    else;
      if iter%trymax<5;   pnew = p - pnew * ((.1+(iter%trymax)/10)*rndu(1,1)) *x;     @ start slow @
      else;               pnew = p - pnew *   (minstep+(1-minstep)*rndu(1,1)) *x; 
      endif;
    endif;
    pnew = minc(UB[doo]'|pnew');
    pnew = maxc(LBw[doo]'|pnew');
    p = pnew;
  endo;
  if maxc(abs(x))>tol;  format /m1 /rds 2,0; "Newton failed with w = " dtable[.,w]' " doo = " doo';; format /m1 /re 8,1; "  returning p= " p' "  foc= " x' ;   endif;
  retp(p);
endp;


proc (2) = calcvalXX(place,w,x,isex,k);
@ This procedure calculates val = EEEV(.,.,.,.)p(.)p(.)p(.), where E @
@ represents sums, and this is the calculation of the 4-firm problem @
@ Vars: place = place of own omega, for calculating value function (v) @
@       w = the vector of omegas; already decoded @
@       x = the vector of investments (nfirms of them) @
@    isex = the vector of exit probabilities (nfirms of them) @
@ Implicit parameter: oldvalue @
@ For efficiency reasons, it outputs the following vector:  @
@ { calcval(k_v+1,w,x), calcval(k_v,w,x) }  @
  local i,valA,valB,d,e,probmask,z1,z2,locmask,unboundedGG,
        p_up,  @ p_down, p of going up/down for all other firms @
        temp,
        pl1,justone, w1, probmask_chksum,
	iEXIT,EXprobmask;  @ for integrating over exit if phiH>0 @

  z1 = zeros(nfirms,1);
  z2 = kmax*ones(nfirms,1);

  /* ddebug = 9; */
  
  unboundedGG = 0;      @ last AND checks whether 1st firm alone at frontier @
  if RLG_wstar==0 AND w[place]==kmax AND ( nfirms==1 );   @ USED to have AND (nfirms==1 OR (w[1]-w[minc(2|nfirms)])>0) @
    unboundedGG = 1;                                      @ but now realize unbounded GG (RLG_wstar=0) does NOT work in oligopoly @
  endif;                                                  @ so use bounded GG (RLG_wstar= -1) to get value and policy functions @
  @ unboundedGG = 1 --> adjust valA since firm at kmax @  @ but  unbounded GG in welf_ma.g to get higher dynamic CS @

  @ Expand mask to allow for non-inclusion of the  place_th  firm (since not integrating over its outcomes)  @
  if nfirms > 1;
    if place == 1;           locmask = zeros(1,two_n)|mask;
    elseif place == nfirms;  locmask = mask|zeros(1,two_n);
    else;                    locmask = mask[1:place-1,.]|zeros(1,two_n)|mask[place:nfirms-1,.];
    endif;
  else;                      locmask = zeros(1,1);
  endif;
  x[place] = 0;
  w[place] = k;
  isex[place] = 0;             @ place firm stays put --> no investment or exit @
  justone = zeros(nfirms,1);
  justone[place] = 1;

  w1 = w;  @ w1 = rev(sortc(w1,1)); @     @ HEY:  w1 is only used in printf statements  @

  i = aeff[w+1].*x;    p_up = i./(1+i);    @ p_down = 1 - p_up; @

  valA = 0;  valB = 0;           @ valA is value if w up, valB if w same @
  probmask_chksum = 0;
  
  @ outer loop over    exit    outcomes @
  @ inner loop over investment outcomes @
  @ two_n = 2^(nfirms-1)  NOT  2^nfirms @

  if maxc(isex)>1; format /m1 /rds 3,0; " w " w';;  format /m1 /rd 7,3; " isex " isex';  endif;

  @ if ddebug==9;  format /m1 /rd 2,0; "locmask (nfirms by 2^(nfirms-1))" locmask; "G: place= " place " w= " w';; format /m1 /rds 16,14;  print " isex= " isex' " p_up= " p_up';   endif; @

  for iEXIT (1,two_n,1);       @ break issued at endfor if phiH==0 --> fixed scrap --> exit handled in optimize()  @
    EXprobmask = prodc(2 .* locmask[.,iEXIT] .* isex + 1 - locmask[.,iEXIT] - isex);
    
  @ if ddebug==9;  format /m1 /rd 2,0;  "G: iEXIT= " iEXIT " EXITmask= " locmask[.,iEXIT]'  " EXprobmask = ";; format /m1 /rds 16,14; print EXprobmask "  p_up= " p_up';   endif; @
  
    if 0;
      format /m1 /rd 6,0;
      " ";
      "two_n = " two_n;
      "size locmask: " rows(locmask) rows(locmask');
      "size    isex: " rows(isex) rows(isex');
      "locmask";  print locmask;
    endif;
  
    @ if nfirms>1;  if w[1]==12 AND w[2]==1 AND place==2;  if iEXIT==1; " ";  endif; format /m1 /rds 3,0; place~9~w'~9~locmask[.,iEXIT]';; format /m1 /rds 8,4; isex'~EXprobmask;  endif;  endif; @

    @ if dop;  format /m1 /rd  3,0;  "iEXIT " iEXIT " w " w' " place~k " place k " locmask " locmask[.,iEXIT]' " EXprobmask>0" (EXprobmask>0);  endif;  @
    
    if EXprobmask>0 OR phiH==0;  @ else skip inner loop since zero probability @
      
      for i (1,two_n,1);
        @ probmask = prodc(mask[.,i] .* p_up + (1 - mask[.,i]) .* p_down); @
	probmask = prodc(2 .* locmask[.,i] .* p_up + 1 - locmask[.,i] - p_up);
	
	@ if ddebug==9;  format /m1 /rd 2,0;  "G: iEXIT= " iEXIT " i= " i " EXITmask= " locmask[.,iEXIT]'  " locmask= " locmask[.,i]' " EX~probmask = ";; format /m1 /rds 12,10; print EXprobmask~probmask;   endif; @
	
	if probmask>0;  @ nonindented if-then @
	  
	  
	d = w+locmask[.,i];

	if phiH /= 0;                        @ phiH >0 --> random scrap --> integrate over exit here in calcval @
	  d = d.*(1-locmask[.,iEXIT]);    @ phiH==0 --> exit handled in optimize before calling calcval      @
	  probmask = probmask*EXprobmask; @ d = 0 for exiting firms @
	endif;
	
	probmask_chksum = probmask_chksum + probmask;    @ should be 1 when done @
	
	temp = rev(sortc(d~justone,1));   @ sorts via column 1 (i.e., d) so can look up in qencode @

	if 0;
	  format /m1 /rd 3,0;
	  "d ~ justone ~ w ~ locmask[.,i] ~ temp[2 cols] ~ isex   place= " place; d~justone~w~locmask[.,i]~temp~isex;
	  temp[99];
	endif;
	
	d = temp[.,1];                    @ 2nd col of temp has 1 in row of firm j, obtained by pl1 = maxindc() @
	e = d - 1;
	@ Check for evaluation of value fn. at -1 @
	e = maxc((e~z1)');

	if RLG_out==0;
	  if e[1]<kmax;   e = e + (e.>0)*(kmax-e[1]);	endif;   @ all frontier firms must hav exited --> move remaining firms up s.t. highest at kmax @
	  if d[1]<kmax;   d = d + (d.>0)*(kmax-d[1]);	endif;   @ all frontier firms must hav exited --> move remaining firms up s.t. highest at kmax @
	endif;
	
	if RLG_no_force_exit;             @ Never bump firms off lowest rung of ladder @
	  @ "1 ";  "e= " e';  "d= " d'; "e+d" e'+(d'.==1); @
	  e = e+(d.==1);                  @ This can move the focal firm into a tie with more firms, @ 
	endif;                            @ but his pl1 based on temp[.,2] still valid since tied firms have same values @

	pl1 = maxindc(temp[.,2]);   @ sumc(d[1:place].>=k) + sumc(d[place:nfirms].>k);@ 
	
        @ if ddebug==9;  format /m1 /rd 8,4;   "G: e = " e'   " B pl1 qencode(e) v "  pl1 qencode(e) oldvalue[qencode(e),pl1];;  format /m1 /rds 16,14;  "  chksum= " probmask_chksum;  endif; @
      
	if RLG_wstar <=0 AND sumc(d.==kmax+1) > 0;  @ RLG: at least one firm beyond "frontier" so use delta = 1.0 --> only use  e  @
	  valB = valB + (                                            oldvalue[qencode(e),pl1] )*probmask;
	else;
	  d = minc((d~z2)');  @ really only needed by PM statespace (ie, if pp_[_WSTAR]>0) but harmless otherwise since previous if already checked d==kmax+1) @
	  if delta>0;
	    valB = valB + ( (1-delta)*oldvalue[qencode(d),pl1] + delta*oldvalue[qencode(e),pl1] )*probmask;
	  else;
	    @ if ddebug==9;  format /m1 /rd 8,4;   "G: d = " d'   " B pl1 qencode(d) v "  pl1 qencode(d) oldvalue[qencode(d),pl1];  endif; @
	    valB = valB +             oldvalue[qencode(d),pl1] * probmask;
	  endif;
	endif;
	
	@ format /m1 /rds 3,0;  "B: place= " place   " w1= " w1'  " d= " d'  " e= " e'  " just1= " justone' "  pl1= " pl1  "  k= " k;  @
	
	
	d = w+locmask[.,i]+justone;       @ +justone is the successful innovation @

	if phiH /= 0;                     @ phiH >0 --> random scrap --> integrate over exit here in calcval @
	  d = d.*(1-locmask[.,iEXIT]);    @ phiH==0 --> exit handled in optimize before calling calcval      @
	endif;
	
	temp = rev(sortc(d~justone,1));
	d = temp[.,1];
	e = d - 1;
	@ Check for evaluation of value fn. at -1 @
	e = maxc((e~z1)');
	
	if RLG_out==0;
	  if e[1]<kmax;   e = e + (e.>0)*(kmax-e[1]);	endif;   @ all frontier firms must have exited --> move remaining firms up s.t. highest at kmax @
	  if d[1]<kmax;   d = d + (d.>0)*(kmax-d[1]);	endif;   @ all frontier firms must have exited --> move remaining firms up s.t. highest at kmax @
	endif;
	
	if RLG_no_force_exit;             @ Never bump firms off lowest rung of ladder @
	  @ "2 ";  "e= " e';  "d= " d'; "e+d" e'+(d'.==1); @
	  e = e+(d.==1);                  @ This can move the focal firm into a tie with more firms, @ 
	endif;                            @ but his pl1 based on temp[.,2] still valid since tied firms have same values @

	pl1 = maxindc(temp[.,2]);   @ sumc(e[1:place].>=k) + sumc(e[place:nfirms].>k);@

	if RLG_wstar <= 0 AND sumc(d.==kmax+1) > 0;  @ RLG: at least one firm beyond "frontier" so use delta = 1.0 --> only use  e  @

	  /* see   eql_ma.approximation_notes */

	  valA = valA + (   /* REMOVED: (1-delta)*() + delta*  */     oldvalue[qencode(e),pl1] )*probmask;

	  @ if ddebug==9;  format /m1 /rd 8,4;   "G: e = " e'   " A pl1 qencode(e) v "  pl1 qencode(e) oldvalue[qencode(e),pl1]  "  here";  endif;  @

	
	  if unboundedGG;    @ only adjust valA (not valB) since approximation needed only when the "place" firm advances beyond kmax.   nfirms==1 required for unboundedGG=1 @

	    @ Unbounded GG approximation occurs here.  Bounded GG (RLG_wstar== -1) and Unbounded GG (RLG_wstar==0) with the nfirms==1 restriction @
	    
	    @ valA = valA + ( ( profit[encode2(e+justone),1] - profit[qencode(e),1] ) *(1-delta) /(1-beta) * 1  )*probmask; @    @ <-- NAILS approximation when 1 firm for any delta @ 
	    valA   = valA + ( ( profit[qencode(e+justone),1] - profit[qencode(e),1] ) *(1-delta) /(1-beta) * 1  )*probmask;     @ <-- NAILS approximation when 1 firm for any delta @

	    if 0 AND nfirms > 1 AND e[nfirms]==24;
	      format /m1 /rds 3,0;  "place= " place   " w1= " w1'  " d= " d'  " e= " e'  " just1= " justone' " e-1= " maxc((e-1~z1)')'  " pl1= " pl1 ;;
	      format /m1 /rds 6,3;  " probmask= " probmask   " v(e)= "  oldvalue[qencode(e),pl1]   "  v(e-1)= " oldvalue[qencode(maxc((e-1~z1)')),pl1];;
	      " v diff=" oldvalue[qencode(e),pl1]-oldvalue[qencode(maxc((e-1~z1)')),pl1];;
	      " profdiff~ /1-beta= "  ( profit[encode2(e+justone),1] - profit[qencode(e),1] )  ( profit[encode2(e+justone),1] - profit[qencode(e),1] )/(1-beta)*(1-delta);
	    endif;
	    
	    @ valA = valA + ( oldvalue[qencode(e),pl1] - oldvalue[qencode(maxc((e-1~z1)')),pl1]  )*probmask; @          @ <-- Nails approximation when 1 firm AND delta = 0 @
	    
	    @ valA = valA + ( profit[encode2(e+justone),1] - profit[qencode(e),1]  + oldvalue[qencode(e),pl1]-oldvalue[qencode(maxc((e-1~z1)')),pl1]  )*probmask; @
	    
	    @ valA = valA + ( profit[encode2(e+justone),1] - profit[qencode(e),1] ) * probmask; @
	    
	    @ valA = valA + ( ( pp[_MKT_SIZE]*pp[_RLG_WSCALE]*(1-1/(1+exp(kmax-1))) ) *(1-delta) /(1-beta) )*probmask; @
	    
	  endif;
	else;
	  d = minc((d~z2)');  @ really only needed by PM statespace (ie, if pp_[_WSTAR]>0) but harmless otherwise since previous if already checked d==kmax+1) @
	  if delta>0;
	    valA = valA + ( (1-delta)*oldvalue[qencode(d),pl1] + delta*oldvalue[qencode(e),pl1] )*probmask;
	  else;
	    valA = valA +             oldvalue[qencode(d),pl1] * probmask;
	  
	  @ if ddebug==9;  format /m1 /rd 8,4;  "G: d = " d'  " A pl1 qencode(d) v "  pl1 qencode(d) oldvalue[qencode(d),pl1];  endif; @

	  endif;
	endif;
	
	endif;  @ if probmask>0 @
	/*
	if dop;
	  format /m1 /rd  3,0;  "iExit~i " iExit i " w " w' " place~k " place k " locmask " locmask[.,i]' " d " d' "  e " e' " pl1 " pl1 " qen(d~e) "  qencode(d) qencode(e);;  
	  format /m1 /rd 12,8;  " valA~B " valA valB  " probmask " probmask " oldvalue[qencode(d~e),pl1] " oldvalue[qencode(d),pl1] oldvalue[qencode(e),pl1] ;
	endif;
	*/
	
	endfor;  @ loop over investment outcomes @
	
      endif;     @ check whether EXprobmask==0 in which case inner loop skipped @
    
      if phiH==0;  break;  endif;   @ exit handled by optimize() before calling calcval --> ignore outer loop @

    endfor;    @ loop over    exit    outcomes @

    if abs(1-probmask_chksum)>1e-10;  
      " ";  "Gauss HEY: probmask_chksum not 1: ~probmask~EXprobmask  " probmask_chksum probmask EXprobmask "  locmask next block:";
      format /m1 /rds 3,0; locmask;  " ";    "place ~ k ~ 9 ~ w " place~k~9~w';; format /m1 /rds 8,4; " isex " isex'  " p_up " p_up';
    endif;
  
  @ if ddebug==9;  format /m1 /re 6,1;  "Gauss code  1-probmask_chksum= " 1-probmask_chksum;;   format /m1 /rds 8,4;  "valA~B " valA valB;  endif; @
  
  retp(valA,valB);
endp;


@  proc (2) = calcvalXX(place,w,x,isex,k);  retp(1,1);  endp; @
 

proc (2) = calcval(place,w,x,isex,k);                                              @  C version of calcval()   Hardcode the desired version by renaming the undesired version  @
@ This procedure calculates val = EEEV(.,.,.,.)p(.)p(.)p(.), where E @
@ represents sums, and this is the calculation of the 4-firm problem @
@ Vars: place = place of own omega, for calculating value function (v) @
@       w = the vector of omegas; already decoded @
@       x = the vector of investments (nfirms of them) @
@    isex = the vector of exit probabilities (nfirms of them) @
@ Implicit parameter: oldvalue @
@ For efficiency reasons, it outputs the following vector:  @
@ { calcval(k_v+1,w,x), calcval(k_v,w,x) }  @
  local i,valA,valB,d,e,probmask,locmask,unboundedGG,
        p_up,  @ p_down, p of going up/down for all other firms @
        temp,
        pl1,probmask_chksum ;
	
  @ Expand mask to allow for non-inclusion of the  place_th  firm (since not integrating over its outcomes)  @
  if nfirms > 1;
    if place == 1;           locmask = zeros(1,two_n)|mask;
    elseif place == nfirms;  locmask = mask|zeros(1,two_n);
    else;                    locmask = mask[1:place-1,.]|zeros(1,two_n)|mask[place:nfirms-1,.];
    endif;
  else;                      locmask = zeros(1,1);
  endif;
  x[place] = 0;
  w[place] = k;
  isex[place] = 0;             @ place firm stays put --> no investment or exit @

  p_up = aeff[w+1].*x;    
  p_up = p_up./(1+p_up);      @ p_down = 1 - p_up; @

  
  /* ddebug = -1;  do while ddebug /= 0 AND ddebug /= 9; */
  
  valA = 0;  valB = 0;         @ valA is value if w up, valB if w same @
  probmask_chksum = 0;

  dllcall calcval( valA, valB, probmask_chksum, place, w, x, isex, locmask, kmax, nfirms, p_up, two_n, phiH, RLG_no_force_exit, RLG_wstar, oldvalue, etable1, multfac1, binomv, delta, RLG_out);
  
  if abs(1-probmask_chksum)>1e-10;  
    format /m1 /re 6,1; " ";  "Ccode HEY: probmask_chksum not 1: 1-chksum= " 1-probmask_chksum;;
    format /m1 /rds 3,0; " ";    "place ~ k ~ 9 ~ w " place~k~9~w';; format /m1 /rds 8,4; " isex " isex'  " p_up " p_up'   "valA~B " valA valB;
    /*
    if ddebug==1;   ddebug = 9;  "setting ddebug=9";     @ will now call gauss version at same state for comparison @
    else;           ddebug = 1;  "setting ddebug=1";     @ will repeat C call with debug turned on @
    endif;
  else;
    ddebug = 0;
   */

  endif;
  
  if 0;  format /m1 /re 6,1;  " C code:  1-probmask_chksum= " 1-probmask_chksum;;  format /m1 /rds 12,8;  "valA~B " valA valB;   endif;

/* endo; */
  
  retp(valA,valB);
endp;


proc encode2(ntuple);    @ already have a  proc encode() in welf_ma.g @
@ This procedure takes a weakly descending n-tuple (n = nfirms)       @
@ with min. elt. 0, max. elt. kmax, and encodes it into an integer    @

retp(1+sumc(binomv[(ntuple+nfirms_oneton)*binomcols+ntuple+1]));   @ RLG: binomv is vecr(binom) and index based on gauss being row-major @
/*
  local code,digit,i;
  code = 1;               @ Coding is from 1 to wmax @
  for i (1,nfirms,1);
    digit = ntuple[i];
    code = code + binom[digit+nfirms+1-i,digit+1];
  endfor;
  retp(code);
*/
endp;


proc qencode(ntuple);
@ This procedure does a quick encode of any n-tuple given in weakly descending order. @
@ Encoding uses a table lookup. Each column of the table consists of an n-tuple;      @
@ the ith column is the ith n-tuple to be decoded.  The table is stored in etable.    @

@ RLG: binomv is vecr(binom) and index based on gauss being row-major                 @
@ <10% slower than old qencode, but does not need etable1 (massive for high nfirms)   @


  if nfirms < maxfirms_qencode;  
    retp(etable1[sumc(ntuple.*multfac1)+1]);
  else;
    retp(1+sumc(binomv[(ntuple+nfirms_oneton)*binomcols+ntuple+1]));    @ about 10% slower than etable1 @
  endif;
  /*
  if nfirms <= encfirm;
    retp(etable1[sumc(ntuple.*multfac1)+1]);
  else;    
    retp(etable1[sumc(ntuple.*multfac1)+1] + etable2[sumc(ntuple.*multfac2)+1]);
  endif;
  */
endp;

/*
proc qdecode(code);
@ This procedure does a quick decode of a previously encoded number into   @
@ a weakly descending n-tuple. Decoding is done using a table lookup. Each @
@ column of the table consists of an n-tuple; the ith column is the ith    @
@ n-tuple to be decoded. The table is stored in the variable "dtable".     @
  retp(dtable[.,code]);    @ Why a function call to only access a matrix?? @
endp;
*/

proc decode(code);
@ This procedure takes a previously encoded number, and decodes it into @
@ a weakly descending n-tuple (n = nfirms)                              @
@ ONLY called to construct dtable (table of w's)                        @

  local ntuple,digit,i;
  code = code-1;
  ntuple = zeros(nfirms,1);
  for i (1,nfirms,1);
    digit = 0;
    do while binom[digit+nfirms-i+2,digit+2] <= code;
      digit=digit+1;
    endo;
    ntuple[i] = digit;
    code = code-binom[digit+nfirms-i+1,digit+1];
  endfor;
  retp(ntuple);
endp;

