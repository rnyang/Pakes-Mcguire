/* ??  compile with:    gcc -g -O2 -o calcval calcval.c -lm -lc -shared -fpic     (or same options, but with  icc  instead of  gcc  */
#define min(x,y)   ((x) < (y) ? (x) : (y))
#define max(x,y)   ((x) > (y) ? (x) : (y))
#define max_nfirms         16  /* enables declarations to create vectors instead of dynamically allocating based on actual nfirm.  possibly faster if a power of 2 */
#define maxfirms_qencode   9   /* maximum nfirms for which etable1 constructed for use in qencode().  This hardcode must match the hardcoded value in eql_ma_foc.g */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* Straight Insertion sort from Numerical recipes in C, 2nd ed p. 331.  Fast for N<20, which will always be the case for us */
/* RLG changed datatypes of arr and brr (and a,b) from float to int, since we are only sorting integers here */
/* RLG also changed the comparison  arr[i] > a  to   arr[i] < a  to sort in reverse order */
 void piksr2(int n, int arr[max_nfirms], int brr[max_nfirms])
 {
   int i,j,a,b;

   /* printf("\n n= %d \n",n); for (j=0;j<n;j++) printf(" %d", arr[j]); */
   
   for (j=1;j<n;j++) {                /* numerical recipes uses j=2; j<=n; j++  and ignores arr[0] element.  WHY?? */
     a=arr[j];
     b=brr[j];
     i=j-1;                           /* numerical recipes uses i > 0 in while() but RLG uses i=0 when changed j   */
     while (i >= 0 && arr[i] < a) {   /* tried descending sort via < a, since RLG flipped arr[i] > a to arr[i] < a */
       arr[i+1]=arr[i];
       brr[i+1]=brr[i];
       i--;
     }
     arr[i+1]=a;
     brr[i+1]=b;
   }

   /* printf("\n");  for (j=0;j<n;j++) printf(" %d", arr[j]);  printf("\n"); */
 }


/* qencode() does a quick encode of any n-tuple given in weakly descending order.
 Encoding uses a table lookup. Each column of the table consists of an n-tuple;      
 the ith column is the ith n-tuple to be decoded.  The table is stored in etable.    

 RLG: binomv is vecr(binom) and index based on gauss being row-major                 
 <10% slower than old qencode, but does not need etable1 (massive for high nfirms)   
*/

static unsigned long qencode( int ntuple[max_nfirms], int nfirms, double *etable1, double *multfac1, int binomcols, double *binomv)

{
  int j;
  unsigned long i;
  
  i = 0;  /* returns value 1 less than the row itself, since C index starts at 0 */
  
  if(nfirms < maxfirms_qencode) {
    for(j=0;j<nfirms;j++)
      i += ntuple[j] * (int)*(multfac1+j);
    i = (unsigned long)*(etable1+i)-1;
  }
  else                           /* in Gauss about 10% slower than etable1 lookup */
    for(j=0;j<nfirms;j++)
      i += (unsigned long)*(binomv + (ntuple[j]+nfirms-j-1)*binomcols + ntuple[j]);
  
  /* printf(" j= %d  arg= %d  ", j, (ntuple[j]+nfirms-j-1)*binomcols + ntuple[j]); */
  
  /* printf("  w = %d %d %d %d  i= %u \n", ntuple[0], ntuple[1], ntuple[2], ntuple[3], i); */
  return(i);
}


/*  this comment has the gauss function in eql_ma_foc.g which calls this c function

proc (2) = calcval_C(place,w,x,isex,k);
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
        pl1,probmask_chksum,
	iEXIT,EXprobmask;  @ for integrating over exit if phiH>0 @

  p_up = aeff[w+1].*x;    
  p_up = p_up./(1+p_up);      @ p_down = 1 - p_up; @
	  
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

  valA = 0;  valB = 0;           @ valA is value if w up, valB if w same @
  probmask_chksum = 0;


... call c    

    if abs(1-probmask_chksum)>1e-10;  
      " ";  "HEY: probmask_chksum not 1: ~probmask~EXprobmask  " probmask_chksum probmask EXprobmask "  locmask next block:";
      format /m1 /rds 3,0; locmask;  " ";    "place ~ 9 ~ w " place~9~w';; format /m1 /rds 8,4; " isex " isex'  " p_up " p_up';
    endif;
    
 */


/* valA, valB, probmask_chksum are the implicit returned values */

int calcval( double *valA, double *valB, double *probmask_chksum, double *dplace, double *dw, double *x, double *isex, double *dlocmask, double *dkmax, double *dnfirms, double *p_up, double *dtwo_n, double *phiH, double *dRLG_no_force_exit, double *dRLG_wstar, double *oldvalue, double *etable1, double *multfac1, double *binomv, double *delta, double *dRLG_out /*, double *debug */)
{
  double EXprobmask, probmask;
  int i,j,pl1,d[max_nfirms],e[max_nfirms], w[max_nfirms], z2[max_nfirms],justone[max_nfirms],temp[max_nfirms],kmax,nfirms,place,iEXIT,two_n,locmask[max_nfirms],EXITmask[max_nfirms],RLG_no_force_exit,RLG_wstar,RLG_out,binomcols,qencode_e,qencode_d;

  /* convert doubles from Gauss to int */
  RLG_no_force_exit = (int)*dRLG_no_force_exit;
  RLG_wstar = (int)*dRLG_wstar;
  RLG_out = (int)*dRLG_out;
  place  = (int)*dplace-1;           /* place-1 since C index starts at 0, Gauss starts at 1 */
  nfirms = (int)*dnfirms;
  two_n = (int)*dtwo_n;
  kmax = (int)*dkmax;
  for(j=0;j<nfirms;j++) {
    w[j] = (int)*(dw+j);
    z2[j] = kmax;
    justone[j] = 0;
  }
  justone[place]=1;
  binomcols = nfirms+kmax+2;  /* used in qencode() if nfirms > maxfirms_qencode */

  /* if (*debug>0)  { printf("\n"); for(j=0;j<nfirms;j++) printf("C: j= %d w= %d isex= %16.14f p_up= %16.14f  place= %2.0f \n", j, w[j], *(isex+j), *(p_up+j), *dplace );} */
  
  /* outer loop over    exit    outcomes */
  /* inner loop over investment outcomes */
  /* two_n = 2^(nfirms-1)  NOT  2^nfirms */
  /* locmask is nfirms by 2^(nfirms-1)   */
  /*    isex is nfirms by 1              */
  for(iEXIT=0;iEXIT<two_n;iEXIT++) {       /* break issued at endfor if phiH==0 --> fixed scrap --> exit handled in optimize()  */
    EXprobmask = 1.0;
    for(j=0;j<nfirms;j++) {
      EXITmask[j] = (int)*(dlocmask+iEXIT+j*two_n);     /* printf(" EXITmask[j= %2d ][iEXIT= %2d ] =  %2d \n", j, iEXIT, EXITmask[j]); */
      EXprobmask *= ( 2.0 * EXITmask[j] * *(isex+j) + 1.0 - (double)EXITmask[j] - *(isex+j) );      /* Gauss: EXprobmask = prodc(2 .* locmask[.,iEXIT] .* isex + 1 - locmask[.,iEXIT] - isex) */
    }
    if(EXprobmask>0.0 || *phiH==0.0) {  /* else skip inner loop since zero probability */
      for(i=0;i<two_n;i++) {
	probmask = 1.0;
	for(j=0;j<nfirms;j++) {                                                         /* Gauss code uses locmask[][] nfirm by two_n; referencing 2nd [] with both iEXIT and i loop index. */
	  locmask[j] = (int)*(dlocmask+i+j*two_n);                                      /* Here; we use EXITmask to be the vector for iEXIT references and locmask for i references */
	  probmask *= (2.0*(double)locmask[j]* *(p_up+j) + 1.0 - (double)locmask[j] - *(p_up+j)); /* --> faster and avoids problem of not knowing how to two_n until run time */
	  /* if (*debug>0)  printf("C: j= %d iEXIT=%d i= %d EXIT~locmask= %d %d EX~probmask = %12.10f %12.10f EX~probmask = %8.1e %8.1e \n", j, iEXIT, i, EXITmask[j], locmask[j], EXprobmask, probmask, EXprobmask, probmask); */
	}
	if(probmask>0.0) {
	  
	  if(*phiH == 0.0) {
	    for(j=0;j<nfirms;j++) d[j] =  w[j]+locmask[j];                              /* phiH==0 --> exit handled in optimize before calling calcval      */
	  }
	  else {
	    probmask *= EXprobmask;                                                     /* d = 0 for exiting firms */
	    for(j=0;j<nfirms;j++)    d[j] = (w[j]+locmask[j])*(1-EXITmask[j]);          /* phiH >0 --> random scrap --> integrate over exit here in calcval */
	    /* if (*debug>0) for(j=0;j<nfirms;j++)   printf("j=%d  w=%d  d= %d  locmask= %d  EXITmask= %d  1-EXITmask=%d \n", j,w[j],d[j],locmask[j], EXITmask[j], 1-EXITmask[j] ); */
	  } 
	  
	  *probmask_chksum += probmask;        /* should be 1 when done.  Returned and checked back in gauss */
	
	  for(j=0;j<nfirms;j++) temp[j] = 0;    
	  temp[place] = 1;

	  piksr2( nfirms, d, temp);           /* descending sort of d, with temp rearranged to match.  temp has 1 in row of place firm */

	  /* if (*debug>0) for(j=0;j<nfirms;j++)   printf("j=%d  w=%d  d= %d temp= %d   chksum= %16.14f \n", j,w[j],d[j],temp[j], *probmask_chksum);  */

	  for(j=0;j<nfirms;j++)
	    e[j] = max(0,d[j]-1);             /* Check for evaluation of value fn. at -1 */

	  if(RLG_out==0) {
	    pl1 = kmax-d[0];  if(pl1) for(j=0;j<nfirms;j++) if(d[j]) d[j] += pl1; else break;   /* all frontier firms must have exited --> move remaining firms up s.t. highest at kmax */
	    pl1 = kmax-e[0];  if(pl1) for(j=0;j<nfirms;j++) if(e[j]) e[j] += pl1; else break;
	  }
	
	  if(RLG_no_force_exit)               /* Never bump firms off lowest rung of ladder */
	    for(j=0;j<nfirms;j++)             /* This can move the focal firm into a tie with more firms, */  
	      if(d[j]==1) e[j]++;             /* but his pl1 based on temp[.,2] still valid since tied firms have same values */
	  
	  pl1=0;
	  while(temp[pl1]==0) pl1++;          /* index of place firm in sorted d */

	  qencode_e = qencode( e, nfirms, etable1, multfac1, binomcols, binomv);

	  /* if (*debug>0) { printf("C: B: e= "); for(j=0;j<nfirms;j++)  printf(" %d", e[j]);    printf(" d="); for(j=0;j<nfirms;j++)  printf(" %d", d[j]);   printf(" w="); for(j=0;j<nfirms;j++)  printf(" %d", w[j]);   printf(" EXITmask="); for(j=0;j<nfirms;j++)  printf(" %d", EXITmask[j]);  printf(" pl1= %d  qencode_e+1= %d oldvalue %8.4f \n", pl1, qencode_e+1, *(oldvalue + qencode_e*(unsigned long)nfirms + (unsigned long)pl1));  }	  */

	  if(RLG_wstar<=0 && d[0]==kmax+1) {  /* RLG: at least one firm beyond "frontier" so use delta = 1.0 --> only use  e  */
	    *valB += *(oldvalue + qencode_e*(unsigned long)nfirms + (unsigned long)pl1) * probmask;
	  }
	  else {
	    for(j=0;j<nfirms;j++)
	      if(d[j]>kmax) d[j]=kmax;      /* really only needed by PM statespace (ie, if pp_[_WSTAR]>0) but harmless otherwise since previous if already checked d==kmax+1) */

	    qencode_d = qencode( d, nfirms, etable1, multfac1, binomcols, binomv );

	    /* if (*debug>0) { printf("C: B: d= "); for(j=0;j<nfirms;j++)  printf(" %d", d[j]);  printf(" pl1= %d  qencode_d+1= %d oldvalue %8.4f \n", pl1, qencode_d+1, *(oldvalue + qencode_d*(unsigned long)nfirms + (unsigned long)pl1)); }  */

	    if(*delta>0.0) *valB += ((1.0-*delta)* *(oldvalue + qencode_d*(unsigned long)nfirms + (unsigned long)pl1 ) + *delta* *(oldvalue + qencode_e*(unsigned long)nfirms + (unsigned long)pl1 )) * probmask;
	    else           *valB +=                *(oldvalue + qencode_d*(unsigned long)nfirms + (unsigned long)pl1 ) * probmask;  /* delta = 0, so only use (1-delta) part, with coeff of 1 */
	  }

	  if(*phiH == 0.0)                                                                 /* +justone is the successful innovation */
	    for(j=0;j<nfirms;j++) d[j] =  w[j]+locmask[j]+justone[j];                      /* phiH==0 --> exit handled in optimize before calling calcval      */
	  else
	    for(j=0;j<nfirms;j++) d[j] = (w[j]+locmask[j]+justone[j])*(1-EXITmask[j]);     /* phiH >0 --> random scrap --> integrate over exit here in calcval */
	  
	  for(j=0;j<nfirms;j++) temp[j] = 0;    
	  temp[place] = 1;

	  /* if (*debug>0 && 0) {  printf("C: pre-A: d= "); for(j=0;j<nfirms;j++)  printf(" %d", d[j]);  printf("\n"); }  */

	  piksr2( nfirms, d, temp);           /* descending sort of d, with temp rearranged to match.  temp has 1 in row of place firm */
	  
	  for(j=0;j<nfirms;j++)
	    e[j] = max(0,d[j]-1);             /* Check for evaluation of value fn. at -1 */
	  
	  if(RLG_out==0) {
	    pl1 = kmax-d[0];  if(pl1) for(j=0;j<nfirms;j++) if(d[j]) d[j] += pl1; else break;   /* all frontier firms must have exited --> move remaining firms up s.t. highest at kmax */
	    pl1 = kmax-e[0];  if(pl1) for(j=0;j<nfirms;j++) if(e[j]) e[j] += pl1; else break;
	  }

	  if(RLG_no_force_exit)               /* Never bump firms off lowest rung of ladder */
	    for(j=0;j<nfirms;j++)             /* This can move the focal firm into a tie with more firms, */  
	      if(d[j]==1) e[j]++;             /* but his pl1 based on temp[.,2] still valid since tied firms have same values */
	  
	  pl1=0;
	  while(temp[pl1]==0) pl1++;          /* index of place firm in sorted d */

	  qencode_e = qencode( e, nfirms, etable1, multfac1, binomcols, binomv);

	  if(RLG_wstar<=0 && d[0]==kmax+1) {  /* RLG: at least one firm beyond "frontier" so use delta = 1.0 --> only use  e  */
	    *valA +=   /* REMOVED: (1-delta)*() + delta*  */   *(oldvalue + qencode_e*(unsigned long)nfirms + (unsigned long)pl1 ) * probmask;

	    /* if (*debug>0) { printf("C: A: e= "); for(j=0;j<nfirms;j++)  printf(" %d", d[j]); printf(" pl1= %d  qencode_e+1= %d oldvalue %8.4f \n", pl1, qencode_e+1, *(oldvalue + qencode_e*(unsigned long)nfirms + (unsigned long)pl1));  } */

	  }
	  else {
	    for(j=0;j<nfirms;j++)
	      if(d[j]>kmax) d[j]=kmax;      /* really only needed by PM statespace (ie, if pp_[_WSTAR]>0) but harmless otherwise since previous if already checked d==kmax+1) */
	    
	    qencode_d = qencode( d, nfirms, etable1, multfac1, binomcols, binomv);
	    
	    /* if (*debug>0) { printf("C: A: d= "); for(j=0;j<nfirms;j++)  printf(" %d", d[j]);  printf(" pl1= %d  qencode_d+1= %d oldvalue %8.4f \n", pl1, qencode_d+1, *(oldvalue + qencode_d*(unsigned long)nfirms + (unsigned long)pl1)); } */
	    
	    if(*delta>0.0) *valA += ((1.0-*delta)* *(oldvalue + qencode_d*(unsigned long)nfirms + (unsigned long)pl1 ) + *delta* *(oldvalue + qencode_e*(unsigned long)nfirms + (unsigned long)pl1 )) * probmask;
	    else           *valA +=                *(oldvalue + qencode_d*(unsigned long)nfirms + (unsigned long)pl1 ) * probmask;  /* delta = 0, so only use (1-delta) part, with coeff of 1 */
	  }
	}    /* if probmask>0 */
	
      }      /* i-loop over two_n industry investment outcome configurations */
	
    }        /* check whether EXprobmask==0 in which case inner loop skipped */
    
    if(*phiH==0)  break;   /* exit handled by optimize() before calling calcval --> ignore outer loop */
    
  }          /* loop over two_n industry exit  outcome configurations */
  /* if (*debug>0) printf("C: 1-chksum = %9.1e  valA~valB= %12.6f %12.6f\n ", 1.0-*probmask_chksum, *valA, *valB); */
  return(1);
}


