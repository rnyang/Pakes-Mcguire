function retp = default
    load pmg.mat
    p=zeros(NPARMS,1);
    p(MAX_FIRMS)= 3;
    p(START_FIRMS)= 1;
    p(EQL_TYPE)= COMPETITION;
    p(IND_TYPE)= QUALITY;
    p(ENTRY_TYPE)= RAN_ENTRY;
    p(ENTRY_LOW)= 0.3;
    p(ENTRY_HIGH)= 0.5;
    p(ENTRY_SUNK)= 0.2;
    p(ENTRY_AT)= 5;
    p(BETA)= 0.925;
    p(DELTA)= 0.7;

    p(INV_MULT)= 3;
    p(INV_COST)= 1;
    p(MC)= 5;
    p(MKT_SIZE)= 5;
    p(KMAX)= 15;
    p(WSTAR)= 12;           % 12 for Rand94 ? %
    p(INTERCEPT)= 3;
    p(FIXED_COST)= 0.0;
    p(GAMMA)= 1;
    p(TAU)= 0.1;

    p(RLG_OUTGOOD)=  1;
    p(RLG_ALPHA)=    1;
    p(RLG_WSCALE)=   1;    %  3 for Rand94 ? %
    p(RLG_WSHIFT)=   0;    % -7 for Rand94 ? %
    p(RLG_SH_CAP)=   1;    % .55, .65 for Rand94 %
    p(RLG_INV)=   0;       % rate at which investment efficiency increases in w %
    p(RLG_y)=  10;         % -alpha*log(RLG_y - price) %
    p(RLG_LEAP)=  0;       % prob entrant starts at KMAX %
    p(RLG_MAXP)=  1.5;     % maximimum monopolist price =  RLG_MAXP*max_duop_price in competitive equilibrium when RLG_y =  0 (ie, linear).  default =  9999 %
    p(RLG_FOC) = 0;        % if 0, use FOC method in  eql_ma_foc.g (where this value is hardcoded)  ONLY if non-FOC fails to converge.   If >0, do not try non-FOC.  If <0, never use FOC %
    p(RLG_MIN_INNOV1)=0;   % minimum innovation rate by frontier firms.  Only used if RLG_WSTAR <= 0 (i.e., not using wstar of original PM) %

    p(SCRAP_VAL)   =  5;   % lower bound of scrap %
    p(RLG_SCRAP_HI)=  6;   % upper bound of scrap value, 0 for fixed scrap %

    p(PROFIT_DONE)= 0;
    p(EQL_DONE)= 0;
%     p(ACTIVE_CFG)= 'default';
    p(ACTIVE_CFG)= 0;
    retp = p;
end % default %