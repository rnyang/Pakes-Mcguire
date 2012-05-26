function [retp1,retp2,retp3] = acronym(et,it)
    load pmg.mat
    load workspace.mat
% output file prefix %
    if et==COMPETITION
       s='';
       etype='Markov Perfect Nash Equilibrium';
    elseif et==MONOPOLY
       s='m';
       etype='Monopoly';
    elseif et==PLANNER
       s='s';
       etype='Social Planner';
    end
    if it==QUALITY
       s=[s,'b'];
       itype='quality';
    elseif it==COST
       s=[s,'c'];
       itype='marginal cost';
    elseif it==CAPACITY
       s=[s,'p'];
       itype='capacity';
    end
    if BatchJob>0  
        s = ['BJ' , num2str(BatchJob) , '_' , s];  
    end
    retp1 = s;
    retp2 = etype;
    retp3 = itype;
end % prefix %