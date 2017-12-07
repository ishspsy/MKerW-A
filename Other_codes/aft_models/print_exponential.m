function print_exponential(time, status, Z,pars, SE, CI, Zscores, pvalues, minimizer,exitflag,gval,fminconhess)


%-------PRINT THE RESULTS-----------------------
disp('Results:')
disp(' ')
tr=repmat('-',1,100);
trr=repmat('=',1,100);
disp(tr)
fprintf(' Sample size:  %0.0f\n ',length(time));

fprintf('Level of right censoring:  %0.4f%%\n ',length(find(status==1))/length(status)*100);
fprintf('Level of left censoring:  %0.4f%%\n ',length(find(status==-1))/length(status)*100);

fprintf('Distribution: Exponential\n');
if minimizer==1
    fprintf(' Minimizer:  fminsearchbnd\n')
else
    fprintf(' Minimizer:  fmincon\n')
end

if fminconhess==1
fprintf('Hessian obtained by fmincon: results may be inaccurate\n');
end



fprintf(' Exitflag:  %0.0f\n ',exitflag)
disp(trr)
disp('Analysis of Parameter Estimates:')

disp(trr)
disp('Parameters    Estimates        SE     95% Confidence Intervals    Z Scores       Chi2      p-values')
disp(tr)




if isempty(find((Z(:,1)==ones(length(Z(:,1)),1))==0, 1))~=1
    parms=[0 pars];
    CIs=[[0 0];CI];
    Zscoress=[NaN;Zscores];
    pvaluess=[NaN;pvalues];
    SEs=[0;SE];
    noint=1;
else    
 
    parms=pars;
    CIs=CI;
    Zscoress=Zscores;
    pvaluess=pvalues;
    SEs=SE;
    noint=0;

end






%---Now use printmyrow:

%noint==0 %--> there is an intercept

for i=1:(length(parms))
    if i==1 && noint==0;
        fprintf('Intercept:    ');
        myrow=[parms(i) SEs(i) CIs(i,1) CIs(i,2) Zscoress(i) Zscoress(i)^2 pvaluess(i)];
        printmyrow(myrow)
        
    elseif i==1 && noint==1;
           fprintf('Intercept:    ');
           fprintf('\n');
    else
        fprintf('g%0.0f:           ',i-1);
        myrow=[parms(i) SEs(i) CIs(i,1) CIs(i,2) Zscoress(i) Zscoress(i)^2 pvaluess(i)];
        printmyrow(myrow)
    end    
end

fprintf('\n')
disp(tr)
fprintf('Scale parameter held fixed at 1')
fprintf('\n')
fprintf('Log-likelihood:');
fprintf(' %0.4f\n ',-gval);






end
