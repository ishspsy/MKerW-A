function [pars covpars SE CI Zscores pvalues gval exitflag hess fminconhess]=aftweibull(time,Z,status,init,minimizer,options,printresults)


%-----Derive Initial Values-------------------------


colZ=size(Z,2);
if isempty(init)~=1
if length(init)~=colZ+1;error('The length of the initial values vector must be equal to the number of columns of Z plus one');end
end


if isempty(init)==1
    y=log(time);
    init=Z\y;
    sigmahat=sqrt((y-Z*init)'*(y-Z*init)/(length(y)-size(Z,2)));
    init=[init' sigmahat];
end


%-----End of derivation of Initial Values-----------

            logL=@(g) -sum((status==0).*(log( 1/g(colZ+1).*evpdf((log(time)-Z*g(1:colZ)')./g(colZ+1),0,1) )  )+...
                           (status==1).*log(1-evcdf((log(time)-Z*g(1:colZ)')./g(colZ+1),0,1))+...
                          (status==-1).*log(evcdf((log(time)-Z*g(1:colZ)')./g(colZ+1),0,1)) );    
            
            
            lb=-Inf.*ones(1,colZ);
            ub=Inf.*ones(1,colZ);
            if isempty(minimizer)==1;minimizer=1;end
            if minimizer==1
                  if isempty(options)==1;options = optimset('MaxFunEvals',8000,'MaxIter',8000);end
                  [pars gval exitflag]=fminsearchbnd(logL,init,[lb 0],[ub Inf],options);
            elseif minimizer==2
                  if isempty(options)==1;options = optimset('LargeScale','off','Algorithm','sqp','MaxFunEvals',8000,'MaxIter',8000);end
                  [pars gval exitflag , ~ , ~, ~, HESSIAN] = fmincon(logL,init,[],[],[],[],[lb 0],[ub Inf],[],options);
            elseif minimizer==3
                  if isempty(options)==1;options = optimset('LargeScale','off','Algorithm','interior-point','MaxFunEvals',8000,'MaxIter',8000);end
                  [pars gval exitflag , ~ , ~, ~, HESSIAN] = fmincon(logL,init,[],[],[],[],[lb 0],[ub Inf],[],options);
            else error('Input argument minimizer must be 1 or 2 or 3 or []');
            end
            %Obtain the hessian: use John D'errico's function.     
            
        
            hess = hessian(logL,pars);
            
            fminconhess=0;
            if minimizer>1
                if sum(sum(isnan(hess)))~=0 || sum(sum(imag(hess)))~=0
                hess=HESSIAN;fminconhess=1;
                end
            end
            covpars=inv(hess);
            SE=sqrt(diag(covpars));
            CI=[pars'-1.96.*SE pars'+1.96.*SE];
            CI(end,:)=[pars(end)/exp(1.96*SE(end)/pars(end)) pars(end)*exp(1.96*SE(end)/pars(end))];
            
            Zscores=pars'./SE;
            pvalues=2*normcdf(-abs(Zscores(1:end-1)));

            
            


if printresults==1   
      print_weibull(time, status, Z,pars, SE, CI, Zscores, pvalues, minimizer,exitflag,gval,fminconhess);
end


end