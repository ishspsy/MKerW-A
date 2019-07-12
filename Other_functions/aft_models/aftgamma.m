function [pars]=aftgamma2(time,Z,status,init,minimizer,options,printresults)

   
%-----Derive Initial Values-------------------------
colZ=size(Z,2);
if isempty(init)~=1
if length(init)~=colZ+2;error('The length of the initial values vector must be equal to the number of columns of Z plus two');end
end


if isempty(init)==1
    y=log(time);
    init=Z\y;
    sigmahat=sqrt((y-Z*init)'*(y-Z*init)/(length(y)-size(Z,2)));
    init=[init' sigmahat];
    
    timenew=time(status==0 | status==1);
    statusnew=status(status==0 | status==1);
    aa=gamfit(timenew,0.05,statusnew);    
    plusinitshape=sqrt(1/aa(1));
    minusinitshape=-sqrt(1/aa(1));
    init=[init plusinitshape];
end

if isempty(init)~=1
     minusinitshape=-init(end);
end

    %-----End of derivation of Initial Values-----------

    
%fwpdf=@(w,d)  abs(d)./(gamma(d.^(-2))) .* (d.^(-2).*exp(d.*w)).^(d.^(-2))  .*(exp(-exp(d.*w).*d.^(-2)));
 %logfwpdf=@(w,d) (log(abs(d))-gammaln(d.^(-2)))+(d.^(-2)).*log(d.^(-2).*exp(d.*w))+(-exp(d.*w).*d.^(-2));
 %logfwpdf=@(w,d) (log(abs(d))-gammaln(d.^(-2)))+(d.^(-2)).*(log(d.^(-2))+(d.*w))+(-exp(d.*w).*d.^(-2));
 logfwpdf=@(w,d) (log(abs(d))-gammaln(d.^(-2)))+(d.^(-2)).*(-2.*log(abs(d))+(d.*w))+(-exp(d.*w).*d.^(-2));
 
 
 Sw=@(w,d)  (d<0).*(gammainc(d.^(-2).*exp(d.*w),d.^(-2)))+...
            (d>0).*(1-((gammainc(d.^(-2).*exp(d.*w),d.^(-2)))));
   
    
 
    %---Just derived finished with the logpdf and survival functions-----
    %--Now build the likelihood:---------------------    
 %            logL=@(g) -sum(status.*(log( 1/g(colZ+1).*fwpdf((log(time)-Z*g(1:colZ)')./g(colZ+1),g(colZ+2)) )  )+...
 %                         (1-status).*log(Sw((log(time)-Z*g(1:colZ)')./g(colZ+1),g(colZ+2))));    
             logL=@(g) -sum((status==0).*(log(1/g(colZ+1))+logfwpdf((log(time)-Z*g(1:colZ)')./g(colZ+1),g(colZ+2)))+...
                            (status==1).*log(Sw((log(time)-Z*g(1:colZ)')./g(colZ+1),g(colZ+2)))+...
                           (status==-1).*log(1-Sw((log(time)-Z*g(1:colZ)')./g(colZ+1),g(colZ+2))));    
 
             lb=-Inf.*ones(1,colZ);
             ub=Inf.*ones(1,colZ);
            if isempty(minimizer)==1;minimizer=1;end  
            if minimizer==1 || isempty(minimizer)==1
                
                  if isempty(options)==1;options = optimset('MaxFunEvals',8000,'MaxIter',8000);end
                  [pars gval exitflag]=fminsearchbnd(logL,init,[lb 0 -Inf],[ub Inf Inf],options);
                  
                  if exitflag~=1;
                        init(end)=minusinitshape;
                        [pars gval exitflag]=fminsearchbnd(logL,init,[lb 0 -Inf],[ub Inf Inf]);
                  end
            elseif minimizer==2
                  if isempty(options)==1;options = optimset('LargeScale','off','Algorithm','sqp','MaxIter',8000,'MaxFunEvals',8000);end
                  [pars gval exitflag , ~ , ~, ~, HESSIAN] = fmincon(logL,init,[],[],[],[],[lb 0 -Inf],[ub Inf Inf],[],options);
                  
                  if exitflag<1
                        init(end)=minusinitshape;
                        [pars,gval,exitflag] = fmincon(logL,init,[],[],[],[],[lb 0 -Inf],[ub Inf Inf],[],options);
                  end
                      
                      
            elseif minimizer==3
                  if isempty(options)==1;options = optimset('LargeScale','off','Algorithm','interior-point','MaxIter',8000,'MaxFunEvals',8000);end
                  [pars gval exitflag , ~ , ~, ~, HESSIAN]  = fmincon(logL,init,[],[],[],[],[lb 0 -Inf],[ub Inf Inf],[],options);
                  
                  
                  if exitflag<1
                        init(end)=minusinitshape;
                        [pars,gval,exitflag] = fmincon(logL,init,[],[],[],[],[lb 0 -Inf],[ub Inf Inf],[],options);
                  end
            else error('Input argument "minimizer" must be set to 1 or 2 or 3 or []');
            end
            %Obtain the hessian use John D'errico's function as they are more accurate.            
 

end