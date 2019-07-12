function [pars]=aftexponential2(time,Z,status,init,minimizer,options,printresults)
  
%-----Derive Initial Values-------------------------
colZ=size(Z,2);
if isempty(init)~=1
if length(init)~=colZ;error('The length of the initial values vector must be equal to the number of columns of Z');end
end

if isempty(init)==1
    y=log(time);
    init=Z\y;
    init=init';
end

    %-----End of derivation of Initial Values-----------

            logL=@(g) -sum((status==0).*(log( 1/1.*evpdf((log(time)-Z*g(1:colZ)')./1,0,1) )  )+...  %event
                           (status==1).*log(1-evcdf((log(time)-Z*g(1:colZ)')./1,0,1))+...           %right cens
                          (status==-1).*log(evcdf((log(time)-Z*g(1:colZ)')./1,0,1)));               %left cens
            lb=-Inf.*ones(1,colZ);
            ub=Inf.*ones(1,colZ);
            if isempty(minimizer)==1;minimizer=1;end
            if minimizer==1
                  if isempty(options)==1;options = optimset('MaxFunEvals',8000,'MaxIter',8000);end
                  [pars gval exitflag]=fminsearchbnd(logL,init,lb,ub,options);
            elseif minimizer==2
                  if isempty(options)==1;options = optimset('LargeScale','off','Algorithm','sqp');end
                  [pars gval exitflag , ~ , ~, ~, HESSIAN] = fmincon(logL,init,[],[],[],[],lb,ub,[],options);
            elseif minimizer==3
                  if isempty(options)==1;options = optimset('LargeScale','off','Algorithm','interior-point');end
                  [pars gval exitflag , ~ , ~, ~, HESSIAN] = fmincon(logL,init,[],[],[],[],lb,ub,[],options);
            else error('Input argument minimizer must be 1 or 2 or 3 or []');
            end
            %Obtain the hessian use John D'errico's function as they are more accurate.            
            

end