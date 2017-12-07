
function printmyrow(myrow)
%-------PRINT THE RESULTS-----------------------

%Print the intercept stuff:    
if myrow(1)<0
    if abs(myrow(1))<10;
        fprintf('%4.5f     ',myrow(1));
    elseif 10<abs(myrow(1)) && abs(myrow(1))<100;
        fprintf('%4.4f     ',myrow(1));
    elseif 100<abs(myrow(1)) && abs(myrow(1))<1000;
        fprintf('%4.3f     ',myrow(1));
    elseif 1000<abs(myrow(1)) && abs(myrow(1))<10000;
        fprintf('%4.2f     ',myrow(1));
    else
        fprintf('%4.1f     ',myrow(1));
    end
    
else
    
    if abs(myrow(1))<10;
        fprintf(' %4.5f     ',myrow(1));
    elseif 10<abs(myrow(1)) && abs(myrow(1))<100;
        fprintf(' %4.4f     ',myrow(1));
    elseif 100<abs(myrow(1)) && abs(myrow(1))<1000;
        fprintf(' %4.3f     ',myrow(1));
    elseif 1000<abs(myrow(1)) && abs(myrow(1))<10000;
        fprintf(' %4.2f     ',myrow(1));
    else
        fprintf(' %4.1f     ',myrow(1));
    end
    
    
end
    
    



if myrow(2)<0
    if abs(myrow(2))<10;
        fprintf('%4.5f     ',myrow(2));
    elseif 10<abs(myrow(2)) && abs(myrow(2))<100;
        fprintf('%4.4f     ',myrow(2));
    elseif 100<abs(myrow(2)) && abs(myrow(2))<1000;
        fprintf('%4.3f     ',myrow(2));
    elseif 1000<abs(myrow(2)) && abs(myrow(2))<10000;
        fprintf('%4.2f     ',myrow(2));
    else
        fprintf('%4.1f     ',myrow(2));
    end
    
else
    
    if abs(myrow(2))<10;
        fprintf(' %4.5f     ',myrow(2));
    elseif 10<abs(myrow(2)) && abs(myrow(2))<100;
        fprintf(' %4.4f     ',myrow(2));
    elseif 100<abs(myrow(2)) && abs(myrow(2))<1000;
        fprintf(' %4.3f     ',myrow(2));
    elseif 1000<abs(myrow(2)) && abs(myrow(2))<10000;
        fprintf(' %4.2f     ',myrow(2));
    else
        fprintf(' %4.1f     ',myrow(2));
    end
        
end
    
    
if myrow(3)<0
    if abs(myrow(3))<10;
        fprintf('%4.5f     ',myrow(3));
    elseif 10<abs(myrow(3)) && abs(myrow(3))<100;
        fprintf('%4.4f     ',myrow(3));
    elseif 100<abs(myrow(3)) && abs(myrow(3))<1000;
        fprintf('%4.3f     ',myrow(3));
    elseif 1000<abs(myrow(3)) && abs(myrow(3))<10000;
        fprintf('%4.2f     ',myrow(3));
    else
        fprintf('%4.1f     ',myrow(3));
    end
else
    
    if abs(myrow(3))<10;
        fprintf(' %4.5f     ',myrow(3));
    elseif 10<abs(myrow(3)) && abs(myrow(3))<100;
        fprintf(' %4.4f     ',myrow(3));
    elseif 100<abs(myrow(3)) && abs(myrow(3))<1000;
        fprintf(' %4.3f     ',myrow(3));
    elseif 1000<abs(myrow(3)) && abs(myrow(3))<10000;
        fprintf(' %4.2f     ',myrow(3));
    else
        fprintf(' %4.1f     ',myrow(3));
    end
    
end
    
    

    
    
if myrow(4)<0
    if abs(myrow(4))<10;
        fprintf('%4.5f     ',myrow(4));
    elseif 10<abs(myrow(4)) && abs(myrow(4))<100;
        fprintf('%4.4f     ',myrow(4));
    elseif 100<abs(myrow(4)) && abs(myrow(4))<1000;
        fprintf('%4.3f     ',myrow(4));
    elseif 1000<abs(myrow(4)) && abs(myrow(4))<10000;
        fprintf('%4.2f     ',myrow(4));
    else
        fprintf('%4.1f     ',myrow(4));
    end
else
    if abs(myrow(4))<10;
        fprintf(' %4.5f     ',myrow(4));
    elseif 10<abs(myrow(4)) && abs(myrow(4))<100;
        fprintf(' %4.4f     ',myrow(4));
    elseif 100<abs(myrow(4)) && abs(myrow(4))<1000;
        fprintf(' %4.3f     ',myrow(4));
    elseif 1000<abs(myrow(4)) && abs(myrow(4))<10000;
        fprintf(' %4.2f     ',myrow(4));
    else
        fprintf(' %4.1f     ',myrow(4));
    end
    
end
    
    if length(myrow)>4
    
if myrow(5)<0
    if abs(myrow(5))<10;
        fprintf('%4.5f     ',myrow(5));
    elseif 10<abs(myrow(5)) && abs(myrow(5))<100;
        fprintf('%4.4f     ',myrow(5));
    elseif 100<abs(myrow(5)) && abs(myrow(5))<1000;
        fprintf('%4.3f     ',myrow(5));
    elseif 1000<abs(myrow(5)) && abs(myrow(5))<10000;
        fprintf('%4.2f     ',myrow(5));
    else
        fprintf('%4.1f     ',myrow(5));
    end
else
    if abs(myrow(5))<10;
        fprintf(' %4.5f     ',myrow(5));
    elseif 10<abs(myrow(5)) && abs(myrow(5))<100;
        fprintf(' %4.4f     ',myrow(5));
    elseif 100<abs(myrow(5)) && abs(myrow(5))<1000;
        fprintf(' %4.3f     ',myrow(5));
    elseif 1000<abs(myrow(5)) && abs(myrow(5))<10000;
        fprintf(' %4.2f     ',myrow(5));
    else
        fprintf('%4.1f     ',myrow(5));
    end
end
    
    
    


if myrow(6)<0
    if myrow(6)<10;
        fprintf('%4.5f     ',myrow(6));
    elseif 10<myrow(6) && myrow(6)<100;
        fprintf('%4.4f     ',myrow(6));
    elseif 100<myrow(6) && myrow(6)<1000;
        fprintf('%4.3f     ',myrow(6));
    elseif 1000<myrow(6) && myrow(6)<10000;
        fprintf('%4.2f     ',myrow(6));
    else
        fprintf('%4.1f     ',myrow(6));
    end
else
    if myrow(6)<10;
        fprintf(' %4.5f     ',myrow(6));
    elseif 10<myrow(6) && myrow(6)<100;
        fprintf(' %4.4f     ',myrow(6));
    elseif 100<myrow(6) && myrow(6)<1000;
        fprintf(' %4.3f     ',myrow(6));
    elseif 1000<myrow(6) && myrow(6)<10000;
        fprintf(' %4.2f     ',myrow(6));
    else
        fprintf(' %4.1f     ',myrow(6));
    end
end
    

    
    fprintf('%0.4f\n',myrow(7));
    
    
    
    
    
    end
    