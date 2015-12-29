mode(1)
//
// Demo of fminbnd.sci
//

halt()   // Press return to continue
 
//Find x in R^6 such that it minimizes:
//f(x)= sin(x1) + sin(x2) + sin(x3) + sin(x4) + sin(x5) + sin(x6)
//-2 <= x1,x2,x3,x4,x5,x6 <= 2
halt()   // Press return to continue
 
//Objective function to be minimised
function y=f(x)
y=0
for i =1:6
y=y+sin(x(i));
end
endfunction
halt()   // Press return to continue
 
//Variable bounds
x1 = [-2, -2, -2, -2, -2, -2];
x2 = [2, 2, 2, 2, 2, 2];
halt()   // Press return to continue
 
//Options
options=list("MaxIter",[1500],"CpuTime", [100],"TolX",[1e-6])
halt()   // Press return to continue
 
//Calling IPopt
[x,fval,exitflag,output,zl,zu] =fminbnd(f, x1, x2, options)
halt()   // Press return to continue
 
halt()   // Press return to continue
 
//Find x in R such that it minimizes:
//f(x)= 1/x^2
//0 <= x <= 1000
halt()   // Press return to continue
 
//Objective function to be minimised
function y=f(x)
y=1/x^2
endfunction
halt()   // Press return to continue
 
//Variable bounds
x1 = [0];
x2 = [1000];
halt()   // Press return to continue
 
//Calling IPopt
[x,fval,exitflag,output,zl,zu] =fminbnd(f, x1, x2)
halt()   // Press return to continue
 
halt()   // Press return to continue
 
The below Problem is an Unbounded problem:
//Find x in R^2 such that it minimizes:
//f(x)= -[(x1-1)^2 + (x2-1)^2]
//-inf <= x1,x2 <= inf
halt()   // Press return to continue
 
//Objective function to be minimised
function y=f(x)
y=-((x(1)-1)^2+(x(2)-1)^2);
endfunction
halt()   // Press return to continue
 
//Variable bounds
x1 = [-%inf , -%inf];
x2 = [];
halt()   // Press return to continue
 
//Options
options=list("MaxIter",[1500],"CpuTime", [100],"TolX",[1e-6])
halt()   // Press return to continue
 
//Calling IPopt
[x,fval,exitflag,output,zl,zu] =fminbnd(f, x1, x2, options)
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
