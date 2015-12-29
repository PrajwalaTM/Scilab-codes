mode(1)
//
// Demo of fminunc.sci
//

halt()   // Press return to continue
 
//Find x in R^2 such that it minimizes rosenbrock function
//f = 100*(x2 - x1^2)^2 + (1-x1)^2
halt()   // Press return to continue
 
//Objective function to be minimised
function y= f(x)
y= 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
endfunction
halt()   // Press return to continue
 
//Starting point
x0=[-1,2];
halt()   // Press return to continue
 
//Options
options=list("MaxIter", [1500], "CpuTime", [500], "Gradient", "ON", "Hessian", "ON");
halt()   // Press return to continue
 
//Gradient of objective function
function y= fGrad(x)
y= [-400*x(1)*x(2) + 400*x(1)^3 + 2*x(1)-2, 200*(x(2)-x(1)^2)];
endfunction
halt()   // Press return to continue
 
//Hessian of Objective Function
function y= fHess(x)
y= [1200*x(1)^2- 400*x(2) + 2, -400*x(1);-400*x(1), 200 ];
endfunction
halt()   // Press return to continue
 
//Calling the Ipopt
[xopt,fopt,exitflag,output,gradient,hessian]=fminunc(f,x0,options,fGrad,fHess)
halt()   // Press return to continue
 
halt()   // Press return to continue
 
halt()   // Press return to continue
 
//Find x in R^2 such that the below function is minimum
//f = x1^2 + x2^2
halt()   // Press return to continue
 
//Objective function to be minimised
function y= f(x)
y= x(1)^2 + x(2)^2;
endfunction
halt()   // Press return to continue
 
//Starting point
x0=[2,1];
halt()   // Press return to continue
 
//Calling the Ipopt
[xopt,fopt]=fminunc(f,x0)
halt()   // Press return to continue
 
halt()   // Press return to continue
 
The below Problem is an Unbounded problem:
//Find x in R^2 such that the below function is minimum
//f = - x1^2 - x2^2
halt()   // Press return to continue
 
//Objective function to be minimised
function y= f(x)
y= -x(1)^2 - x(2)^2;
endfunction
halt()   // Press return to continue
 
//Starting point
x0=[2,1];
//Options
options=list("MaxIter", [1500], "CpuTime", [500], "Gradient", "ON", "Hessian", "ON");
halt()   // Press return to continue
 
//Gradient of objective function
function y= fGrad(x)
y= [-2*x(1),-2*x(2)];
endfunction
halt()   // Press return to continue
 
//Hessian of Objective Function
function y= fHess(x)
y= [-2,0;0,-2];
endfunction
halt()   // Press return to continue
 
//Calling the Ipopt
[xopt,fopt,exitflag,output,gradient,hessian]=fminunc(f,x0,options,fGrad,fHess)
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
