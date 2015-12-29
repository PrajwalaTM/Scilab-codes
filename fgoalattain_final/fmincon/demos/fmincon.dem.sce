mode(1)
//
// Demo of fmincon.sci
//

halt()   // Press return to continue
 
//Find x in R^2 such that it minimizes:
//f(x)= -x1 -x2/3
//x0=[0,0]
//constraint-1 (c1): x1 + x2 <= 2
//constraint-2 (c2): x1 + x2/4 <= 1
//constraint-3 (c3): x1 - x2 <= 2
//constraint-4 (c4): -x1/4 - x2 <= 1
//constraint-5 (c5): -x1 - x2 <= -1
//constraint-6 (c6): -x1 + x2 <= 2
//constraint-7 (c7): x1 + x2 = 2
halt()   // Press return to continue
 
//Objective function to be minimised
function y=f(x)
y=-x(1)-x(2)/3;
endfunction
halt()   // Press return to continue
 
//Starting point, linear constraints and variable bounds
x0=[0 , 0 , 0];
A=[1,1 ; 1,1/4 ; 1,-1 ; -1/4,-1 ; -1,-1 ; -1,1];
b=[2;1;2;1;-1;2];
Aeq=[1,1];
beq=[2];
lb=[];
ub=[];
halt()   // Press return to continue
 
//Options
options=list("GradObj", "ON", "Hessian", "ON","GradCon", "OFF");
halt()   // Press return to continue
 
//Gradient of objective function
function y= fGrad(x)
y= [-1,-1/3];
endfunction
halt()   // Press return to continue
 
//Hessian of lagrangian
function y= lHess(x,obj,lambda)
y= obj*[0,0;0,0]
endfunction
halt()   // Press return to continue
 
//Calling IPopt
[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,lHess)
halt()   // Press return to continue
 
halt()   // Press return to continue
 
//Find x in R^3 such that it minimizes:
//f(x)= x1*x2 + x2*x3
//x0=[0.1 , 0.1 , 0.1]
//constraint-1 (c1): x1^2 - x2^2 + x3^2 <= 2
//constraint-2 (c2): x1^2 + x2^2 + x3^2 <= 10
halt()   // Press return to continue
 
//Objective function to be minimised
function y=f(x)
y=x(1)*x(2)+x(2)*x(3);
endfunction
halt()   // Press return to continue
 
//Starting point, linear constraints and variable bounds
x0=[0.1 , 0.1 , 0.1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];
halt()   // Press return to continue
 
//Nonlinear constraints
function [c,ceq]=nlc(x)
c = [x(1)^2 - x(2)^2 + x(3)^2 - 2 , x(1)^2 + x(2)^2 + x(3)^2 - 10];
ceq = [];
endfunction
halt()   // Press return to continue
 
//Options
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "Hessian", "ON","GradCon", "ON");
halt()   // Press return to continue
 
//Gradient of objective function
function y= fGrad(x)
y= [x(2),x(1)+x(3),x(2)];
endfunction
halt()   // Press return to continue
 
//Hessian of the Lagrange Function
function y= lHess(x,obj,lambda)
y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,-2,0;0,0,2] + lambda(2)*[2,0,0;0,2,0;0,0,2]
endfunction
halt()   // Press return to continue
 
//Gradient of Non-Linear Constraints
function [cg,ceqg] = cGrad(x)
cg=[2*x(1) , -2*x(2) , 2*x(3) ; 2*x(1) , 2*x(2) , 2*x(3)];
ceqg=[];
endfunction
halt()   // Press return to continue
 
//Calling IPopt
[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,lHess,cGrad)
halt()   // Press return to continue
 
halt()   // Press return to continue
 
The below Problem is an Unbounded problem:
//Find x in R^3 such that it minimizes:
//f(x)= -(x1^2 + x2^2 + x3^2)
//x0=[0.1 , 0.1 , 0.1]
//  x1 <= 0
//  x2 <= 0
//  x3 <= 0
halt()   // Press return to continue
 
//Objective function to be minimised
function y=f(x)
y=-(x(1)^2+x(2)^2+x(3)^2);
endfunction
halt()   // Press return to continue
 
//Starting point, linear constraints and variable bounds
x0=[0.1 , 0.1 , 0.1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[0,0,0];
halt()   // Press return to continue
 
//Options
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "Hessian", "OFF","GradCon", "OFF");
halt()   // Press return to continue
 
//Calling IPopt
[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,[],options)
halt()   // Press return to continue
 
halt()   // Press return to continue
 
The below Problem is an Infeasible problem:
//Find x in R^3 such that in minimizes:
//f(x)=x1*x2 + x2*x3
//x0=[1,1,1]
//constraint-1 (c1): x1^2 <= 1
//constraint-2 (c2): x1^2 + x2^2 <= 1
//constraint-3 (c3): x3^2 <= 1
//constraint-4 (c4): x1^3 = 0.5
//constraint-5 (c5): x2^2 + x3^2 = 0.75
// 0 <= x1 <=0.6
// 0.2 <= x2 <= inf
// -inf <= x3 <= 1
halt()   // Press return to continue
 
//Objective function to be minimised
function y=f(x)
y=x(1)*x(2)+x(2)*x(3);
endfunction
halt()   // Press return to continue
 
//Starting point, linear constraints and variable bounds
x0=[1,1,1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0.2,-%inf];
ub=[0.6 %inf,1];
halt()   // Press return to continue
 
//Nonlinear constraints
function [c,ceq]=nlc(x)
c=[x(1)^2-1,x(1)^2+x(2)^2-1,x(3)^2-1];
ceq=[x(1)^3-0.5,x(2)^2+x(3)^2-0.75];
endfunction
halt()   // Press return to continue
 
//Options
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "Hessian", "ON","GradCon", "ON");
halt()   // Press return to continue
 
//Gradient of objective function
function y= fGrad(x)
y= [x(2),x(1)+x(3),x(2)];
endfunction
halt()   // Press return to continue
 
//Hessian of the Lagrange Function
function y= lHess(x,obj,lambda)
y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,0,0;0,0,0] + lambda(2)*[2,0,0;0,2,0;0,0,0] +lambda(3)*[0,0,0;0,0,0;0,0,2] + lambda(4)*[6*x(1),0,0;0,0,0;0,0,0] + lambda(5)*[0,0,0;0,2,0;0,0,2];
endfunction
halt()   // Press return to continue
 
//Gradient of Non-Linear Constraints
function [cg,ceqg] = cGrad(x)
cg = [2*x(1),0,0;2*x(1),2*x(2),0;0,0,2*x(3)];
ceqg = [3*x(1)^2,0,0;0,2*x(2),2*x(3)];
endfunction
halt()   // Press return to continue
 
//Calling IPopt
[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,lHess,cGrad)
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
