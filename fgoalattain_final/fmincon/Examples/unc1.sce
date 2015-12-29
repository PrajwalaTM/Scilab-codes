//function [c]=item(x)
//c(1)=x(1)^2 + x(2)^2 -1;
//c(2)=x(1)^3+x(2)^3;
//endfunction
//pt=[1,2];
//[a,b]=numderivative(item,pt)




//fun = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
// x(1) + 2x(2) <=1
// 2x(1) + x(2) = 1
// 0<=x(1)<=2
// 0<=x(2)<=2
//(x(1)-1/3)^2 + (x(2)-1/3)^2 - (1/3)^2<=0
//-x(1)^2<=0
// Gradient= [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1)), 200*(x(2)-x(1)^2)];
// Hessian= [1200*x(1)^2, -400*x(1);-400*x(1), 200 ];

function y=f(x)
	y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

x0=[-1,2];

options=list("MaxIter", [1500], "CpuTime", [500], "Gradient", "ON", "Hessian", "ON");

function y= funG(x)
	y= [-400*x(1)*x(2) + 400*x(1)^3 + 2*x(1)-2, 200*(x(2)-x(1)^2)];
endfunction

function y= funH(x)
	y= [1200*x(1)^2- 400*x(2) + 2, -400*x(1);-400*x(1), 200 ];
endfunction

[x,fval,exitflag,output,grad,hessian] =fminunc(f, x0,options,funG,funH)
