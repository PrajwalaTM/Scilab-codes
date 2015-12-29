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

function y=_f(x)
	y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

x0=[1,2];
A=[5,7];
b=[1];
Aeq=[2,1];
beq=1
no_nlic=[1];
lb=[0,2];
ub=[1,6];

function [y]=_nlc(x)
	y(1)=(x(1)-1/3)^2 + (x(2)-1/3)^2 - (1/3)^2
	y(2)=-x(1)^2
	y(3)=x(2)
endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "HessObj", "ON","GradCon", "OFF");
function y= _funG(x)
	y= [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1)), 200*(x(2)-x(1)^2)];
endfunction

function y= _funH(x)
	y= [1200*x(1)^2, -400*x(1);-400*x(1), 200 ];
endfunction

function [y]= _conG(x)
	y(1)=2*(x(1)-1/3)
	y(2)=2*(x(1)-1/3)
	y(3)=-2*x(1)
	y(4)=x(2)
	y(5)=0
	y(6)=1
endfunction

function y= _conH(x)
	y(1)=[2,0;0,2];
	y(2)=[-2,0;0,0]
endfunction
//[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(_f, x0,[],[] ,Aeq ,beq,[],[], no_nlic, _nlc, options, _funG, _funH);
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(_f, x0,[] ,[])
