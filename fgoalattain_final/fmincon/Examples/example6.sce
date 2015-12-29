function y=_f(x)
	y=0
	for i =1:10
		y=y+(%e^x(i))*(-1)^i;
	end	
	
endfunction

for i =1:10
	x1(i)=-2;
	x2(i)=2;
end	

[xopt,fopt,exitflag,output,zl,zu] = fminbnd(_f,x1,x2)
	
