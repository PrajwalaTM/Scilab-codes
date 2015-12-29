function y=_f(x)
	y=0
	for i =1:100
		y=y+x(i)^2;
	end	
	
endfunction

for i =1:100
	x1(i)=-1;
	x2(i)=1;
end	

[xopt,fopt,exitflag,output,zl,zu] = fminbnd(_f,x1,x2)

