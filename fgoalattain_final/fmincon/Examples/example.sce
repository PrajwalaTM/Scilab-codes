  //      //Objective function to be minimised
      function y= f(x)
     	     y= -x(1)^2 - x(2)^2;
        endfunction
  //
  //	  //Starting point  
        x0=[2,1];
  //     //Options
       options=list("MaxIter", [1500], "CpuTime", [500], "Gradient", "ON", "Hessian", "ON");
  //
  //	 //Gradient of objective function
       function y= fGrad(x)
     	     y= [-2*x(1),-2*x(2)];
       endfunction
  //
  //	 //Hessian of Objective Function
       function y= fHess(x)
     	     y= [-2,0;0,-2];
       endfunction
  //
  //      //Calling the Ipopt  
        [xopt,fopt,exitflag,output,gradient,hessian]=fminunc(f,x0,options,fGrad,fHess)
