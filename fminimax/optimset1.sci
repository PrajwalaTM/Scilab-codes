// Calling Sequence
  //   options = optimset ()
  //   options = optimset ( funname )
  //   options = optimset ( key , value )
  //   options = optimset ( key1 , value1 , key2 , value2 , ... )
  //   options = optimset ( options , key , value )
  //   options = optimset ( options , key1 , value1 , key2 , value2 , ... )
  //
function options = optimset (varargin)
	[lhs,rhs]=argn()
	if(rhs==0) then
	options=optimset_new();
	return
	else if(rhs==1) then
	method=varargin(1)
	options=optimset_method(method)
	return
	end
        else if mod(rhs,2)<>0 then
	options=varargin(1)
        t = typeof(options);
    if t<>"st" then
      errmsg = msprintf(gettext("%s: Odd number of arguments : the first argument is expected to be a struct, but is a %s"), "optimset", t);
      error(errmsg)
    end
    else
    options = optimset_new ();
    end
    if mod(rhs,2) <>2 then 
    nvar = 1;
    else
    nvar = 0;
    end
    nkeys = rhs/2;
    for i=1:nkeys
    nvar = nvar + 1;
    key = varargin(nvar);
    nvar = nvar + 1;
    prot = funcprot();
    funcprot(0);
    value = varargin(nvar);
    funcprot(prot);
    options = optimset_configure (options,key,value);
  end
endfunction
function options=optimset_new()
options = struct(...
      "MinAbsMax"   ,[], ...
      "GradObj"     , [], ...
      "GradConstr"  , [] , ...
      "MaxFunEvals" , [] ...
      );
endfunction

function options = optimset_configure ( options , key , value )
    select key
    case "MinAbsMax" then
      options.MinAbsMax = value;
    case "MaxFunEvals" then
      options.MaxFunEvals = value;
    case "GradObj" then
      options.GradObj = value;
    case "GradConstr" then
      options.GradConstr = value;
    else
      errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "optimset", key)
      error(errmsg)
    end
endfunction
function options = optimset_method ( method )
    options = optimset_new ();
    select method
    case "fminsearch" then
      options = optimset_configure ( options , "Display" , "notify" );
      options = optimset_configure ( options , "MaxFunEvals" , "200*numberofvariables" );
      options = optimset_configure ( options , "MaxIter" , "200*numberofvariables" );
      options = optimset_configure ( options , "TolFun" , 1.e-4 );
      options = optimset_configure ( options , "TolX" , 1.e-4 );

    case "fmincon" then
      options = optimset_configure ( options , "Display" , "final" );
      options = optimset_configure ( options , "TolFun" , 1.e-6 );
      options = optimset_configure ( options , "TolCon" , 1.e-6 );
      options = optimset_configure ( options , "Algorithm" , "ipopt" );
    
    case "fminimax" then
      options = optimset_configure ( options , "MinAbsMax","min(abs(max))");
      options = optimset_configure ( options , "MaxFunEvals" , "200*numberofvariables" );
      options = optimset_configure ( options , "GradConstr" , "200*numberofvariables" );
      options = optimset_configure ( options , "GradObj" , 1.e-4 );
    else
      errmsg = msprintf(gettext("%s: No default options available: the function ''%s'' does not exist on the path."), "optimset", method)
      error(errmsg)
    end
endfunction

