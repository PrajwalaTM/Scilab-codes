function options = optimset (varargin)
  //   Configures and returns an updated optimization data structure.
  //
  // Calling Sequence
  //   options = optimset ()
  //   options = optimset ( funname )
  //   options = optimset ( key , value )
  //   options = optimset ( key1 , value1 , key2 , value2 , ... )
  //   options = optimset ( options , key , value )
  //   options = optimset ( options , key1 , value1 , key2 , value2 , ... )
  //
  // Description
  //   Configures and returns an updated optimization data structure.
  //
  // Authors
  // Michael Baudin, DIGITEO, 2010

  [lhs,rhs]=argn();
  if rhs == 0 then
    options = optimset_new ();
    return
  elseif rhs==1 then
    //
    //   options = optimset ( funname )
    // If there is only one argument, it is expected to be the
    // name of a method.
    //
    method = varargin(1);
    options = optimset_method ( method );
    return
  end
  // Set the options variable
  if modulo(rhs,2)<>0 then
    //
    //   options = optimset ( options , key , value )
    //   options = optimset ( options , key1 , value1 , key2 , value2 , ... )
    // If the number of arguments is odd,
    // the first argument is expected to be a optimset struct.
    //
    options = varargin(1);
    t1 = typeof(options);
    if t1<>"st" then
      errmsg = msprintf(gettext("%s: Odd number of arguments : the first argument is expected to be a struct, but is a %s"), "optimset", t1);
      error(errmsg)
    end
  else
    //   options = optimset ( key , value )
    //   options = optimset ( key1 , value1 , key2 , value2 , ... )
    // Number of input argument is even.
    options = optimset_new ();
  end
  // Set ivar : index of input variable.
  // The variable ivar allows to make a loop over input arguments.
  if ( modulo(rhs,2)<>0 ) then
    ivar = 1;
  else
    ivar = 0;
  end
  //
  // Process key,values as pairs.
  //
  nbkeys = rhs/2;
  for i=1:nbkeys
    ivar = ivar + 1;
    key = varargin(ivar);
    ivar = ivar + 1;
    // Use funcprot to enable the set of a function into the variable "value".
    // If not, a warning message is triggered, when a double value 
    // is stored into "value" after a function has already been 
    // stored in it.
    prot = funcprot();
    funcprot(0);
    value = varargin(ivar);
    funcprot(prot);
    options = optimset_configure (options,key,value);
  end
endfunction
//
// optimset_configure --
//   Configure the given key with the given value.
// Arguments
//   options : a optimset struct
//   key : a string, the key to configure
//   value : the value
//
function options = optimset_configure ( options , key , value )
    select key
    case "Display" then
      if ( ...
      ( value == "off" ) | ...
      ( value == "iter"  ) | ...
      ( value == "final" ) | ...
      ( value == "notify" ) ...
      ) then
        options.Display = value;
      else
        errmsg = msprintf(gettext("%s: Unrecognized value ''%s'' for ''Display'' option."), "optimset", value )
        error(errmsg)
      end
    case "FunValCheck" then
      options.FunValCheck = value;
    case "MaxFunEvals" then
      options.MaxFunEvals = value;
    case "MaxIter" then
      options.MaxIter = value;
    case "OutputFcn" then
      options.OutputFcn = value;
    case "PlotFcns" then
      options.PlotFcns = value;
    case "TolFun" then
      options.TolFun = value;
    case "TolX" then
      options.TolX = value;
    case "Algorithm" then
      options.Algorithm = value;
    case "TolCon" then
      options.TolCon = value;
    case "GradObj" then
      options.GradObj = value;
    case "GradConstr" then
      options.GradConstr = value;
    else
      errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "optimset", key)
      error(errmsg)
    end
endfunction
//
// optimset_new --
//   Returns an empty optimset struct
//
function options = optimset_new ()
    options = struct(...
      "Display"     , [], ...
      "FunValCheck" , [], ...
      "MaxFunEvals" , [], ...
      "MaxIter"     , [], ...
      "OutputFcn"   , [], ...
      "PlotFcns"    , [], ...
      "TolFun"      , [], ...
      "TolX"        , [], ...
      "TolCon"      , [], ...
      "GradObj"     , [], ...
      "Algorithm"     , [], ...
      "GradConstr"  , [] ...
      );
endfunction
//
// optimset_method --
//   Returns an optimset struct which content is the default 
//   settings for the given method.
// Arguments
//   method : a string, the name of the method
//
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

    case "fgoalattain" then
      options = optimset_configure ( options , "Display" , "final" );
      options = optimset_configure ( options , "TolFun" , 1.e-6 );
      options = optimset_configure ( options , "TolCon" , 1.e-6 );
      options = optimset_configure ( options , "GradConstr" , "off" );
      options = optimset_configure ( options , "GradObj" , "off" );

    else
      errmsg = msprintf(gettext("%s: No default options available: the function ''%s'' does not exist on the path."), "optimset", method)
      error(errmsg)
    end
endfunction

