// Copyright (C) 2008-2009 - INRIA - Michael Baudin
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function val = optimget (varargin)
  //   Queries an optimization structure.
  //
  // Calling Sequence
  //   val = optimget ( options , key )
  //   val = optimset ( options , key , value )
  //
  // Description
  //   The key is searched in the list of existing fields,
  //   ignoring the case.
  //   Only leading characters of the field name are sufficient to
  //   make the search successful.
  //
  // Authors
  // Michael Baudin, DIGITEO, 2010

  [lhs,rhs]=argn();
  if rhs<=1 | rhs >=4 then
    errmsg = error(msprintf(gettext("%s: Wrong number of arguments : %d expected while %d given"),"optimget", 2,rhs));
    error(errmsg);
  end
  options = varargin(1);
  // Search the field by index
  key = varargin(2);
  fields = [
    "Display"
    "FunValCheck"
    "MaxFunEvals"
    "MaxIter"
    "OutputFcn"
    "PlotFcns"
    "TolFun"
    "TolX"
    "TolCon"
    "Algorithm"
    "GradObj"
    "GradConstr"
    "MinAbsMax"
    ];
  // Search for the given key in the list of available fields.
  // Use a regexp which ignores the case.
  regstr = "/" + key + "/i";
  [r,w] = grep ( fields , regstr , "r" );
  opsize = size(r,2)
  if opsize<>1 then
    matching = strcat(fields(r(1:opsize))," ");
    errmsg = error(msprintf(gettext("%s: Ambiguous property name %s matches several fields : %s"),"optimget", key , matching));
    error(errmsg);
  end
  // Get the matching field
  name = fields(r(1));
  val = options(name);
  // When the value is given and the field is empty, return the value.
  if rhs==3 then
    if val == [] then
      val = varargin(3)
    end
  end
endfunction
