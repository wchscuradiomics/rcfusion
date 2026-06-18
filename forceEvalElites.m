function [state, options, changed] = forceEvalElites(options, state, flag)
changed = false;
switch flag
  case {'init', 'iter'}    
    state.EvalElites = true;
    changed = true; 
end
end