function [g,l,a] = mygini(pop,val,makeplot)
    % check arguments
    assert(nargin >= 2, 'gini expects at least two arguments.')
   
    
    assert(numel(pop) == numel(val), ...
        'gini expects two equally long vectors (%d ~= %d).', ...
        size(pop,1),size(val,1))
    pop = [0;pop(:)]; val = [0;val(:)];     % pre-append a zero
    isok = all(~isnan([pop,val]'))';        % filter out NaNs
    if sum(isok) < 2
        warning('gini:lacking_data','not enough data');
        g = NaN; l = NaN(1,4);
        return;
    end
    pop = pop(isok); val = val(isok);
    
    %assert(all(pop>=0) && all(val>=0), ...
     %   'gini expects nonnegative vectors (neg elements in pop = %d, in val = %d).', ...
      %  sum(pop<0),sum(val<0))
    
    % process input
    z = val .* pop;
    [~,ord] = sort(val);
    pop    = pop(ord);     z    = z(ord);
    pop    = cumsum(pop);  z    = cumsum(z);
    relpop = pop/pop(end); relz = z/z(end);
    
    
      g = 1 - sum((relz(1:end-1)+relz(2:end)) .* diff(relpop));
    
      
      l = [relpop,relz];
    a = [pop,z];
    if makeplot
        plot(relpop,relz)
        hold on
        plot([0 1],[0 1],'k')
    end
    
end