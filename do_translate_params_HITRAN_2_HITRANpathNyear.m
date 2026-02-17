function toptsx = do_translate_params_HITRAN_2_HITRANpathNyear(topts);

%% check to see if earlier driver codes are using topts.HITRAN instead of topts.HITRANpathNyear

toptsx = topts;

if isfield(topts,'HITRAN') & isfield(topts,'HITRANpathNyear')
  fprintf(1,'oh oh topts has conflicting fields HITRANpathNyear and HITRAN')
  topts
  error('please fix by getting rid of topts.HITRAN since I dunno which field you really want')
end

allowedparams = {'HITRAN'};
optvar = fieldnames(topts);
for i = 1 : length(optvar)
 if (length(intersect(allowedparams,optvar{i})) == 1)
   disp('translating field "HITRAN" in topts to "HITRANpathNyear" ')
   toptsx = rmfield(toptsx,'HITRAN');
   toptsx.HITRANpathNyear = topts.HITRAN;
  end
end
