function result = vap_run_cases(P, T_L, H, Qdot, p)
%VAP_RUN_CASES Evaluate an array while preserving every status and residual.
% H is either the same shape as Qdot or one height per row.
if ~isequal(size(H),size(Qdot))
    if isvector(H) && numel(H)==size(Qdot,1)
        H=repmat(H(:),1,size(Qdot,2));
    else
        error('vap:Shape','H must match Qdot or contain one height per row.');
    end
end
sz=size(Qdot);
for name={'mdot','T_out','x_out','pressure_residual','energy_residual'}
    result.(name{1})=nan(sz);
end
result.ok=false(sz);
result.regime=cell(sz); result.status=cell(sz); result.message=cell(sz);
for k=1:numel(Qdot)
    r=vaporizer_thermosiphon(P,T_L,H(k),Qdot(k),p);
    result.ok(k)=r.ok; result.regime{k}=r.regime;
    result.status{k}=r.status; result.message{k}=r.message;
    if r.ok
        for name={'mdot','T_out','x_out','pressure_residual','energy_residual'}
            result.(name{1})(k)=r.(name{1});
        end
    end
end
result.failed_count=nnz(~result.ok);
if result.failed_count>0
    warning('vap:FailedCases','%d/%d cases failed; values remain NaN. Inspect status and message.',result.failed_count,numel(Qdot));
end
end
