function [Fac_est, Xrec, amplitudes] = solve_parafac(method,X,d,options)

% SOLVE_PARAFAC   Compute a solution to the R-D PARAFAC model
%
% Syntax:
%    FACTORS = SOLVE_PARAFAC(METHOD,X,d[,options])
%    [FACTORS,Xrec] = SOLVE_PARAFAC(METHOD,X,d[,options])
%
% Input:
%    METHOD  - Which method shall be used to solve the PARAFAC problem?
%       'Jointdiag' - reduce the problem onto joint diagonalization of matrices
%       'MALS'      - multilinear alternating least squares ("plain vanilla")
%       'PARAFAC'   - PARAFAC algorithm from the N-way toolbox
%       'PARAFACnn' - PARAFAC algorithm from the N-way toolbox, with
%                     non-negativity constraint
%       'COMFAC'    - COMFAC (Sidiropoulos, Bro)
%       'DTLD'      - direct trilinear decomposition (Sanchez Kowalski)
%       'GRAM'      - generalized rank annihilation method
%    X       - (noisy) measurement tensor to decompose
%    d       - number of components to estimate
%    options - struct specifying additional options. May contain:
%         options.speedup - integer specifying speedup options. Default is
%               0 (no speedup), higher values indicate more speedups.
%         options.donormalize - if set, all factors are normalized such
%               that their first row is equal to one. Defaults to one (1).
%               Note that currently, this option is ignored for
%               'JointDiag'.
%         options.jdoptions - a struct that is passed to
%               solve_parafac_jd_Rd, possibly containing more options for
%               this function.
%
% Output:
%    Factors - length-R cell array containing the estimates for each of the
%         R factors.
%    Xrec    - optional second output argument containing the reconstructed
%         tensor (to compute reconstruction errors).
%
% Author:
%    Florian Roemer, Communications Resarch Lab, TU Ilmenau
% Date:
%    Dec 2007

DEFAULT_speedup = 0;
DEFAULT_donormalize = 0;
if nargin < 4
    options = struct('speedup',DEFAULT_speedup);
end

if ~isfield(options,'speedup')
    options.speedup = DEFAULT_speedup;
end
if ~isfield(options,'donormalize')
    options.donormalize = DEFAULT_donormalize;
end

if strcmp(lower(method),'jdobc bp')
    method = 'jdobc';
    options.speedup = 100;
end

switch lower(method)
    case 'jdobc'
        [S,U] = hosvd(X);
        [Sc,Uc] = cuthosvd(S,U,d);
        if options.speedup > 1
            [Fac_est,Residuals] = solve_parafac_jd_RD_onlybc(Sc,Uc,d);%,options.jdoptions);
            Fac_est = eliminatefactors(Fac_est,Residuals,options.speedup);
        else
            Fac_est = solve_parafac_jd_RD_onlybc(Sc,Uc,d);%,options.jdoptions);
        end
        %Fac_est0 = Fac_est;
        %Fac_est = optassignment_Rd(Fac_est);
        if nargout > 1
            [which,Xrec] = bestmatch_Rd_new(X, Fac_est);
        else
            which = bestmatch_Rd_new(X, Fac_est);
        end            
        for r = 1:length(Fac_est)
            Fac_est{r} = Fac_est{r}{which(r)};
        end
%         if options.donormalize ~= 2
%             [amplitudes,Fac_est] = estimate_amplitudes(X, normalizefactors(Fac_est));
%         end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    case 'jdobc0'
        [S,U] = hosvd(X);
        [Sc,Uc] = cuthosvd(S,U,d);
        options.jdoptions.condtolerance = 0;
        if options.speedup > 1
            [Fac_est,Residuals] = solve_parafac_jd_RD_onlybc(Sc,Uc,d,options.jdoptions);
            Fac_est = eliminatefactors(Fac_est,Residuals,options.speedup);
        else
            Fac_est = solve_parafac_jd_RD_onlybc(Sc,Uc,d,options.jdoptions);
        end
        %Fac_est0 = Fac_est;
        %Fac_est = optassignment_Rd(Fac_est);
        if nargout > 1
            [which,Xrec] = bestmatch_Rd_new(X, Fac_est);
        else
            which = bestmatch_Rd_new(X, Fac_est);
        end            
        for r = 1:length(Fac_est)
            Fac_est{r} = Fac_est{r}{which(r)};
        end
%         if options.donormalize ~= 2
%             [amplitudes,Fac_est] = estimate_amplitudes(X, normalizefactors(Fac_est));
%         end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end        
    case 'jdoc'
        [S,U] = hosvd(X);
        [Sc,Uc] = cuthosvd(S,U,d);
        if options.speedup > 1
            [Fac_est,Residuals] = solve_parafac_jd_RD_conly(Sc,Uc,d,options.jdoptions);
            Fac_est = eliminatefactors(Fac_est,Residuals,options.speedup);
        else
            Fac_est = solve_parafac_jd_RD_conly(Sc,Uc,d);%,options.jdoptions);
        end
        %Fac_est0 = Fac_est;
        %Fac_est = optassignment_Rd(Fac_est);
        if nargout > 1
            [which,Xrec] = bestmatch_Rd_new(X, Fac_est);
        else
            which = bestmatch_Rd_new(X, Fac_est);
        end            
        for r = 1:length(Fac_est)
            Fac_est{r} = Fac_est{r}{which(r)};
        end
%         if options.donormalize ~= 2
%             [amplitudes,Fac_est] = estimate_amplitudes(X, normalizefactors(Fac_est));
%         end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    case 'jointdiag_normrv'
        if ~isfield(options,'jdoptions')
            options.jdoptions = struct;
        end
        if options.speedup > 0
            options.jdoptions.compute_combined = 0;
        end
        
        [S,U] = hosvd(X);
        [Sc,Uc] = cuthosvd(S,U,d);
        Fac_est = solve_parafac_jd_RD_normrv(Sc,Uc,d,options.jdoptions);
        %Fac_est0 = Fac_est;
        Fac_est = optassignment_Rd_sign(Fac_est);
        if nargout > 1
            [which,Xrec] = bestmatch_Rd_new(X, Fac_est);
        else
            which = bestmatch_Rd_new(X, Fac_est);
        end            
        for r = 1:length(Fac_est)
            Fac_est{r} = Fac_est{r}{which(r)};
        end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
        
    case 'jointdiag'
        if ~isfield(options,'jdoptions')
            options.jdoptions = struct;
        end
        if options.speedup > 0
            options.jdoptions.compute_combined = 0;
        end
%         if ~isfield(options.jdoptions,'donormalize')
%             options.jdoptions.donormalize = options.donormalize;
%         end
      
        
        [S,U] = hosvd(X);
        [Sc,Uc] = cuthosvd(S,U,d);
        if options.speedup > 1
            [Fac_est,Residuals, thecond_all] = solve_parafac_jd_RD(Sc,Uc,d,options.jdoptions);
            Fac_est = eliminatefactors(Fac_est,Residuals,options.speedup);
        else
           [Fac_est,Residuals, thecond_all] = solve_parafac_jd_RD(Sc,Uc,d,options.jdoptions);
        end
        save fac_est_closed_form_flor;        
        %Fac_est0 = Fac_est;
        %Fac_est = optassignment_Rd(Fac_est);
        if nargout > 1
            [which,Xrec] = bestmatch_Rd_new(X, Fac_est);
        else
            which = bestmatch_Rd_new(X, Fac_est);
        end            
        for r = 1:length(Fac_est)
            Fac_est{r} = Fac_est{r}{which(r)};
        end
        if options.donormalize ~= 2
            [amplitudes,Fac_est] = estimate_amplitudes(X, normalizefactors(Fac_est));
        end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    case 'mals'
        maxit = 10000;
        if options.speedup>0
            maxit = ceil(maxit / (options.speedup+1));
        end
        Fac_est = plainvanilla_tals_Rd(X,d,[],maxit);
        if options.donormalize
            Fac_est = normalizefactors(Fac_est);        
        end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    case 'parafac'
        if ~isreal(X)
            error('PARAFAC method only applicable to real-valued problems.');
        end
    	Fac_est = parafac(X,d,[0,0,0,0,NaN]);
        if options.donormalize
            Fac_est = normalizefactors(Fac_est);
        end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    case 'parafacnn'
        if ~isreal(X)
            error('PARAFAC method only applicable to real-valued problems.');
        end
    	Fac_est = parafac(X,d,[0,0,0,0,NaN],2*ones(1,ndims(X)));
        if options.donormalize
            Fac_est = normalizefactors(Fac_est);
        end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    case 'comfac'
        if ndims(X) > 3
            error('COMFAC method only applicable to 3-D problems.');
        end
        [A, B, C] = comfac(X,d);
        Fac_est = {A,B,C};
        if options.donormalize
            Fac_est = normalizefactors(Fac_est);
        end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    case 'dtld'
        if ndims(X) > 3
            error('DTLD method only applicable to 3-D problems.');
        end
        M = size(X);
        if isreal(X)
            [A, B, C] = dtld(X,d,1);
        else
            [A, B, C] = cdtld(reshape(X,M(1),prod(M(2:3))),M,d,1);
        end
        Fac_est = {A,B,C};
        if options.donormalize
            Fac_est = normalizefactors(Fac_est);
        end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    case 'gram'
        if ndims(X) > 3
            error('GRAM method only applicable to 3-D problems.');
        end
        if isreal(X)
            [Aest,Best,Cest] = gram(X(:,:,1),X(:,:,2),d);
            Cest = unfolding(X,3)*pinv(krp(Aest,Best).');
        else
            [Aest,Best,Cest] = cgram(X(:,:,1),X(:,:,2),d);
            Cest = unfolding(X,3)*pinv(krp(Aest,Best).');
        end
        Fac_est = {Aest,Best,Cest};
        if options.donormalize
            Fac_est = normalizefactors(Fac_est);
        end
        if nargout > 1
            Xrec = build_ten_Rd(Fac_est);
        end
    otherwise
        error('Unkown method.');
end


