function C = QAMconst(M)
%QAMconst -  QAM Constellation
%
%   Returns a matrix with a 'regular'
%   QAM constellation
%
%   Usage
%      C = QAMconst(M) : M = (2^L) > 1
%

    if ( M <= 1 || mod(log(M),log(2)) ~= 0)
        warning('input must be a power of two');
    else
        p = 2^ceil(log(M)/(2*log(2)));
        q = 2^floor(log(M)/(2*log(2)));
        C = repmat((-p+1:2:p-1),q,1)+1j*repmat((q-1:-2:-q+1)',1,p);
    end

end
