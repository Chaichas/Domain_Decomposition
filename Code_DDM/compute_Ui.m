function [Uisub] = compute_Ui(Ksub, Fsub, Ub, blockAP, Ns)
    blockUb = blockAP'*Ub;
    Uisub = {};
    last = 0;
    for s = 1:Ns
        if (s == 1) || (s == Ns)
            Nb = 1;
        else
            Nb = 2;
        end
        Uisub{s} = inv(Ksub{s}(1:end-Nb,1:end-Nb))...
            *(Fsub{s}(1:end-Nb,1)-...
              Ksub{s}(1:end-Nb,end-Nb+1:end)*blockUb(last+1:last+Nb,1));
    last = last + Nb;
    end
end