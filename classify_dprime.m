function label = classify_dprime(d, thresh)
    if isnan(d)
        label = "NaN";
    elseif abs(d) < thresh
        label = "nonresponsive";
    elseif d > 0
        label = "enhanced";
    else
        label = "suppressed";
    end
end