function entropy = entropyMeasurements(A)

    overlapMatrix = abs(A*A').^2;
    ovs = reshape(overlapMatrix,[numel(overlapMatrix) 1]);
    
    entovs = -sum(ovs.*log(ovs+eps));
    
    entropy = entovs;
    
    
    