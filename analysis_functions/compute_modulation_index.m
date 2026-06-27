function mi = compute_modulation_index(vals_base, vals_drug)
    % MI = (drug - baseline) / (|drug| + |baseline|), bounded [-1, 1].
    valid = ~isnan(vals_base) & ~isnan(vals_drug);
    vb = vals_base(valid); vd = vals_drug(valid);
    denom = abs(vd) + abs(vb);
    ok = denom > 0;
    mi = (vd(ok) - vb(ok)) ./ denom(ok);
end
