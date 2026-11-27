function d=rp_stage2_distance_complete(c,z)
%RP_STAGE2_DISTANCE_COMPLETE Base bank and width-boundary extensions.
d=rp_stage2_distance(c,z); d=rp_stage2_extend(c,d);
d=rp_stage2_apply_refinements(c,d);
end
