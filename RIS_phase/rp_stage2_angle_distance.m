function result=rp_stage2_angle_distance(c,d)
%RP_STAGE2_ANGLE_DISTANCE Serial convenience wrapper over resumable angle cases.
cases=cell(1,numel(c.study.angleGridDeg));
for j=1:numel(cases), cases{j}=rp_stage2_angle_case(c,d,j); end
result=rp_stage2_merge_angles(c,cases,d.distanceM);
end
