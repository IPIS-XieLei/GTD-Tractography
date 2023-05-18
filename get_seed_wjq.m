function curr_seeds = get_seed_wjq(seedmask)

curr_seeds = zeros(1,3);
 points_seeds_curr = 1;
for in_curr = 1: size(seedmask,1)
     %  index = round(size(seedmask,1)*rand());
        if in_curr<=size(seedmask,1) && in_curr~=0
            
            for  seednum = 1:1
               seed = seedmask(in_curr,:)+[rand(),rand(),rand()]-[0.5,0.5,0.5];
               curr_seeds(points_seeds_curr,:) =  seed;
               points_seeds_curr = points_seeds_curr + 1;
            end            
            
        end
end

 

end