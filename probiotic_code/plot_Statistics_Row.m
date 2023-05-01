function sb_id = plot_Statistics_Row(n)
    b = 4; % base size
    sb_grid = [b b*n];
    
    % pattern of grid (MUST ONLY BE ONE ROW)
    pathp = repmat([0 0 1 1; 0 0 1 1; 0 0 0 0; 0 0 0 0],1,n);
    patv = repmat([0 1 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0],1,n);
    path = repmat([0 0 0 0; 0 0 0 0; 0 0 1 1; 0 0 0 0],1,n);
    
    tot_grid = reshape(1:sb_grid(1)*sb_grid(2),[sb_grid(2),sb_grid(1)])';
    
    hp_area = tot_grid(pathp == 1);
    v_area = tot_grid(patv == 1);
    h_area = tot_grid(path == 1);
    
    m1 = length(hp_area)/n;
    m2 = length(v_area)/n;
    sb_id = cell(n,1);
    for i = 1:n
        sp = i*m1 - (m1-1);
        ep = m1*i;
        hp_in = sort(hp_area(sp:ep))';
    
        sp = i*m2 - (m2-1);
        ep = m2*i;
        v_in = sort(v_area(sp:ep))';
        h_in = sort(h_area(sp:ep))';
    
        sb_id(i) = {{sb_grid,hp_in,v_in,h_in}};
    end

end