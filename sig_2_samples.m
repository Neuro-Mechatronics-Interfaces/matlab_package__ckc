function Y = sig_2_samples(data, big_grid_dims, small_grid_dims)
%SIG_2_SAMPLES Back-convert from DEMUSE "big grid" SIG field to `samples`-like rows-are-channels data array.
arguments
    data
    big_grid_dims (1,2) = [nan nan];
    small_grid_dims (1,2) = [nan nan];
end
if isstruct(data)
    data = data.SIG;
end
Y = zeros(numel(data),numel(data{1,1}));
if any(isnan(big_grid_dims))
    big_grid_dims = size(data);
end
if any(isnan(small_grid_dims))
    ROWS = repmat((1:big_grid_dims(1))', big_grid_dims(2), 1);
    COLS = repelem(1:big_grid_dims(2), 1, big_grid_dims(1))';
else
    [ROWS,COLS] = ckc.grid_of_grids_to_columns(big_grid_dims,small_grid_dims);
end
for ii = 1:numel(data)
    if ~isempty(data{ROWS(ii),COLS(ii)})
        Y(ii,:) = data{ROWS(ii),COLS(ii)};
    end
end
end