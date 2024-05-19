function Y = sig_2_samples(data, big_grid_dims, small_grid_dims)
%SIG_2_SAMPLES Back-convert from DEMUSE "big grid" SIG field to `samples`-like rows-are-channels data array.
arguments
    data struct
    big_grid_dims (1,2) {mustBePositive, mustBeInteger} = [16 16]
    small_grid_dims (1,2) {mustBePositive, mustBeInteger} = [8 8]
end
Y = zeros(numel(data.SIG),numel(data.SIG{1,1}));
[ROWS,COLS] = ckc.grid_of_grids_to_columns(big_grid_dims,small_grid_dims);
for ii = 1:numel(data.SIG)
    if ~isempty(data.SIG{ROWS(ii),COLS(ii)})
        Y(ii,:) = data.SIG{ROWS(ii),COLS(ii)};
    end
end
end