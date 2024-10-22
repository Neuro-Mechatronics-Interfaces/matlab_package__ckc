function snips = extended_2_snips(Ye, pulseSamples, extFactor)
%EXTENDED_2_SNIPS  Converts extended data array into snippets tensor based on sample indices.
r = size(Ye,1);
n = numel(pulseSamples);
snips = reshape(Ye(:,pulseSamples),extFactor,r/extFactor,n);
end