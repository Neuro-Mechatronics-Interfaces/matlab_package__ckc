function sortVector = sortByRecruitmentOrder(MUPulses, options)
%SORTBYRECRUITMENTORDER  Returns indexing vector to sort by recruitment
%
% Syntax:
%   sortVector = ckc.sortByRecruitmentOrder(MUPulses, 'Name', value, ...);
%
% Inputs:
%   MUPulses (1,:) cell -- Array of pulse instants
%
% Options:
%   'SortOrder' {mustBeMember(options.SortOrder, {'ascend','descend'})} = 'ascend';
%
% Output:
%   sortVector - 1 x nMUAPs vector indicating ordering to sort by recruitment
%
% See also: Contents, ckc.plotIDR

arguments
    MUPulses (1,:) cell
    options.SortOrder {mustBeMember(options.SortOrder, {'ascend','descend'})} = 'ascend';
end

r = numel(MUPulses);
firstOccurrence = nan(1,r);
for ii = 1:r
    if numel(MUPulses{ii}) > 0
        firstOccurrence(ii) = min(MUPulses{ii});
    else
        firstOccurrence(ii) = inf;
    end
end
[~,sortVector] = sort(firstOccurrence, options.SortOrder);

end