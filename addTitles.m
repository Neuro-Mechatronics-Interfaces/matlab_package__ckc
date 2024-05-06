function addTitles(L, Title, Subtitle, TitleColor, SubtitleColor, options)
%ADDTITLES  Adds title and subtitle text (optionally) to graphic handle and specifies the color/font.
%
% Syntax:
%   ckc.addTitles(L, Title, Subtitle, TitleColor, SubtitleColor);
%   [hTitle, hSubtitle] = ckc.addTitles(...);
%   ckc.addTitles(...,'FontName', value);
%
% Inputs: 
%     L - Graphics handle to target
%     Title {mustBeTextScalar} = "";
%     Subtitle {mustBeTextScalar} = "";
%     TitleColor (1,3) {mustBeInRange(TitleColor,0,1)} = [0 0 0];
%     SubtitleColor (1,3) {mustBeInRange(SubtitleColor,0,1)} = [0.65 0.65 0.65];
%
% Options:
%     options.FontName {mustBeTextScalar} = 'Tahoma';

arguments
    L
    Title {mustBeTextScalar} = "";
    Subtitle {mustBeTextScalar} = "";
    TitleColor (1,3) {mustBeInRange(TitleColor,0,1)} = [0 0 0];
    SubtitleColor (1,3) {mustBeInRange(SubtitleColor,0,1)} = [0.65 0.65 0.65];
    options.FontName {mustBeTextScalar} = 'Tahoma';
end
if strlength(Title) > 0
    title(L, strrep(Title,"_","\_"), 'FontName', options.FontName, 'Color', TitleColor);
end

if strlength(Subtitle) > 0
    subtitle(L, strrep(Subtitle,"_","\_"), 'FontName', options.FontName, 'Color', SubtitleColor);
end

end