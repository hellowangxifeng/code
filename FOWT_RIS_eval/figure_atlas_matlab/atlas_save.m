function atlas_save(fig,name)
%ATLAS_SAVE Save a 300-dpi PNG and an editable FIG without embedded captions.
if ~isfolder('figures'), mkdir('figures'); end
exportgraphics(fig,fullfile('figures',[name '.png']),'Resolution',300,'BackgroundColor','white');
savefig(fig,fullfile('figures',[name '.fig']));
close(fig);
end
