function t = rd_transport(c,photons)
%RD_TRANSPORT Copy of c with a given Monte Carlo packet count.
t = c;
t.mc.photons = round(photons);
t.mc.chunkSize = max(1,min(round(c.dynamic.chunkSize),t.mc.photons));
end
