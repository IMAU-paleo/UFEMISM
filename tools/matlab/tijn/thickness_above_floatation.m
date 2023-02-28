function TAF = thickness_above_floatation( Hi, Hb, SL)
seawater_density = 1028;
ice_density      = 910;
TAF = Hi - max( 0, (SL - Hb) * (seawater_density / ice_density));
end