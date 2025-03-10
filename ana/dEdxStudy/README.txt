Hi Maria,

you'll find here some "reduced files" with a first selection on the quality of the gammas, and to have a pi0 reconstructed + I've removed some heavy branches to make the processing faster
/exp/minerva/data/users/noeroy/REDUCED_MAD_ELECTRONSCALE/
for instance that file would have enough stat to start with /exp/minerva/data/users/noeroy/REDUCED_MAD_ELECTRONSCALE/MC_1N/tmp_1N.root

Then,  for the gammas, here are the branches you might want to have a look at:

gamma1_E = the calibrated energy of the gamma, however it uses the wrong splines, it's close enough to make first tests but for precision later, you'll want to use that formula for the energy:
(1.326 * (gamma1_evis_trkr + 2.205 * gamma1_evis_ecal + (4 * 2.205 - 1) * gamma1_evis_scal_UV + (2 * 2.205 - 1) *gamma1_evis_scal_X))
which applies the calibration splines for the MC to the original visible energy

gamma1_dEdx which is the computed dEdx, i'll have to check why it's the default -9999 sometimes even when you have a gamma reconstructed. 

Angles (wrt to beam) are gamma1_phi, gamma1_theta

The truth variable would be 
truth_gamma1_E
truth_gamma1_phi
truth_gamma1_theta

Same for gamma2.


The cuts I like to yuse are gamma1_E>0 && gamma2_E>0 (to be extra certain),  and to have good gammas, you can also try 
(gamma1_E+gamma2_E > 400. || pi0_openingAngle > 20) That one makes sure that you don't have low energy colinear gammas, otherwise you  most likely badly reconstruct the pi0s.
