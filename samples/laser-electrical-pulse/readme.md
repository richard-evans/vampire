# Laser and electrical pulse

The `laser-electrical-pulse` program on a single nickel spin at 300 K: an
ultrafast laser pulse heats the electrons (two-temperature model), then a
Gaussian spin-transfer torque pulse arrives 500 fs after the laser peak.
The timing is explained in the comments in `input`; the material and
spin-transfer torque parameters are in `Ni-stt.mat`.

```bash
cd samples/laser-electrical-pulse
../../vampire-serial
```

`output` has columns for time, electron and phonon temperature, the fractional
electric field strength and the magnetisation. Plot the magnetisation and
electron temperature against time to see the demagnetisation and the effect of
the torque pulse.
