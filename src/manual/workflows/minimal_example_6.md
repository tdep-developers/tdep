title: Minimal example VI
author: Olle Hellman

In the previous sections we dealt with calculations at a single temperature. Sometimes you want properties as a function of temperature. In `some folder` I prepared some input for ScF<sub>3</sub> @note Or something else.

#### Phonons dispersions as a function of temperature

Make them plot in all three folders.

Point out that you should only run extractforceconstants once.

Plot all three.

#### Interpolate

For some reason you realise that you where not interested in the phonon dispersions at x,y,or z K, for reasons unknown you want to determine the dispersions at XX K instead. The naive way would be to do a new calculation at that temperature, but that takes a lot of time. It's better to interpolate.

Interpolating the frequencies themselves could be done, but it's not great. Empirically, I have found out that if you interpolate the forceconstants instead, you will have better luck.

This is a little involved, so we will break it into steps.

Provided you followed the instructions and have run extract forceconstants in all three directories, you might have noted that there are files called `outfile.phi_secondorder`. As [described](../../program/extract_forceconstants.html), these are a list of the irreducible forceconstants. Copy this to an infile

```
cp outfile.phi_secondorder infile.phi_secondorder
```
And run

```
./extract_forceconstants -r
```
Note that this runs much faster: by invoking the `-r` flag, it will read the values from file. Now, edit `infile.phi_secondorder` and change some random value to 100, run

```
./extract_forceconstants -r
phonon_dispersion_relations
gnuplot outfile.dispersion_relations.gnuplot
```
The dispersions will probably look really strange. The point of this exercise was to understand how to modify an irreducible component and get a new `outfile.forceconstant` out.

Since we know how to translate irreducible components to forceconstants, all we need to to is interpolate each irreducible component separately, translate to forceconstants, and calculate dispersions. This way you can get the dispersions at any temperature in the interval.

I like to use matlab for this step, but it should be just as easy in python or some other scripting language. I will not go into details, since I assume you know how to interpolate things. You should construct a script that

1.	Reads all N `outfile.phi_secondorder`
2.	For each of the values in `outfile.phi_secondorder`, create an interpolation that can evaluate that value at any temperature. A second order polynomial works great.
3.	At a temperature of interest, print a new `infile.phi_secondorder` with the interpolated values and run `extract_forceconstants -r`. This gives you a forceconstant at this temperature.
4.	Run `phonon_dispersion_relations` with the new forceconstant, and store `outfile.dispersion_relations` somewhere.
5.	Repeat this for say 50 temperatures in the valid range. Plot the results somehow.

If successful, it should look something like this:

@todo Some picture

Or alternatively, just tracking the TO mode at gamma:

@todo Some other picture.


