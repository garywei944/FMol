import textwrap

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

from mcts import parameters as p


def wavelength_to_rgb(wavelength, gamma=0.8):
    """
    taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range

    :param wavelength: wavelength
    :param gamma: gamma
    :return: color (r, g, b, a)
    """
    wavelength = float(wavelength)
    if 380 <= wavelength <= 750:
        a = 1.
    else:
        a = 0.5
    if wavelength < 380:
        wavelength = 380.
    if wavelength > 750:
        wavelength = 750.
    if 380 <= wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        r = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        g = 0.0
        b = (1.0 * attenuation) ** gamma
    elif 440 <= wavelength <= 490:
        r = 0.0
        g = ((wavelength - 440) / (490 - 440)) ** gamma
        b = 1.0
    elif 490 <= wavelength <= 510:
        r = 0.0
        g = 1.0
        b = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif 510 <= wavelength <= 580:
        r = ((wavelength - 510) / (580 - 510)) ** gamma
        g = 1.0
        b = 0.0
    elif 580 <= wavelength <= 645:
        r = 1.0
        g = (-(wavelength - 645) / (645 - 580)) ** gamma
        b = 0.0
    elif 645 <= wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        r = (1.0 * attenuation) ** gamma
        g = 0.0
        b = 0.0
    else:
        r = 0.0
        g = 0.0
        b = 0.0
    return r, g, b, a


def plot_wl(data, smiles):
    """
    Generate a plot of the SMILES with the wavelength and oscillator strength peaks

    :param data: data dict
    :type data: dict
    :param smiles: the SMILES to plot
    :type smiles: str
    :return: None
    """
    if data[smiles]["valid"]:
        clim = (350, 780)
        norm = plt.Normalize(*clim)
        wl = np.arange(clim[0], clim[1] + 1, 2)
        colorlist = list(zip(norm(wl), [wavelength_to_rgb(w) for w in wl]))
        spectralmap = matplotlib.colors.LinearSegmentedColormap.from_list("spectrum", colorlist)

        plt.subplots(1, 1, figsize=(8, 4), tight_layout=True)
        wl = []
        f = []
        for l in data[smiles]['dft']:
            wl.append(l['nm'])
            f.append(l['f'])
        wl = wl[::-1]
        f = f[::-1]
        wavelengths = np.linspace(200, 1000, 1000)
        out = []
        for w, f_ in zip(wl, f):
            if w > 1000:
                out.append((w, f_))
            else:
                c = wavelength_to_rgb(w)
                plt.bar(w, f_, 3, color=c)
        if out:
            text = ""
            for w, f_ in out:
                text += " +" + str(w) + " nm, " + str(f_) + "\n"
            plt.text(800, 0.55, text)
        y = np.linspace(-0.05, 0.6, 1000)
        x_, y_ = np.meshgrid(wavelengths, y)

        extent = (np.min(wavelengths), np.max(wavelengths), np.min(y), np.max(y))

        plt.imshow(x_, clim=clim, extent=extent, cmap=spectralmap, aspect='auto')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Oscillator Strength')
        s = "".join(p.config['long_prefix']) + smiles
        s = '{:4d}'.format(data[smiles]["id"]) + " : " + s
        s = textwrap.fill(s, 50)
        plt.title(s)
        plt.fill_between(wavelengths, 0.6, color='w')

        plt.savefig("../data_out/" + p.config["configuration_name"] + "/plot/" + str(data[smiles]["id"]) + '_wl.png', dpi=200)

        plt.show()
