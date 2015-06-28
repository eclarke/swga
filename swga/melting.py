'''
melting.py

Calculates the melting temperatures of a given nucleotide sequence.

The melting temperature equations, including correction for salt concentrations,
derives from information found on this page:
  https://www.idtdna.com/Calc/Analyzer/Home/Definitions#MeltTemp

The values for the delta H (dH) and delta S (dS) of each nucleotide pair, and
the initial adjustments to dH and dS are derived from:
  Allawi and SantaLucia (1997), Biochemistry 36 : 10581-10594.

Based heavily on code by the following authors:
- Sebastian Bassi <sbassi@genesdigitales.com>
- Greg Singer <singerg@tcd.ie>
- Nicolas Le Novere <lenov@ebi.ac.uk>
- Calvin Morrison <mutantturkey@gmail.com>
'''
from __future__ import division
from math import sqrt, log
import click
import signal


def _is_sym(s):
    """Returns True if s is symmetric (same as rev. complement)"""
    r = {'A': 'T',
         'T': 'A',
         'G': 'C',
         'C': 'G'}
    return s == ''.join([r[i] for i in s][::-1])


def _overcount(st, p):
    ocu = 0
    x = 0
    while 1:
        try:
            i = st.index(p, x)
        except ValueError:
            break
        ocu += 1
        x = i + 1
    return ocu


def _tercorr(st):
    _dh = 0
    _ds = -1.4 if _is_sym(st) else 0
    start = st[0]
    end = st[-1]

    if start == 'G' or start == 'C':
        _dh += 0.1
        _ds -= 2.8
    elif start == 'A' or start == 'T':
        _dh += 2.3
        _ds += 4.1

    if end == 'G' or end == 'C':
        _dh += 0.1
        _ds -= 2.8
    elif end == 'A' or end == 'T':
        _dh += 2.3
        _ds += 4.1
    return _dh, _ds


def Tm(s, DNA_c=5000.0, Na_c=10.0, Mg_c=20.0, dNTPs_c=10.0, correction=True):
    '''
    Returns the DNA/DNA melting temp using nearest-neighbor thermodynamics.

    This function returns better results than EMBOSS DAN because it uses updated
    thermodynamics values and takes into account initialization parameters from
    the work of SantaLucia (1998).

    Corrects for mono- and divalent cation concentrations.

    Arguments:
    - DNA_c:   DNA concentration [nM]
    - Na_c:    Na+ concentration [mM]
    - Mg_c:    Mg2+ concentration [mM]
    - dNTPs_c: dNTP concentration [mM]
    - correction: correct for cation concentration?
    '''

    R = 1.987    # Universal gas constant (cal/(K*mol))
    s = s.upper()
    dh, ds = _tercorr(s)
    k = DNA_c * 1e-9

    # Adapted from Table 1 in Allawi and SantaLucia (1997).
    # delta H (kcal/mol)
    dh_coeffs = {"AA": -7.9,  "TT": -7.9,
                 "AT": -7.2,
                 "TA": -7.2,
                 "CA": -8.5,  "TG": -8.5,
                 "GT": -8.4,  "AC": -8.4,
                 "CT": -7.8,  "AG": -7.8,
                 "GA": -8.2,  "TC": -8.2,
                 "CG": -10.6,
                 "GC": -9.8,
                 "GG": -8.0,  "CC": -8.0}

    # delta S (eu)
    ds_coeffs = {"AA": -22.2,  "TT": -22.2,
                 "AT": -20.4,
                 "TA": -21.3,
                 "CA": -22.7,  "TG": -22.7,
                 "GT": -22.4,  "AC": -22.4,
                 "CT": -21.0,  "AG": -21.0,
                 "GA": -22.2,  "TC": -22.2,
                 "CG": -27.2,
                 "GC": -24.4,
                 "GG": -19.9,  "CC": -19.9}

    # Multiplies the number of times each nuc pair is in the sequence by the
    # appropriate coefficient, then returns the sum of all the pairs
    dh = dh + \
        sum(_overcount(s, pair) * coeff for pair, coeff in dh_coeffs.items())
    ds = ds + \
        sum(_overcount(s, pair) * coeff for pair, coeff in ds_coeffs.items())

    fgc = len(filter(lambda x: x == 'G' or x == 'C', s)) / float(len(s))

    # Melting temperature
    tm = (1000 * dh) / (ds + (R * log(k)))

    if not correction:
        return tm - 273.15

    MNa = Na_c * 1e-3
    MMg = Mg_c * 1e-3
    MdNTPs = dNTPs_c * 1e-3

    # Free magnesium concentration
    Ka = 3e4  # association constant in biological buffers
    D = (Ka * MdNTPs - Ka * MMg + 1)**2 + (4 * Ka * MMg)
    Fmg = (-(Ka * MdNTPs - Ka * MMg + 1) + sqrt(D)) / (2 * Ka)

    cation_ratio = sqrt(Fmg) / MNa if MNa > 0 else 7.0

    if cation_ratio < 0.22:
        tm = 1 / (
            (1 / tm) +
            ((4.29 * fgc - 3.95) * log(MNa) + 0.94 * log(MNa)**2) * 1e-5)
    else:
        a = 3.92
        d = 1.42
        g = 8.31
        Fmg = MMg
        if cation_ratio < 6.0:
            a = a * (0.843 - 0.352 * sqrt(MNa) * log(MNa))
            d = d * \
                (1.279 - 4.03 * log(MNa) * 1e-3 - 8.03 * log(MNa)**2 * 1e-3)
            g = g * (0.486 - 0.258 * log(MNa) + 5.25 * log(MNa)**3 * 1e-3)
        tm = 1 / (
            (1 / tm) +
            (a - 0.911 * log(Fmg) + fgc * (6.26 + d * log(Fmg)) +
             1 / (2 * (len(s) - 1)) * (-48.2 + 52.5 * log(Fmg) +
                                       g * log(Fmg)**2)) * 1e-5)

    return tm - 273.15


def _init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def _Tm(primer, **kwargs):
    primer.tm = Tm(primer.seq, **kwargs)
    return primer


@click.command()
@click.argument("sequence")
@click.option("-c", "--correction", type=bool, default=True)
@click.option("--dna", default=5000.0)
@click.option("--na", default=10.0)
@click.option("--mg", default=20.0)
@click.option("--dntp", default=10.0)
def main(sequence, dna, na, mg, dntp, correction):
    print Tm(sequence, DNA_c=dna, Na_c=na, Mg_c=mg, dNTPs_c=dntp, correction=correction)


if __name__ == "__main__":
    main()
