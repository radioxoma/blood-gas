#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Arterial blood gas interpreter.

Eugene Dvoretsky, 2015-05-20


Check abg
    if metabolic alcalosis:
        check anion gap
Check A-a gradient? [< 2.6 kPa (20 mmHg)]

[The Computing Techniques](http://www.acid-base.com/computing.php)

$ md5sum *
41f8e7d98fcc26ea2319ac8da72ed8cd ABL800 Reference Manual English US.pdf
a00e5f337bd2c65e513fda1202827c6a ABL800 Operators Manual English US.pdf

I do not perform any algebraic optimization.
"""

from __future__ import absolute_import
from __future__ import division
try:
    from uncertainties import umath as math
except ImportError:
    import math


# Arterial blood reference
norm_pH = (7.35, 7.45)
norm_pCO2 = (35., 45.)  # 4.7-6.0 kPa
norm_HCO3 = (22., 26.)
norm_pO2 = (80., 100.)
live_pH = (6.8, 7.8)  # Live borders


kPa = 0.133322368  # kPa to mmHg, 1 mmHg = 0.133322368 kPa


def calculate_anion_gap(Na, Cl, HCO3act, K=0.0, albuminum=None):
    """Calculate serum 'Anion Gap' or 'Anion Gap (K+)'.

    Помогает при выявлении противонаправленных метаболических процессов.
    Напрмиер потеря Cl (алкалоз) и лактат-ацидоз.

    May be known as SID [1], AG.

    Corresponds to phosphates, sulphates, proteins.
    High gap: acute kidney injury, lactate, ketoacidosis, salicylate ->
        secondary loss of HCO3− which is a buffer, without a concurrent
        increase in Cl− for electroneutrality equilibrium support.
    Low gap: increase in Cl−.


    Examples
    --------

    To reproduce 'Radiometer ABL800 Flex' Anion Gap calculation:

    >>> abg.calculate_anion_gap(
        Na=173, Cl=77, HCO3act=abg.calculate_hco3p(pH=6.656, pCO2=27.9))
    93.0681487508615


    References
    ----------

    [1] Kostuchenko S.S., ABB in the ICU, 2009, p. 59
    [2] Patrick J Neligan MA MB FCARCSI, Clifford S Deutschman MS MD FCCM
        Acid base balance in critical care medicine
    [3] https://en.wikipedia.org/wiki/Anion_gap

    :param float Na: Serum sodium, mmol/L.
    :param float Cl: Serum chloride, mmol/L.
    :param float HCO3act: Serum actual bicarbonate (HCO3(P)), mmol/L.
    :param float K: Serum potassium, mmol/L.
        If not given returns AG, otherwise AG(K). Serum potassium value
        usually low and frequently omitted. Normal reference values different
        for potassium and non-potassium anion gap.
    :param float albuminum: Protein correction, g/dL. If not given,
        hypoalbuminemia leads to lower anion gap.
    :return:
        Anion gap mEq/L.
    :rtype: float
    """
    anion_gap = (Na + K) - (Cl + HCO3act)
    if albuminum is not None:
        # Protein correction. Normal albuminum 2.5 g/dL see [1] p. 61.
        anion_gap += 2.5 * (4.4 - albuminum)
    return anion_gap


def calculate_mosm(Na, glucosae):
    """Calculate serum osmolarity (mOsm).


    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-41, p. 277, equation 48.

    :param float Na: mmol/L
    :param float glucosae: mmol/L
    :return:
        Serum osmolarity, mmol/kg.
    :rtype: float
    """
    return 2 * Na + glucosae


def calculate_hco3(pH, pCO2):
    """Concentration of HCO3 in plasma (actual bicarbonate).

    May be known as HCO3act, cHCO3(P)).

    Generic approximation of Henderson-Hasselbalch equation.
    Good for water solutions, lower precision in lipemia cases.


    References
    ----------

    [1] http://www-users.med.cornell.edu/~spon/picu/calc/basecalc.htm

    :param float pH:
    :param float pCO2: mmHg
    :return:
        HCO3act mmol/L.
    :rtype: float
    """
    # 0.03 - CO2 solubility coefficient mmol/L/hg
    # 6.1 - dissociation constant for H2CO3
    return 0.03 * pCO2 * 10 ** (pH - 6.1)


def calculate_hco3p(pH, pCO2):
    """Concentration of HCO3 in plasma (actual bicarbonate).

    May be known as HCO3act, cHCO3(P)) [1].

    Sophisticated calculation by [1]. Good for water solutions,
    lower precision in lipemia cases.


    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-28, p. 264, equation 4.
    [2] Siggaard-Andersen O, Wimberley PD, Fogh-Andersen N, Gøthgen IH.
        Measured and derived quantities with modern pH and blood gas equipment:
        calculation algorithms with 54 equations. Scand J Clin Lab Invest 1988;
        48, Suppl 189: 7-15. p. 11, equations 6, 7.

    :param float pH:
    :param float pCO2: kPa
    :return:
        cHCO3(P) mmol/L.
    :rtype: float
    """
    pKp = 6.125 - math.log10(1 + 10 ** (pH - 8.7))  # Dissociation constant
    # 1 mmHg = 0.133322368 kPa
    # aCO2(P) solubility coefficient for 37 °С = 0.230 mmol/L/kPa
    # pCO2 *= 0.133322368  # mmHg to kPa
    return 0.230 * pCO2 * 10 ** (pH - pKp)


def calculate_hco3pst(pH, pCO2, ctHb, sO2):
    """Standard Bicarbonate, the concentration of HCO3- in the plasma
    from blood which is equilibrated with a gas mixture with
    pCO2 = 5.33 kPa (40 mmHg) and
    pO2 >= 13.33 kPa (100 mmHg) at 37 °С.

    May be known as cHCO3(P,st)) [1].

    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-29, p. 265, equation 9.

    :param float pH:
    :param float pCO2: kPa
    :param float ctHb: Concentration of total hemoglobin in blood, mmol/L
    :param float sO2: Fraction of saturated hemoglobin, fraction.
    :return:
        cHCO3(P,st) mmol/L.
    :rtype: float
    """
    a = 4.04 * 10 ** -3 + 4.25 * 10 ** -4 * ctHb
    Z = calculate_cbase(pH, pCO2, ctHb=ctHb) - 0.3062 * ctHb * (1 - sO2)
    return 24.47 + 0.919 * Z + Z * a * (Z - 8)


def calculate_be(pH, pCO2, HCO3act):
    """Calculate base excess (BE), Siggaard Andersen approximation.

    Synonym for cBase(Ecf) and SBE?


    References
    ----------

    [1] http://www-users.med.cornell.edu/~spon/picu/calc/basecalc.htm
    [2] http://www.acid-base.com/computing.php

    :param float pH:
    :param float pCO2: mmHg
    :return:
        Base excess, mEq/L.
    :rtype: float
    """
    return 0.9287 * HCO3act + 13.77 * pH - 124.58


def calculate_cbase(pH, pCO2, ctHb=3):
    """Calculate base excess.

    Calculate standard base excess, known as SBE, cBase(Ecf) or
    actual base excess known as ABE, cBase(B).


    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-29, p. 265, equation 5.

    [2] Siggaard-Andersen O. The acid-base status of the blood.
        4th revised ed. Copenhagen: Munksgaard, 1976.
        http://www.cabdirect.org/abstracts/19751427824.html

    :param float pH:
    :param float pCO2: kPa
    :param float ctHb: Concentration of total hemoglobin in blood, mmol/L
        If not given, calculate cBase(Ecf), otherwise cBase(B).
    :return:
        Standard base excess (SBE) or actual base excess (ABE), mEq/L.
    :rtype: float
    """
    # pCO2 *= 0.133322368  # mmHg to kPa
    a = 4.04 * 10 ** -3 + 4.25 * 10 ** -4 * ctHb
    pHHb = 4.06 * 10 ** -2 * ctHb + 5.98 - 1.92 * 10 ** (-0.16169 * ctHb)
    log_pCO2Hb = -1.7674 * (10 ** -2) * ctHb + 3.4046 + 2.12 * 10 ** (
        -0.15158 * ctHb)
    pHst = pH + math.log10(5.33 / pCO2) * (
        (pHHb - pH) / (log_pCO2Hb - math.log10(7.5006 * pCO2)))
    cHCO3_533 = 0.23 * 5.33 * 10 ** ((pHst - 6.161) / 0.9524)
    # There is no comments, as there is no place for weak man
    cBase = 0.5 * ((8 * a - 0.919) / a) + 0.5 * math.sqrt(
        (((0.919 - 8 * a) / a) ** 2) - 4 * ((24.47 - cHCO3_533) / a))
    return cBase


def calculate_hct(ctHb):
    """Calculate hematocrit.

    The ratio between the volume of erythrocytes and the volume of whole blood.
    See [4] for interpretating caveats.


    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-30, p. 266, equation 13.
    [2] http://www.derangedphysiology.com/php/Arterial-blood-gases/
        Haematocrit-is-a-derived-measurment-in-the-blood-gas-analyser.php
    [3] http://www.ncbi.nlm.nih.gov/pubmed/2128562
    [4] http://www.derangedphysiology.com/php/Arterial-blood-gases/
        Haematocrit-is-a-derived-measurment-in-the-blood-gas-analyser.php

    :param float ctHb: Concentration of total hemoglobin in blood, mmol/L.
    :return:
        Hematocrit, fraction (not %).
    :rtype: float
    """
    # ctHb(mmol/L) == ctHb(g/dL) / 1.61140 == ctHb(g/dL) * 0.62058
    # By [1] p. 6-14 or 6-49.
    return 0.0485 * ctHb + 8.3 * 10 ** -3


def calculate_pHT(pH, t):
    """pH of blood at patient temperature.


    Examples
    --------

    >>> calculate_pHT(6.919, 39.6)
    6.8891689
    >>> calculate_pHT(7.509, 38.6)
    7.4845064


    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-28, p. 264, equation 1.

    :param float pH:
    :param float t: Body temperature, °C.
    :return:
        pH of blood at given temperature.
    :rtype: float
    """
    # 37 °C is temperature of measurement in device
    return pH - (0.0146 + 0.0065 * (pH - 7.40)) * (t - 37)


def calculate_pCO2T(pCO2, t):
    """Partial pressure of CO2 in blood at patient temperature.


    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-28, p. 264, equation 3.

    :param float pCO2: kPa or mmHg (sic!)
    :param float t: Body temperature, °C.
    :return:
        Partial pressure of CO2 at given temperature, kPa or mmHg.
    :rtype: float
    """
    return pCO2 * 10 ** (0.021 * (t - 37))


def calculate_ctO2(pO2, sO2, FCOHb, FMetHb, ctHb):
    """Total oxygen concentration of blood (O2 content).

    Also known as ctO2(B).


    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-35, p. 271, equation 27.

    :param float pO2: kPa
    :param float sO2: fraction
    :param float FCOHb: fraction
    :param float FMetHb: fraction
    :param float ctHb: mmol/L
    :return:
        O2 content, mmol/L.
    :rtype: float
    """
    # May be negative FMetHb replased with zero?
    # if FMetHb < 0:
    #     FMetHb = 0

    # ctHb(mmol/L) == ctHb(g/dL) / 1.61140 == ctHb(g/dL) * 0.62058
    # Vol % == 2.241 (mmol/L)
    # By [1] p. 6-14 or 6-49.
    # pO2 *= 0.133322368  # mmHg to kPa, 1 mmHg = 0.133322368 kPa
    # O2 solubility coefficient in blood at 37 °С.
    alphaO2 = 9.83 * 10 ** -3  # mmol/L/kPa
    return alphaO2 * pO2 + sO2 * (1 - FCOHb - FMetHb) * ctHb


def calculate_pO2_FO2_fraction(pO2, FO2):
    """Oxygen tension ratio of arterial blood and the fraction of oxygen
    in dry inspired air.

    Also known as pO2(a)/FO2(I).

    https://en.wikipedia.org/wiki/Fraction_of_inspired_oxygen#PaO2.2FFiO2_ratio
    Eq. 17

    :param float pO2: kPa
    :param float FO2: Fraction of oxygen in dry inspired air, fraction.
    :return:
        pO2(a)/FO2(I), mmHg.
    :rtype: float
    """
    pO2 /= kPa
    return pO2 / FO2


def calculate_Ca74(pH, Ca):
    """Ionized calcium at pH 7.4.

    References
    ----------

    [1] Radiometer ABL800 Flex Reference Manual English US.
        chapter 6-41, p. 277, equation 45.

    :param float pH: Due to biological variations this function can only be
        used for a pH value in the range 7.2-7.6 [1].
        Radiometer device returns '?' with pH=6.928, Ca=1.62.
    :param float Ca: mmol/L
    :return:
        cCa2+(7.4), mmol/L.
    :rtype: float
    """
    if not 7.2 <= pH <= 7.4:
        raise ValueError(
            "Can calculate only for pH 7.2-7.6 due to biological variations")
    return Ca * (1 - 0.53 * (7.4 - pH))


def expected_pH(pCO2, status='acute'):
    """Calculate expected pH for given pCO2 (chronic or acute patient status).


    References
    ----------

    [1] Kostuchenko S.S., ABB in the ICU, 2009, p. 55.

    :param float pCO2: mmHg
    :return:
        Expected pH.
    :rtype: float
    """
    st = {'acute': 0.008, 'chronic': 0.003}
    return 7.4 + st[status] * (40.0 - pCO2)


def abg(pH, pCO2):
    """Evaluate arterial blood gas status.

    http://en.wikipedia.org/wiki/Arterial_blood_gas

    :param float pH:
    :param float pCO2: mmHg
    :return:
        Opinion.
    :rtype: unicode
    """
    # https://www.kernel.org/doc/Documentation/CodingStyle
    # The answer to that is that if you need more than 3 levels of
    # indentation, you're screwed anyway, and should fix your program.

    def check_metabolic(pH, pCO2):
        """Check metabolic status by expected pH level.

        Does this pH and pCO2 means hidden metabolic process?
        """
        guess = ''
        magic_number = 0.07
        ex_pH = expected_pH(pCO2)
        if abs(pH - ex_pH) > magic_number:
            if pH > ex_pH:
                guess += "background metabolic alcalosis, "
            else:
                guess += "background metabolic acidosis, "
        return "%(guess)sexpected pH %(ex_pH).2f" % {
            'guess': guess, 'ex_pH': ex_pH}

    if norm_pH[0] <= pH <= norm_pH[1]:  # pH is normal or compensated
        # Don't calculating expected CO2/pH values because both values are
        # normal or represent two opposed processes (no need for searching
        # hidden one)
        # pCO2 is abnormal, checking slight pH shifts
        if pCO2 < norm_pCO2[0]:
            # Low (respiratory alcalosis)
            if pH >= 7.41:
                return "Respiratory alcalosis, full comp. by metabolic acidosis"
            else:
                return "Metabolic acidosis, full comp. by CO2 alcalosis"
        elif pCO2 > norm_pCO2[1]:
            # High (respiratory acidosis)
            if pH <= 7.39:  # pH almost acidotic
                # Classic "chronic" COPD gas
                return "Respiratory acidosis, full comp. by metabolic alcalosis. COPD?"
            else:
                return "Metabolic alcalosis, full comp. by CO2 acidosis"
        else:
            return "ABG normal"
    else:
        # pH decompensation
        if pCO2 < norm_pCO2[0]:  # Low (respiratory alcalosis)
            # Достаточно ли такого изменения pH, чтобы дать такой pCO2?
            if pH < norm_pH[0]:
                # Check anion gap here?
                return "Metabolic acidosis, partial comp. by CO2 alcalosis [check BDG]"
            elif pH > norm_pH[1]:
                return "Respiratory alcalosis (%s)" % check_metabolic(pH, pCO2)
        elif pCO2 > norm_pCO2[1]:
            if pH < norm_pH[0]:
                return "Respiratory acidosis (%s)" % check_metabolic(pH, pCO2)
            elif pH > norm_pH[1]:
                return "Metabolic alcalosis, partial comp. by CO2 acidosis [check Na, Cl, albumin]"
        else:
            # Normal pCO2 (35 <= pCO2 <= 45 normal)
            if pH < norm_pH[0]:
                return "Metabolic acidosis, no respiratory comp."
            elif pH > norm_pH[1]:
                return "Metabolic alcalosis, no respiratory comp."


def abg2(pH, pCO2, HCO3=None):
    """Calculate expected ABG values.


    References
    ----------
    [1] Kostuchenko S.S., ABB in the ICU, 2009

    :param float pH:
    :param float pCO2: mmHg
    :param float HCO3: mEq/L, standartized. Evaluated automatically if not
        provided
    :return:
        Opinion.
    :rtype: str
    """
    if HCO3 is None:
        HCO3 = calculate_hco3(pH, pCO2)

    # Assess respiratory problem
    print("y = ΔpH/ΔpCO2×100 = %.2f" % ((7.4 - pH) / (pCO2 - 40.0) * 100))

    # Expected pH (acute, chronic)
    print("pH\t\tby Genderson\texpected %.2f .. %.2f acute-chronic" % (
        7.4 + 0.008 * (40 - pCO2), 7.4 + 0.003 * (40 - pCO2)))

    # Winter's formula (acidosis, alkalosis)
    print("pCO2\tby Winter (x)\texpected %.1f±2 .. %.1f±1.5"
          " acidisis-alkalosis" % (1.5 * HCO3 + 8, 0.7 * HCO3 + 20))


def describe(pH, pCO2):
    HCO3act = calculate_hco3p(pH, pCO2)
    info = "pH %.2f, pCO2 %.2f, HCO3act %.2f, BE %+.2f" % (pH, pCO2, HCO3act, calculate_be(pH, pCO2, HCO3act))
    print(info)
    print(abg(pH, pCO2))
    abg2(pH, pCO2)


def test():
    variants = (
        # pH, pCO2, HCO3, comment
        (7.46, 35., 27., "Metabolic Alkalosis"),
        (7.30, 35., 20., "Metabolic Acidosis"),
        (7.47, 32., 23., "Respiratory Alkalosis"),
        (7.4,  40.01, 25., "Normal"),
        (7.1, 28., 14., "Костюченко 1 вариант"),
        (7.1, 68., None, "Костюченко 2 вариант"),
        (7.3, 68., None, "Костюченко 3 вариант"),
        (7.37, 58., None, "Untagged"),
        (7.07, 53.2, None, "Resp. acidosis + metabolic"),
        (7.39, 47., None, "Resp. acidosis, met. alcalosis, full comp. (COPD)")
    )
    for v in variants:
        describe(pH=v[0], pCO2=v[1])
        print("[%s]\n" % v[3])


if __name__ == '__main__':
    test()
