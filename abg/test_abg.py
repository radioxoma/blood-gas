#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
import math
import pandas as pd
from uncertainties import ufloat
import abg
import odc

kPa = 0.133322368
# kPa = 0.133322  # By Radiometer


def diff(first, second):
    """Is float fits in the uncertainity of ufloat?

    :param ufloat first: Etalon
    :param ufloat first: Calculated value
    """

    sum_dev = first.s + second.s
    if abs(second.std_score(first)) > 1:
        # print('Warning! Bad std_score %s (got %s, expected %s)' % (
        #    second.std_score(first), second, first))
        # print("%f == %f + %f" % (sum_dev, second.s, first.s))
        # print('')
        pass
    if abs(first.n - second.n) <= sum_dev:
        return True
    print('ERROR! std_score %s (got %s, expected %s)' % (
        second.std_score(first), second, first))
    print('Values Δ %s' % (first.n - second.n))
    return False


def main_test():
    data = pd.read_csv('samples.csv')
    good = 0
    err = 0
    for idx, r in data.iterrows():
        # We know input with full precision
        patient = r['Patient']
        FO2 = float(r['FO2']) / 100  # Fraction
        T = r['Temp']

        ID = r['ID']
        ctHb = ufloat(r['ctHb'], 0.1) * 0.62058  # mmol/L
        Hct = ufloat(r['Hct'], 0.1) / 100  # Fraction
        Na = ufloat(r['cNa'], 1)
        Cl = ufloat(r['cCl'], 1)
        pH = ufloat(r['pH'], 0.001)
        SBE = ufloat(r['SBE'], 0.1)
        ABE = ufloat(r['ABE'], 0.1)
        AnionGap = ufloat(r['AnionGap'], 0.1)
        mOsm = ufloat(r['mOsm'], 0.1)
        glucosae = ufloat(r['cGlu'], 0.1)
        pHT = ufloat(r['pHT'], 0.001)

        pCO2T = ufloat(r['pCO2T'], 0.1) * kPa  # mmHg to kPa
        pCO2 = ufloat(r['pCO2'], 0.1) * kPa  # mmHg to kPa
        # pO2 = ufloat(r['pO2'], 0.1) * kPa  # mmHg to kPa
        pO2 = ufloat(r['pO2'], 1) * kPa  # mmHg to kPa, precision may wary
        p50 = ufloat(r['p50'], 0.01) * kPa  # mmHg to kPa
        pO2T = ufloat(r['pO2T'], 0.1) * kPa  # mmHg to kPa

        ctO2 = ufloat(r['ctO2'], 0.1) / 2.241  # mmol/L
        sO2 = ufloat(r['sO2'], 0.1) / 100
        FCOHb = ufloat(r['FCOHb'], 0.1) / 100
        FMetHb = ufloat(r['FMetHb'], 0.1) / 100
        HCO3st = ufloat(r['HCO3st'], 0.1)
        RespIdx = ufloat(r['RespIdx'], 1)

        ## Oxymetry
        # Hct
        if not diff(Hct, abg.calculate_hct(ctHb)):
            print('Hct mismatch')
            print(patient)

        ## Electrolytes
        # Anion gap
        if not diff(AnionGap, abg.calculate_anion_gap(
                Na, Cl, abg.calculate_hco3p(pH, pCO2))):
            print('AnionGap mismatch')
            print(patient)

        # ## Temperature
        # pHT
        if not diff(pHT, abg.calculate_pHT(pH, T)):
            print('pHT mismatch')
            print(patient)

        # pCO2T
        if not diff(pCO2T, abg.calculate_pCO2T(pCO2, T)):
            print('pCO2T mismatch')
            print(patient)

        # pO2T
        # if not diff(pO2T, abg.calculate_pO2T(pO2, t)):
        #     print(patient)

        # O2 status
        # ctO2
        if not diff(ctO2, abg.calculate_ctO2(pO2, sO2, FCOHb, FMetHb, ctHb)):
            print('ctO2 mismatch %s, %s\n' % (patient, ID))

        # p50
        # if sO2 <= 0.97:
        #     odc_model = odc.ODC()
        #     odc_model.fit(sO2=sO2, pO2=pO2 * kPa, pCO2=pCO2, pH=pH,
        #     FCOHb=FCOHb, FMetHb=FMetHb)
        #     if diff(p50, odc_model.eval_p50() / kPa):
        #         print('p50 mismatch %s' % patient)
        #         print(p50, odc_model.eval_p50() / kPa)
        #         err += 1
        #     else:
        #         good += 1

        # pO2(a)/FO2(I)
        if not math.isnan(RespIdx.n):
            if not diff(RespIdx, abg.calculate_pO2_FO2_fraction(pO2, FO2)):
                print('pO2(a)/FO2(I) mismatch')
                print(patient + '\n')
            # else:
                # print("Good pO2(a)/FO2(I)")

        # ## Acid-base status
        # HCO3st
        if not diff(HCO3st, abg.calculate_hco3pst(pH, pCO2, ctHb, sO2)):
            print('HCO3st mismatch')
            print(patient, ID)
        # SBE
        if not diff(SBE, abg.calculate_cbase(pH, pCO2)):
            print('SBE mismatch')
            print(patient)
        # ABE
        if not diff(ABE, abg.calculate_cbase(pH, pCO2, ctHb)):
            print('ABE mismatch')
            print(patient)
        # mOsm
        if not diff(mOsm, abg.calculate_mosm(Na, glucosae)):
            print('mOsm mismatch')
            print(patient)
    # print("Good %d / err %d, total: %d" % (good, err, good + err))


if __name__ == '__main__':
    main_test()
