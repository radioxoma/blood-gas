#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Eugene Dvoretsky 2015, Vitebsk state medical university.

Рассмотреть вариант вынос ac в отдельную функцию.
p50,c   `a`
p50(st) `a6`, т.к. ac при НУ == 0.
pO2(T)  `ac` с модификацией

Not shure about parameters "by default": Radiometer sometimes declares it
as zero (e.g. ODC reference position), sometimes no (default values).
"""

from __future__ import absolute_import
from __future__ import division
try:
    from uncertainties import umath as math
except ImportError:
    import math

epsilon = 0.0001  # Precision of Newton-Raphson algorithm

norm_p50 = (24, 28)  # mmHg or ~26.6
norm_p50st = norm_p50
# p50st = 3.578  # kPa (26.84 mmHg)
k_0 = 0.5343
h_0 = 3.5
T_0 = 37  # °C
s_0 = 0.867
p_00 = 7  # kPa

FHbF = 0  # Указан только в одном месте, а там все нули.
cDPG = 5  # mmol\L
# FCOHb = 0.004  # Default
# FMetHb = 0.004  # Default
# FCOHb = 0  # Standard
# FMetHb = 0  # Standard


class ODC(object):

    """Oxygemoglobin dissiciation curve (ODC) model.

    .. warning:: This all about reverse-engineering. I can't guarantee that
        algorithm exactly the same as used by Radiometer ABL800 Flex analyzer.
        You should not rely on this code for making life-threatening decisions.

    Usage:

        * Make instance of the class.
        * `fit()` curve for data measured in blood of specific patient.
        * Now you can use other object methods.

    If you interesting about the ABG machine internals, these sources give you
    some inspiration:

        * Check [1] chapter 6-44, p. 280, equation 46-47 for basic curve model.
        * [2] is fundamental paper. Read it!
        * Read [3] for insight about Radiometer [1] and this class internals.


    References
    ----------

    .. [1] Radiometer ABL800 Flex Reference Manual English US.
    .. [2] Clin Lab Invest 1990; 50, Suppl 203: 75-86. Available as AS107. 16.
        Siggaard-Andersen O, Wimberley PD, Gøthgen IH, Siggaard-Andersen M. A
        mathematical model of the hemoglobin-oxygen dissociation curve of human
        blood and of the oxygen partial pressure as a function of temperature.
        Clin Chem 1984; 30: 1646-51.
    .. [3] http://www.derangedphysiology.com/php/Arterial-blood-gases/p50.php
    """

    def fit(
            self, sO2, pO2, pCO2, pH,
            T=37, FCOHb=0.004, FMetHb=0.004, p50st=None):
        """Evaluate ODC position for specific blood sample.


        .. note:: The reference position of the ODC was chosen to be the one
        that corresponds to the default value for p50(st) = 3.578 kPa,
        which is traditionally considered the most likely value of p50 for
        adult humans under standard conditions, namely:

                pH = 7.40
                pCO2 = 5.33  # kPa
                FCOHb, FMetHb, FHbF = (0, 0, 0)
                cDPG = 5  # mmol/L

        Pass measured parameters to this function.


        .. note:: For some devices [1]:
            If the sO2 value for establishing the ODC is greater than 0.97,
            the calculation of the some parameter is not performed unless
            the p50(st) value is keyed in.

        Enter measured parameters to fit curve in it.
        """
        self.FCOHb = FCOHb
        self.FMetHb = FMetHb
        self.sO2 = sO2
        self.pO2 = pO2
        self.pCO2 = pCO2
        self.pH = pH
        self.T = T
        self.p50st = p50st

        # Eq. 46.3
        self.y_0 = math.log(s_0 / (1 - s_0))
        #######################################################################
        # Determining actual displacement 'a'. It includes:
        #     * 'ac' - an guess, by measured parameters
        #     * 'a6' - an additional shift
        a1 = -0.88 * (pH - 7.40)
        a2 = 0.048 * math.log(pCO2 / 5.33)  # 5.33 pCO2
        a3 = -0.7 * FMetHb
        a4 = (0.06 - 0.02 * FHbF) * (cDPG - 5)
        a5 = -0.25 * FHbF
        ac = a1 + a2 + a3 + a4 + a5

        if sO2 <= 0.97 and not p50st:  # I
            # Рассчитать P0S0 по измеренным значениям
            # На основе измеренных параметров рассчитать сдвиг 'ac'
            # Использовать 46.3, 46.4?
            # Определить 'a6' кривой reference position при которых она
            # будет проходить через измеренную точку P0S0
            P0 = pO2 + (pO2 / sO2) * (FCOHb / (1 - FCOHb - FMetHb))  # 46.9
            x_measured = math.log(P0)
            S0 = (sO2 * (1 - FCOHb - FMetHb) + FCOHb) / (1 - FMetHb)  # 46.11
            y_measured = math.log(S0 / (1 - S0))
            # print('P0', P0, 'S0', S0)
            # print('x_measured', x_measured, 'y_measured', y_measured)

            # Newtom-Rapson method
            # http://web.mit.edu/10.001/Web/Course_Notes/NLAE/node6.html
            a = 0  # Start value, as described in paper
            while True:
                x_0i = eval_x_0(a=a, T=T)
                # n ~ 2.7 according to paper
                y_i = haldane_odc(x=x_measured, x_0=x_0i, y_0=self.y_0, a=a)
                if abs(y_measured - y_i) < epsilon:
                    break
                n = haldane_odc_diff(x=x_measured, x_0=x_0i, y_0=self.y_0, a=a)
                a = a + (y_measured - y_i) / (
                    -n + math.tanh(k_0 * (x_measured - x_0i)))  # Brackets
            self.ac = ac
            self.a = a
            self.a6 = a - ac
        else:
            if p50st is not None:  # II
                # Experimental for pO2(T)?
                # 46.9, sO2 = 0.5
                P0 = p50st + (p50st / 0.5) * (FCOHb / (1 - FCOHb - FMetHb))
                x_measured = math.log(P0)
                # 46.11
                S0 = (0.5 * (1 - FCOHb - FMetHb) + FCOHb) / (1 - FMetHb)
                y_measured = math.log(S0 / (1 - S0))
                # *Кривая при стандартных условиях p50st для данного пациента*
                # Рассчитать точку P0S0 по давлению p50st, сатурации 0.5
                # Итеративно определаить `a6` (без учёта ac) при котором кривая
                #     проходиn через рассчитанную точку P0S0
                a6 = 0  # Start value, as described in paper
                while True:
                    x_0i = eval_x_0(a=a6, T=T)
                    # n ~ 2.7 according to paper
                    y_i = haldane_odc(
                        x=x_measured, x_0=x_0i, y_0=self.y_0, a=a6)
                    if abs(y_measured - y_i) < epsilon:
                        break
                    n = haldane_odc_diff(
                        x=x_measured, x_0=x_0i, y_0=self.y_0, a=a6)
                    a6 = a6 + (y_measured - y_i) / (
                        -n + math.tanh(k_0 * (x_measured - x_0i)))  # Brackets
                # *Расчёт кривой p50act*
                # К рассчитанному для p50st `a6` прибавить рассчитанный по
                # измеряемым параметрам сдвиг 'ac'
                # So `a6` it's shift from reference to keyed standard cond.
                # `ac` is shift from standard conditions to patient body cond.
                self.a6 = a6
                self.ac = ac
                self.a = a6 + ac
            else:  # III, ошибочный или зашкаливающий sO2, p50st неизвестно.
                # Из-за зашкалиающего pO2 кривая будет расчитана приблизительно
                # Рассчитать по измеряемым параметрам (pH, pCO2, FCOHb, FMetHb,
                #    FHbF) сдвиг 'ac'
                # Кривая пациента приблизительно соответсвует reference-кривой,
                #    сдвинутой на рассчитанный 'ac'
                # a = ac  # `a6` не нужно определять
                self.ac = ac
                self.a = ac
                self.a6 = 0

    def fit_standard(self, p50st=3.578, *args, **kwargs):
        """Not shure about *args/**kwargs trick.

        p50st from 6-30, p. 266

        Fixme: if p50st given, no need for patient sO2, pO2.
        """
        self.fit(*args, p50st=p50st, **kwargs)

    def eval_pressure(self, sO2, A, T):
        """Calculate O2 pressure by saturation.

        P = ODC(S,A,T)

        [16] Siggaard-Andersen O, Wimberley PD, Göthgen IH,
        Siggaard-Andersen M.
        A mathematical model of the hemoglobin-oxygen dissociation curve
        of human blood and of the oxygen partial pressure as a function
        of temperature. Clin Chem 1984; 30: 1646-51.

        [18] Siggaard-Andersen O, Siggaard-Andersen M.
        The oxygen status algorithm: a computer program for calculating and
        displaying pH and blood gas data.
        Scand J Clin Lab Invest 1990; 50, Suppl 203: 29-45.

        :param float sO2: measured hemoglobin saturation, fraction (46.2)
        :param float A: an curve displacement along axis x (46.5).
        :param float T: Body temperature, °C (46.7).
        :return:
            p, partial O2 pressure at given sO2 and T conditions for current
            ODC, kPa. No hemoglobin corrections performed.
        :rtype: float
        """
        # 46.2
        y = math.log(sO2 / (1 - sO2))
        # Newtom-Rapson iterative method
        x_0 = eval_x_0(a=A, T=T)
        x = x_0  # Start value, as described in paper
        while True:
            y_i = haldane_odc(x=x, x_0=x_0, y_0=self.y_0, a=A)
            if abs(y - y_i) < epsilon:
                # print("FOUND x", x)
                break
            # n ~ 2.7 according to paper
            n = haldane_odc_diff(x, x_0, self.y_0, A)
            x = x + (y - y_i) / n
            # print(y, y_i)
        # print('y', y, 'x', x)
        # print('\n%s y ~ \n%s y_hal\n' % (
        #     y, haldane_odc(x, x_0, self.y_0, A)))
        p = math.exp(x)  # Reverse 46.1
        return p

    def eval_saturation(self, pO2, A, T):
        """Calculate saturation by curve.

        (Saturation can be measured by multiwavelength hemoximetry.)

        S = ODC(P,A,T)


        References
        ----------

        .. [1] Radiometer ABL800 Flex Reference Manual English US.
            chapter 6-44, p. 280, equation 46-47.

        :param float pO2: measured partial O2 pressure, kPa (46.1)
        :param float A: an curve displacement along axis x (46.5).
        :param float T: Body temperature, °C (46.7).
        :return:
            s, hemoglobin saturation at given pO2 and T conditions for current
            ODC, fraction. No hemoglobin corrections performed.
        :rtype: float
        """
        x_0 = eval_x_0(a=A, T=T)
        # 46.1
        x = math.log(pO2)
        y = haldane_odc(x=x, x_0=x_0, y_0=self.y_0, a=A)
        s = 1 / (math.exp(-y) + 1)  # Reverse 46.2
        return s

    def eval_p50(self):
        """Partial pressure of oxygen at half saturation (sO2 50 %) in blood.

        If sO2 > 97 % for Siggaard-Andersen Oxygen Status Algorithm can be
        non reliable.


        References
        ----------

        .. [1] Radiometer ABL800 Flex Reference Manual English US.
            chapter 6-32, p. 268, equation 19.

        :return:
            p50, kPa.
        :rtype: float
        """
        S = (0.5 * (1 - self.FCOHb - self.FMetHb) + self.FCOHb) / (
            1 - self.FMetHb)
        P = self.eval_pressure(sO2=S, A=self.a, T=37)
        return P / (1 + (self.FCOHb / 0.5 * (1 - self.FCOHb - self.FMetHb)))

    def eval_p50st(self):
        """Partial pressure of oxygen at half saturation (sO2 50 %) in blood
        and standard conditions.

        Not tested with sO2 > 97 %.
        p50st can be input or derived parameter.


        References
        ----------

        .. [1] Radiometer ABL800 Flex Reference Manual English US.
            chapter 6-33, p. 269, equation 21.

        :return:
            p50(st), kPa.
        :rtype: float
        """
        # If `a = ac + a6` (`ac == 0` at standard conditions), then `a == a6`
        # pH = 7.4; FCOHb = 0; FMetHb = 0; FHbF = 0; pCO2 = 5.33  # kPa
        return self.eval_pressure(sO2=0.5, A=self.a6, T=37)

    def eval_pO2T(self, ctHb, T):
        # FCOHb=0.004, FMetHb=0.004
        """Not implemented yet. Fixme: must return same O2 at 37 degrees.

        pO2 of blood at patient temperature.

        Standard ODC shifted with patient values is used.


        Возможно некоторые параметры жёстко привязаны к текущей модели кривой.

        Note that pO2T value at 37 °C must be the same as pO2.


        References
        ----------

        .. [1] Radiometer ABL800 Flex Reference Manual English US.
            chapter 6-30, p. 266, equation 14.

        :param float ctHb:
        :param float T: Body temperature, °C.
        :return:
            pO2(T), kPa
        :rtype: float
        """
        def calc_tiT(pO2, sO2, T):
            """O2 content. Based on eq. 19 from paper.
            """
            alphaO2 = 9.83 * 10 ** -3 * math.exp(
                -1.15 * 10 ** -2 * (T - 37) + 2.1 * 10 ** -4 * (T - 37) ** 2)
            return ctHb * (1 - self.FCOHb - self.FMetHb) * sO2 + alphaO2 * pO2

        # def calc_tiT_diff(pO2, sO2, T):
        #     # I'm not shure if it is a derivative
        #     alphaO2 = 9.83 * 10 ** -3 * math.exp(
        #         -1.15 * 10 ** -2 * (T - 37) + 2.1 * 10 ** -4 * (T - 37) ** 2)
        #     return 1 + alphaO2  # * pO2

        def calc_tiT_diff(pO2, sO2, T, A):
            # Derivative from paper
            alphaO2 = 9.83 * 10 ** -3 * math.exp(
                -1.15 * 10 ** -2 * (T - 37) + 2.1 * 10 ** -4 * (T - 37) ** 2)
            x = math.log(pO2)
            x_0 = eval_x_0(a=A, T=T)
            n = haldane_odc_diff(x=x, x_0=x_0, y_0=self.y_0, a=A)
            return alphaO2 + ctHb * n * (1 - sO2) / pO2

        dpHdT = -1.46 * 10 ** -2 - 6.5 * 10 ** -3 * (self.pH - 7.4)
        A_37 = self.ac - 1.04 * dpHdT * (T - 37)
        P_37 = self.pO2 + (self.pO2 / self.sO2) * (self.FCOHb / (
            1 - self.FCOHb - self.FMetHb))  # 46.9
        S_37 = self.eval_saturation(pO2=P_37, A=A_37, T=37)
        # P_37 = self.pO2
        # S_37 = self.sO2
        t_37 = calc_tiT(pO2=P_37, sO2=S_37, T=37)
        Ai = self.ac - 1.04 * dpHdT * (T - 37)
        P = 3  # kPa, by paper
        counter = 0
        while True:
            counter += 1
            Si = self.eval_saturation(pO2=P, A=Ai, T=T)
            sO2iT = (Si * (1 - self.FMetHb) - self.FCOHb) / (
                1 - self.FCOHb - self.FMetHb)
            pO2iT = P / (
                1 + (self.FCOHb / (sO2iT * (1 - self.FCOHb - self.FMetHb))))
            tiT = calc_tiT(pO2=pO2iT, sO2=sO2iT, T=T)
            # print(t_37 - tiT)
            if (t_37 - tiT) < epsilon:
                break
            n = calc_tiT_diff(pO2=pO2iT, sO2=sO2iT, T=T, A=Ai)
            P = P + (t_37 - tiT) / n

            if counter > 50:
                break
                raise ValueError("Infinite iteration")
        return pO2iT

    # def test_pO2T(self, ctHb, T):
    #     P_37 = self.pO2 + (self.pO2 / self.sO2) * (self.FCOHb / (
    #         1 - self.FCOHb - self.FMetHb))  # 46.9
    #     sO2iT = self.sO2
    #     print(self.pO2, self.sO2)
    #     return P_37 / (
    #         1 + (self.FCOHb / (sO2iT * (1 - self.FCOHb - self.FMetHb))))


def eval_x_0(a, T):
    """Will be calculated multiple times to allow other functions get
    temperature as parameter.

    :param float T: Celsus temperature.
    Return `x_0` for `haldane_odc` or `haldane_odc_diff`.
    """
    b = 0.055 * (T - T_0)  # Eq. 46.7, temperature shift along `x`
    return math.log(p_00) + a + b  # Eq. 46.4


def haldane_odc(x, x_0, y_0, a):
    """Oxygen dissotiation curve equation.

    Return `y` coordinate (saturation).
    """
    h = h_0 + a  # Eq. 46.6
    return y_0 + (x - x_0) + h * math.tanh(k_0 * (x - x_0))  # Eq. 46


def haldane_odc_diff(x, x_0, y_0, a):
    """Differential of `haldane_odc` (dy/dx Hill slope).
    """
    h = h_0 + a  # Eq. 46.6
    return 1 + h * k_0 * (1 - math.tanh(k_0 * (x - x_0)) ** 2)


def main_test():
    # sO2 = 97.2 / 100
    # pO2 = 63.5 * 0.133322368
    # pCO2 = 63.5 * 0.133322368
    # pH = 7.390
    # T = 37
    # FCOHb = 5.2 / 100
    # FMetHb = -0.5 / 100
    # ctHb = 11.9 * 0.62058  # mmol/L
    # # tiT ~ 7.261  # 37 celsus
    # pO2T = 63.5

    # Иванова
    # sO2 = 45.3 / 100
    # pO2 = 33.7 * 0.133322368
    # pCO2 = 68.6 * 0.133322368
    # pH = 6.919
    # T = 39.6
    # FCOHb = 1.6 / 100
    # FMetHb = 0.7 / 100
    # ctHb = 12.3 * 0.62058  # mmol/L
    # pCO2T = 77.8 * 0.133322368  # 10.372 kPa
    # pO2T = 40.3 * 0.133322368

    # Неизвестная 1
    sO2 = 100.3 / 100
    # pO2 = 454 * 0.133322368
    pO2 = 454 * 0.133322368
    pCO2 = 26.2 * 0.133322368
    pH = 6.945
    T = 37
    FCOHb = 2.8 / 100
    FMetHb = -0.2 / 100
    ctHb = 11.4 * 0.62058  # mmol/L
    pCO2T = 26.2 * 0.133322368
    pO2T = 454 * 0.133322368

    odc = ODC()
    # odc.fit(sO2=sO2, pO2=pO2, pCO2=pCO2, pH=pH, FCOHb=FCOHb, FMetHb=FMetHb)
    odc.fit_standard(sO2=sO2, pO2=pO2, pCO2=pCO2, pH=pH, FCOHb=FCOHb, FMetHb=FMetHb)

    # print("--> eval_odc_saturation = %s" % eval_odc_saturation(
    #     sO2=sO2, pO2=pO2, pCO2=pCO2, pH=pH, T=T, FCOHb=FCOHb, FMetHb=FMetHb))
    # print("Expected: y = 0.910810\n\n")

    # p50 = odc.eval_p50()
    # assert(round(p50, 6) == round(3.51851472023, 6))
    # print("calculate_p50 %s kPa (%s mmHg)" % (p50, p50 / 0.133322368))
    # print("Expected                         26.3910308001 mmHg")
    # print('eval_p50st %f' % (odc.eval_p50st() / 0.133322368))
    print('eval_pO2T %f' % (odc.eval_pO2T(ctHb=ctHb, T=T) / 0.133322368))
    print("expected  %s" % (pO2T / 0.133322368))

    # print("Test expected %s" % (odc.pO2 / 0.133322368))
    # print('test_pO2T %f' % (odc.test_pO2T(ctHb=ctHb, T=T) / 0.133322368))


if __name__ == '__main__':
    main_test()
