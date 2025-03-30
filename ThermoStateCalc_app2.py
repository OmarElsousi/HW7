import sys
from ThermoStateCalc2 import Ui__frm_StateCalculator
from pyXSteam.XSteam import XSteam
from PyQt5.QtWidgets import QWidget, QApplication
from UnitConversion import UC
from scipy.optimize import fsolve


class thermoState:
    def __init__(self, p=None, t=None, v=None, u=None, h=None, s=None, x=None):
        """
        This is a class I use for storing a thermodynamic state. Calling setState requires you
        to specify two independent thermodynamic properties. One ambiguity exists if you specify
        both psat and tsat. In that case I assume two-phase with x=0.5.

        :param p:
        :param t:
        :param v:
        :param u:
        :param h:
        :param s:
        :param x:
        """
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.region = "saturated"
        self.p = p
        self.t = t
        self.v = v
        self.u = u
        self.h = h
        self.s = s
        self.x = x

    def getSaturationProperties(self, p):
        """Return saturated vf, vg, hf, hg, sf, sg at pressure p."""
        vf = self.steamTable.vL_p(p)
        vg = self.steamTable.vV_p(p)
        hf = self.steamTable.hL_p(p)
        hg = self.steamTable.hV_p(p)
        sf = self.steamTable.sL_p(p)
        sg = self.steamTable.sV_p(p)
        return vf, vg, hf, hg, sf, sg

    def computeProperties(self):
        """
        Assumes p and t are already calculated (or region is identified).
        If region == "two-phase", we interpolate based on quality x.
        Else, we look up single-phase properties directly from the steam tables.
        """
        if self.region == "two-phase":
            # Interpolate based on quality
            self.u = self.steamTable.uL_p(self.p) + self.x * (
                self.steamTable.uV_p(self.p) - self.steamTable.uL_p(self.p)
            )
            self.h = self.steamTable.hL_p(self.p) + self.x * (
                self.steamTable.hV_p(self.p) - self.steamTable.hL_p(self.p)
            )
            self.s = self.steamTable.sL_p(self.p) + self.x * (
                self.steamTable.sV_p(self.p) - self.steamTable.sL_p(self.p)
            )
            self.v = self.steamTable.vL_p(self.p) + self.x * (
                self.steamTable.vV_p(self.p) - self.steamTable.vL_p(self.p)
            )
        else:
            # Single-phase: use p-T directly
            self.u = self.steamTable.u_pt(self.p, self.t)
            self.h = self.steamTable.h_pt(self.p, self.t)
            self.s = self.steamTable.s_pt(self.p, self.t)
            self.v = self.steamTable.v_pt(self.p, self.t)
            # For convenience, set quality to 1.0 if "super-heated vapor"; otherwise 0.0
            self.x = 1.0 if self.region == "super-heated vapor" else 0.0

    def setState(self, stProp1, stProp2, stPropVal1, stPropVal2, SI=True):
        """
        Calculates the thermodynamic state variables based on two specified properties.
        We have p, t, v, h, u, s, x as possible properties, giving 21 distinct pairings.

        The big block of if/elif logic below detects which pair of properties we're given,
        figures out whether it's sub-cooled, superheated, or two-phase, sets self.region,
        and calculates self.p, self.t, self.x, etc. At the end, we call computeProperties().
        """
        # Switch between SI or English
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)

        # For clarity, store lower-cased property labels
        SP = [stProp1.lower(), stProp2.lower()]
        f1 = float(stPropVal1)
        f2 = float(stPropVal2)

        # ----------------------------------------------------------------------
        # The logic below is taken from your original code, with indentation fixes
        # and duplicated lines removed. Each property combination is handled.
        # ----------------------------------------------------------------------

        if SP[0] == 'p' or SP[1] == 'p':
            oFlipped = SP[0] != 'p'
            SP1 = SP[0] if oFlipped else SP[1]
            self.p = f1 if not oFlipped else f2
            tSat = self.steamTable.tsat_p(self.p)
            # case 1: p & t
            if SP1 == 't':
                self.t = f2 if not oFlipped else f1
                epsilon = 1e-3
                if self.t >= tSat + epsilon:
                    self.region = "super-heated vapor"
                elif self.t <= tSat - epsilon:
                    self.region = "sub-cooled liquid"
                else:
                    self.region = "saturated"

            # case 2: p & v
            elif SP1 == 'v':
                self.v = f2 if not oFlipped else f1
                vf = round(self.steamTable.vL_p(self.p), 5)
                vg = round(self.steamTable.vV_p(self.p), 5)
                vgf = vg - vf
                if self.v < vf or self.v > vg:
                    self.region = "sub-cooled liquid" if self.v < vf else "super-heated vapor"
                    # Solve for T
                    dt = 1.0 if self.v > vg else -1.0
                    fn_v = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn_v, [tSat + dt])[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.v - vf) / vgf
                    self.t = tSat

            # case 3: p & u
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                uf = round(self.steamTable.uL_p(self.p), 5)
                ug = round(self.steamTable.uV_p(self.p), 5)
                ugf = ug - uf
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled liquid" if self.u < uf else "super-heated vapor"
                    dt = 1.0 if self.u > ug else -1.0
                    fn_u = lambda T: self.u - self.steamTable.u_pt(self.p, T)
                    self.t = fsolve(fn_u, [tSat + dt])[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.u - uf) / ugf
                    self.t = tSat

            # case 4: p & h
            elif SP1 == 'h':
                self.h = f2 if not oFlipped else f1
                hf = self.steamTable.hL_p(self.p)
                hg = self.steamTable.hV_p(self.p)
                hgf = hg - hf
                if self.h < hf or self.h > hg:
                    self.region = "sub-cooled liquid" if self.h < hf else "super-heated vapor"
                    self.t = self.steamTable.t_ph(self.p, self.h)
                else:
                    self.region = "two-phase"
                    self.x = (self.h - hf) / hgf
                    self.t = tSat

            # case 5: p & s
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                sgf = sg - sf
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                    self.t = self.steamTable.t_ps(self.p, self.s)
                else:
                    self.region = "two-phase"
                    self.x = (self.s - sf) / sgf
                    self.t = tSat

            # case 6: p & x
            elif SP1 == 'x':
                self.region = "two-phase"
                self.x = f2 if not oFlipped else f1
                self.t = tSat

        elif SP[0] == 't' or SP[1] == 't':
            oFlipped = SP[0] != 't'
            SP1 = SP[0] if oFlipped else SP[1]
            self.t = f1 if not oFlipped else f2
            pSat = self.steamTable.psat_t(self.t)

            # case 7: t & v
            if SP1 == 'v':
                self.v = f2 if not oFlipped else f1
                vf = self.steamTable.vL_p(pSat)
                vg = self.steamTable.vV_p(pSat)
                vgf = vg - vf
                if self.v < vf or self.v > vg:
                    self.region = "sub-cooled liquid" if self.v < vf else "super-heated vapor"
                    dp = -0.1 if self.v > vg else 0.1
                    fn_tv = lambda P: self.v - self.steamTable.v_pt(P, self.t)
                    self.p = fsolve(fn_tv, [pSat + dp])[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.v - vf) / vgf
                    self.p = pSat

            # case 8: t & u
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                uf = self.steamTable.uL_p(pSat)
                ug = self.steamTable.uV_p(pSat)
                ugf = ug - uf
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled liquid" if self.u < uf else "super-heated vapor"
                    dp = 0.1 if self.u > ug else -0.1
                    fn_tu = lambda P: self.u - self.steamTable.u_pt(P, self.t)
                    self.p = fsolve(fn_tu, [pSat + dp])[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.u - uf) / ugf
                    self.p = pSat

            # case 9: t & h
            elif SP1 == 'h':
                self.h = f2 if not oFlipped else f1
                hf = self.steamTable.hL_p(pSat)
                hg = self.steamTable.hV_p(pSat)
                hgf = hg - hf
                if self.h < hf or self.h > hg:
                    self.region = "sub-cooled liquid" if self.h < hf else "super-heated vapor"
                    self.p = self.steamTable.p_th(self.t, self.h)
                else:
                    self.region = "two-phase"
                    self.x = (self.h - hf) / hgf
                    self.p = pSat

            # case 10: t & s
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                sf = self.steamTable.sL_p(pSat)
                sg = self.steamTable.sV_p(pSat)
                sgf = sg - sf
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                    self.p = self.steamTable.p_ts(self.t, self.s)
                else:
                    self.region = "two-phase"
                    self.x = (self.s - sf) / sgf
                    self.p = pSat

            # case 11: t & x
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                self.p = pSat

        elif SP[0] == 'v' or SP[1] == 'v':
            oFlipped = SP[0] != 'v'
            SP1 = SP[0] if oFlipped else SP[1]
            self.v = f1 if not oFlipped else f2

            # case 12: v & h
            if SP1 == 'h':
                self.h = f2 if not oFlipped else f1

                def fn_vh(P):
                    hf = self.steamTable.hL_p(P)
                    hg = self.steamTable.hV_p(P)
                    hgf = hg - hf
                    vf = self.steamTable.vL_p(P)
                    vg = self.steamTable.vV_p(P)
                    vgf = vg - vf
                    # if h in [hf, hg], we can have two-phase, otherwise single-phase
                    if hf <= self.h <= hg:
                        x_ = (self.h - hf) / hgf
                        return self.v - (vf + x_ * vgf)
                    else:
                        # single-phase approach
                        return self.v - self.steamTable.v_ph(P, self.h)

                self.p = fsolve(fn_vh, [1.0])[0]
                tsat = self.steamTable.tsat_p(self.p)
                vf = self.steamTable.vL_p(self.p)
                vg = self.steamTable.vV_p(self.p)
                if self.v < vf:
                    self.region = "sub-cooled liquid"
                    dt = -1
                    fn_t = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn_t, [tsat + dt])[0]
                elif self.v > vg:
                    self.region = "super-heated vapor"
                    dt = 1
                    fn_t2 = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn_t2, [tsat + dt])[0]
                else:
                    self.region = "two-phase"
                    vgf = vg - vf
                    hf = self.steamTable.hL_p(self.p)
                    hg = self.steamTable.hV_p(self.p)
                    hgf = hg - hf
                    self.t = tsat
                    xh = (self.h - hf) / hgf
                    xv = (self.v - vf) / vgf
                    # ideally xh == xv, we could pick one or average
                    self.x = xh

            # case 13: v & u
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1

                def fn_vu(PT):
                    # We'll pass a guess for (p, t)
                    p_, t_ = PT
                    uf = self.steamTable.uL_p(p_)
                    ug = self.steamTable.uV_p(p_)
                    vf = self.steamTable.vL_p(p_)
                    vg = self.steamTable.vV_p(p_)
                    # two-phase check
                    if uf <= self.u <= ug:
                        # we might be in two-phase or single-phase
                        x_ = (self.u - uf) / (ug - uf)
                        v_calc = vf + x_ * (vg - vf)
                        return [self.v - v_calc, 0.0]  # second eq doesn't matter if two-phase
                    else:
                        # single-phase approach
                        return [
                            self.v - self.steamTable.v_pt(p_, t_),
                            self.u - self.steamTable.u_pt(p_, t_),
                        ]

                # Start guess
                props = fsolve(fn_vu, [1.0, 100.0])
                self.p, self.t = props[0], props[1]
                uf = self.steamTable.uL_p(self.p)
                ug = self.steamTable.uV_p(self.p)
                if self.u < uf:
                    self.region = "sub-cooled liquid"
                elif self.u > ug:
                    self.region = "super-heated vapor"
                else:
                    self.region = "two-phase"
                    x_ = (self.u - uf) / (ug - uf)
                    self.x = x_

            # case 14: v & s
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1

                def fn_vs(PT):
                    p_, t_ = PT
                    sf = self.steamTable.sL_p(p_)
                    sg = self.steamTable.sV_p(p_)
                    vf = self.steamTable.vL_p(p_)
                    vg = self.steamTable.vV_p(p_)
                    if sf <= self.s <= sg:
                        x_ = (self.s - sf) / (sg - sf)
                        v_calc = vf + x_ * (vg - vf)
                        return [
                            self.v - v_calc,
                            0.0
                        ]
                    return [
                        self.v - self.steamTable.v_pt(p_, t_),
                        self.s - self.steamTable.s_pt(p_, t_)
                    ]

                props = fsolve(fn_vs, [1.0, 100.0])
                self.p, self.t = props[0], props[1]
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                if self.s < sf:
                    self.region = "sub-cooled liquid"
                elif self.s > sg:
                    self.region = "super-heated vapor"
                else:
                    self.region = "two-phase"
                    self.x = (self.s - sf) / (sg - sf)

            # case 15: v & x
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.region = "two-phase"

                def fn_vx(p_):
                    vf = self.steamTable.vL_p(p_)
                    vg = self.steamTable.vV_p(p_)
                    return self.v - (vf + self.x * (vg - vf))

                self.p = fsolve(fn_vx, [1.0])[0]
                self.t = self.steamTable.tsat_p(self.p)

        elif SP[0] == 'h' or SP[1] == 'h':
            oFlipped = SP[0] != 'h'
            SP1 = SP[0] if oFlipped else SP[1]
            self.h = f1 if not oFlipped else f2

            # case 16: h & u
            if SP1 == 'u':
                self.u = f2 if not oFlipped else f1

                def fn_hu(PT):
                    p_, t_ = PT
                    uf = self.steamTable.uL_p(p_)
                    ug = self.steamTable.uV_p(p_)
                    hf = self.steamTable.hL_p(p_)
                    hg = self.steamTable.hV_p(p_)
                    if uf <= self.u <= ug:
                        x_ = (self.u - uf) / (ug - uf)
                        h_calc = hf + x_ * (hg - hf)
                        return [self.h - h_calc, 0.0]
                    return [
                        self.h - self.steamTable.h_pt(p_, t_),
                        self.u - self.steamTable.u_pt(p_, t_)
                    ]

                props = fsolve(fn_hu, [1.0, 100.0])
                self.p, self.t = props[0], props[1]
                uf = self.steamTable.uL_p(self.p)
                ug = self.steamTable.uV_p(self.p)
                if self.u < uf:
                    self.region = "sub-cooled liquid"
                elif self.u > ug:
                    self.region = "super-heated vapor"
                else:
                    self.region = "two-phase"
                    x_ = (self.u - uf) / (ug - uf)
                    self.x = x_

            # case 17: h & s
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1

                def fn_hs(PT):
                    p_, t_ = PT
                    sf = self.steamTable.sL_p(p_)
                    sg = self.steamTable.sV_p(p_)
                    hf = self.steamTable.hL_p(p_)
                    hg = self.steamTable.hV_p(p_)
                    if sf <= self.s <= sg:
                        x_ = (self.s - sf) / (sg - sf)
                        h_calc = hf + x_ * (hg - hf)
                        s_calc = sf + x_ * (sg - sf)
                        return [self.h - h_calc, self.s - s_calc]
                    return [
                        self.h - self.steamTable.h_pt(p_, t_),
                        self.s - self.steamTable.s_pt(p_, t_)
                    ]

                props = fsolve(fn_hs, [1.0, 100.0])
                self.p, self.t = props[0], props[1]
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                if self.s < sf:
                    self.region = "sub-cooled liquid"
                elif self.s > sg:
                    self.region = "super-heated vapor"
                else:
                    self.region = "two-phase"
                    x_ = (self.s - sf) / (sg - sf)
                    self.x = x_

            # case 18: h & x
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.region = "two-phase"

                def fn_hx(p_):
                    hf = self.steamTable.hL_p(p_)
                    hg = self.steamTable.hV_p(p_)
                    return self.h - (hf + self.x * (hg - hf))

                self.p = fsolve(fn_hx, [1.0])[0]
                self.t = self.steamTable.tsat_p(self.p)

        elif SP[0] == 'u' or SP[1] == 'u':
            oFlipped = SP[0] != 'u'
            SP1 = SP[0] if oFlipped else SP[1]
            self.u = f1 if not oFlipped else f2

            # case 19: u & s
            if SP1 == 's':
                self.s = f2 if not oFlipped else f1

                def fn_us(PT):
                    p_, t_ = PT
                    sf = self.steamTable.sL_p(p_)
                    sg = self.steamTable.sV_p(p_)
                    uf = self.steamTable.uL_p(p_)
                    ug = self.steamTable.uV_p(p_)
                    if sf <= self.s <= sg:
                        x_ = (self.s - sf) / (sg - sf)
                        # we won't try to solve 2-phase with that eq alone, but let's try
                        # single-phase approach next
                        return [
                            self.u - (uf + x_*(ug - uf)),
                            0.0
                        ]
                    return [
                        self.u - self.steamTable.u_pt(p_, t_),
                        self.s - self.steamTable.s_pt(p_, t_)
                    ]

                props = fsolve(fn_us, [1.0, 100.0])
                self.p, self.t = props[0], props[1]
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                if self.s < sf:
                    self.region = "sub-cooled liquid"
                elif self.s > sg:
                    self.region = "super-heated vapor"
                else:
                    self.region = "two-phase"
                    x_ = (self.s - sf) / (sg - sf)
                    uf = self.steamTable.uL_p(self.p)
                    ug = self.steamTable.uV_p(self.p)
                    # if needed, compare self.u with uf + x_*(ug-uf)
                    self.x = x_

            # case 20: u & x
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.region = "two-phase"

                def fn_ux(p_):
                    uf = self.steamTable.uL_p(p_)
                    ug = self.steamTable.uV_p(p_)
                    return self.u - (uf + self.x * (ug - uf))

                self.p = fsolve(fn_ux, [1.0])[0]
                self.t = self.steamTable.tsat_p(self.p)

        elif SP[0] == 's' or SP[1] == 's':
            oFlipped = SP[0] != 's'
            SP1 = SP[0] if oFlipped else SP[1]
            self.s = f1 if not oFlipped else f2

            # case 21: s & x
            if SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0, 1)
                self.region = "two-phase"

                def fn_sx(p_):
                    sf = self.steamTable.sL_p(p_)
                    sg = self.steamTable.sV_p(p_)
                    return self.s - (sf + self.x * (sg - sf))

                self.p = fsolve(fn_sx, [1.0])[0]
                self.t = self.steamTable.tsat_p(self.p)

        # After setting self.p, self.t, self.region, and possibly self.x,
        # calculate all properties from the steam table
        self.computeProperties()

    def between(self, x, low, high):
        return low <= x <= high

    def clamp(self, x, low, high):
        return max(low, min(x, high))

    def __sub__(self, other):
        """
        Overload the minus operator to get property differences.
        NOTE: Replaced 'other.timeData' with 'other.t' to fix a likely typo.
        """
        delta = thermoState()
        delta.p = self.p - other.p
        delta.t = self.t - other.t
        delta.h = self.h - other.h
        delta.u = self.u - other.u
        delta.s = self.s - other.s
        delta.v = self.v - other.v
        return delta


class main_window(QWidget, Ui__frm_StateCalculator):
    def getShortProp(self, fullText):
        """
        Convert the user's combo-box text (e.g. 'Pressure (P)') into a short property key ('p').
        """
        fullText = fullText.lower()
        if "pressure" in fullText:
            return "p"
        if "temperature" in fullText:
            return "t"
        if "enthalpy" in fullText:
            return "h"
        if "entropy" in fullText:
            return "s"
        if "energy" in fullText:
            return "u"
        if "volume" in fullText:
            return "v"
        if "quality" in fullText:
            return "x"
        return ""

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)  # Default to SI units
        self.currentUnits = 'SI'
        self.setUnits()
        self.SetupSlotsAndSignals()
        self.show()

    def SetupSlotsAndSignals(self):
        """
        Connect GUI buttons and combo box events to their respective callbacks.
        """
        self._pb_Calculate.clicked.connect(self.calculateProperties)
        self._rdo_SI.clicked.connect(self.setUnits)
        self._rdo_English.clicked.connect(self.setUnits)
        self._cmb_Property1.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property2.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property1_2.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property2_2.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property1.currentIndexChanged.connect(self.updateUnitLabels)
        self._cmb_Property2.currentIndexChanged.connect(self.updateUnitLabels)
        self._cmb_Property1_2.currentIndexChanged.connect(self.updateUnitLabels)
        self._cmb_Property2_2.currentIndexChanged.connect(self.updateUnitLabels)

    def updateUnitLabels(self):
        """
        Update the unit labels for each property, depending on user selection (e.g. bar vs psi).
        """
        SI = self._rdo_SI.isChecked()

        def get_units(prop_text):
            if 'Pressure' in prop_text:
                return self.p_Units
            elif 'Temperature' in prop_text:
                return self.t_Units
            elif 'Energy' in prop_text:
                return self.u_Units
            elif 'Enthalpy' in prop_text:
                return self.h_Units
            elif 'Entropy' in prop_text:
                return self.s_Units
            elif 'Volume' in prop_text:
                return self.v_Units
            elif 'Quality' in prop_text:
                return ""
            else:
                return ""

        self._lbl_Property1_Units.setText(get_units(self._cmb_Property1.currentText()))
        self._lbl_Property2_Units.setText(get_units(self._cmb_Property2.currentText()))
        self._lbl_Property1_2_Units.setText(get_units(self._cmb_Property1_2.currentText()))
        self._lbl_Property2_2_Units.setText(get_units(self._cmb_Property2_2.currentText()))

    def setUnits(self):
        """
        Set the display units for all GUI elements based on radio buttons (SI or English).
        Updates internal steam table accordingly and converts displayed numeric values.
        """
        SI = self._rdo_SI.isChecked()
        newUnits = 'SI' if SI else 'EN'
        UnitChange = self.currentUnits != newUnits
        self.currentUnits = newUnits

        if SI:
            self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
            self.p_Units = "bar"
            self.t_Units = "C"
            self.u_Units = "kJ/kg"
            self.h_Units = "kJ/kg"
            self.s_Units = "kJ/kg*C"
            self.v_Units = "m^3/kg"
        else:
            self.steamTable = XSteam(XSteam.UNIT_SYSTEM_FLS)
            self.p_Units = "psi"
            self.t_Units = "F"
            self.u_Units = "btu/lb"
            self.h_Units = "btu/lb"
            self.s_Units = "btu/lb*F"
            self.v_Units = "ft^3/lb"

        if UnitChange:
            # Convert and update for State 1
            self.convertAndUpdateField(self._cmb_Property1, self._le_Property1, SI)
            self.convertAndUpdateField(self._cmb_Property2, self._le_Property2, SI)
            # Convert and update for State 2
            self.convertAndUpdateField(self._cmb_Property1_2, self._le_Property1_2, SI)
            self.convertAndUpdateField(self._cmb_Property2_2, self._le_Property2_2, SI)

        self.updateUnitLabels()

    def convertAndUpdateField(self, combo, lineedit, toSI):
        # Convert the value in the lineedit from EN to SI or vice versa
        try:
            val = float(lineedit.text())
        except ValueError:
            return  # If it's not a valid float, skip
        prop = combo.currentText().lower()

        if 'pressure' in prop:
            val = val * UC.psi_to_bar if toSI else val * UC.bar_to_psi
        elif 'temperature' in prop:
            val = UC.F_to_C(val) if toSI else UC.C_to_F(val)
        elif 'energy' in prop:
            val = val * UC.btuperlb_to_kJperkg if toSI else val * UC.kJperkg_to_btuperlb
        elif 'enthalpy' in prop:
            val = val * UC.btuperlb_to_kJperkg if toSI else val * UC.kJperkg_to_btuperlb
        elif 'entropy' in prop:
            val = val * UC.btuperlbF_to_kJperkgC if toSI else val * UC.kJperkgC_to_btuperlbF
        elif 'volume' in prop:
            val = val * UC.ft3perlb_to_m3perkg if toSI else val * UC.m3perkg_to_ft3perlb
        # Quality (x) is unitless

        lineedit.setText(f"{val:.3f}")

    def calculateProperties(self):
        """
        Read the user inputs for both states, compute their properties,
        display each state, and display the difference.
        """
        SI = self._rdo_SI.isChecked()

        # State 1
        prop1_1 = self.getShortProp(self._cmb_Property1.currentText())
        prop2_1 = self.getShortProp(self._cmb_Property2.currentText())
        val1_1 = float(self._le_Property1.text())
        val2_1 = float(self._le_Property2.text())

        # State 2
        prop1_2 = self.getShortProp(self._cmb_Property1_2.currentText())
        prop2_2 = self.getShortProp(self._cmb_Property2_2.currentText())
        val1_2 = float(self._le_Property1_2.text())
        val2_2 = float(self._le_Property2_2.text())

        # Compute state 1
        try:
            s1 = thermoState()
            s1.setState(prop1_1, prop2_1, val1_1, val2_1, SI)
            self._lbl_StateProperties.setText(self.makeLabel(s1))
        except Exception as e:
            self._lbl_StateProperties.setText("State 1 ERROR: " + str(e))
            return

        # Compute state 2
        try:
            s2 = thermoState()
            s2.setState(prop1_2, prop2_2, val1_2, val2_2, SI)
            self._lbl_StateProperties_2.setText(self.makeLabel(s2))
        except Exception as e:
            self._lbl_StateProperties_2.setText("State 2 ERROR: " + str(e))
            return

        # Display differences
        self._lbl_State.setText("State 1")
        self._lbl_State_2.setText("State 2")
        self._lbl_StateChange.setText("State Change")
        self._lbl_StateChangeProperties_3.setText(self.makeDeltaLabel(s1, s2))

    def computeState(self, prop1, prop2, val1, val2, SI=True):
        """
        A helper method (optional usage) for computing a state if you only had p & t, for example.
        """
        st = thermoState()
        steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)

        def setter(p, t):
            st.p = p
            st.t = t
            st.u = steamTable.u_pt(p, t)
            st.h = steamTable.h_pt(p, t)
            st.s = steamTable.s_pt(p, t)
            st.v = steamTable.v_pt(p, t)
            # If T > Tsat => superheated
            st.x = 1.0 if t > steamTable.tsat_p(p) else 0.0
            st.region = "super-heated" if st.x == 1.0 else "sub-cooled"

        if sorted([prop1, prop2]) == ['p', 't']:
            setter(val1 if prop1 == 'p' else val2, val2 if prop2 == 't' else val1)
        else:
            st.region = "unknown"
        return st

    def makeLabel(self, state):
        """
        Format a multi-line string for the displayed results of a single state.
        """
        stProps = f"Region = {state.region}"
        stProps += f"\nPressure = {state.p:.3f} ({self.p_Units})"
        stProps += f"\nTemperature = {state.t:.3f} ({self.t_Units})"
        stProps += f"\nInternal Energy = {state.u:.3f} ({self.u_Units})"
        stProps += f"\nEnthalpy = {state.h:.3f} ({self.h_Units})"
        stProps += f"\nEntropy = {state.s:.3f} ({self.s_Units})"
        stProps += f"\nSpecific Volume = {state.v:.6f} ({self.v_Units})"
        stProps += f"\nQuality = {state.x:.3f}"
        return stProps

    def makeDeltaLabel(self, s1, s2):
        """
        Format a string showing the change in thermodynamic properties from state 1 to state 2.
        """
        stDelta = "Property change:"
        stDelta += f"\nT2 - T1 = {s2.t - s1.t:.3f} {self.t_Units}"
        stDelta += f"\nP2 - P1 = {s2.p - s1.p:.3f} {self.p_Units}"
        stDelta += f"\nh2 - h1 = {s2.h - s1.h:.3f} {self.h_Units}"
        stDelta += f"\nu2 - u1 = {s2.u - s1.u:.3f} {self.u_Units}"
        stDelta += f"\ns2 - s1 = {s2.s - s1.s:.3f} {self.s_Units}"
        stDelta += f"\nv2 - v1 = {s2.v - s1.v:.6f} {self.v_Units}"
        return stDelta


def main():
    """
    Start the PyQt5 application and show the main window.
    """
    app = QApplication(sys.argv)
    win = main_window()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
