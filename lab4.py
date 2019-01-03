import ROOT
import numpy as np


Culon = ROOT.TF1("Culon", "    [0]*[1]*14.4/x*[2] ")
dCulon = ROOT.TF1("dCulon", "  [0]*[1]*14.4/x* (  [2]/[3] - [4]/x )  ")

fEcran = ROOT.TF1("fEcran", "[0]*exp(-[1]*x) + [2]*exp(-[3]*x) + [4]*exp(-[5]*x) + [6]*exp(-[7]*x)")
dfEcran = ROOT.TF1("dfEcran", "-([0]*[1]*exp(-[1]*x) + [2]*[3]*exp(-[3]*x) + [4]*[5]*exp(-[5]*x) + [6]*[7]*exp(-[7]*x) )")

f = ROOT.TF1("f", "1-[0]/[1]-([2]/x)*([2]/x) ")
fstrih = ROOT.TF1("fstrih", "-[0]/[1]+2*[2]*[2]/(x*x*x)")

ecran_params = [0.1818, 3.2, 0.5099, 0.9432, 0.2802, 0.4029, 0.0282, 0.2016]
for iparam in range(8):
    dfEcran.SetParameter(iparam, ecran_params[iparam])
    fEcran.SetParameter(iparam, ecran_params[iparam])


def sintheta(Z1, Z2, M1, M2, E, p, method=0):

    a = 0.88534 * 0.529 / (pow(Z1, 0.23) + pow(Z2, 0.23))
    a = a * (-0.003 * (pow(Z1, 0.5) + pow(Z2, 0.5)) + 0.98)
    if p == 0: p = 0.00001
    B = p / a
    Ec = E / (1 + M1 / M2)
    rm = a

    Culon.SetParameter(0, Z1)
    Culon.SetParameter(1, Z2)
    Culon.SetParameter(2, fEcran(rm / a))

    dCulon.SetParameter(0, Z1)
    dCulon.SetParameter(1, Z2)
    dCulon.SetParameter(3, a)

    f.SetParameter(1, Ec)
    f.SetParameter(2, p)
    f.SetParameter(0, Culon(rm))

    fstrih.SetParameter(1, Ec)
    fstrih.SetParameter(2, p)

    # for gggg in range(10):
    while(abs(f(rm))>0.0001):
        Culon.SetParameter(2, fEcran(rm/a))
        dCulon.SetParameter(2, dfEcran(rm/a))
        dCulon.SetParameter(4, fEcran(rm/a))
        f.SetParameter(0, Culon(rm))
        fstrih.SetParameter(0, dCulon(rm))
        rm = rm-f(rm)/fstrih(rm)

    ro = -2 * (Ec - Culon(rm)) / dCulon(rm)
    Rm = rm / a
    Rc = ro / a
    Esm = M2 / (M1 + M2) * E
    epsilon = a * Esm / (Z1 * Z2 * 14.4)
    alfa = 1 + 0.75984 * pow(epsilon, -0.5)
    beta = (5.71974 + pow(epsilon, 0.5)) / (6.14171 + pow(epsilon, 0.5))
    gamma = (9.5217 + epsilon) / (6.2612 + epsilon)
    A = 2 * alfa * epsilon * pow(B, beta)
    G = gamma / (pow((1 + A * A), 0.5) - A)

    delta = A * (Rm - B) / (1 + G)
    if method == 0:  # coulomb
        return 1/(1+2*epsilon*B*2*epsilon*B)

    # semi-empirical
    if method == 1:
        if (B + Rc + delta) / (Rm + Rc) < 1 and (B + Rc + delta) / (Rm + Rc) > -1:
            theta = 2 * np.arccos((B + Rc + delta) / (Rm + Rc))
        else:
            # print B
            return 1/(1+2*epsilon*B*2*epsilon*B)
        return np.sin(theta/2)*np.sin(theta/2)

    if B <= 1:
        theta = (1.012 / B) * np.exp(-0.356 * B) / epsilon
    if B > 1:
        theta = 0.707 * pow(B, -1.78) / epsilon

    # impulse
    if method == 2:
        return np.sin(theta / 2)**2


if __name__ == "__main__":
    Z1 = 1.
    Z2 = 3.
    M1 = 3.016
    M2 = 6.015

    E = [1, 10, 100, 1000, 10000]
    n = 100
    Bmax = 2. - 0.001

    X = np.zeros(n)
    Y = np.zeros(n)
    impY = np.zeros(n)
    couY = np.zeros(n)

    a = 0.88534 * 0.529 / (pow(Z1, 0.23) + pow(Z2, 0.23))
    a = a * (-0.003 * (pow(Z1, 0.5) + pow(Z2, 0.5)) + 0.98)

    text_C = "Z_{1}=%1i, M_{1}=%1.1f amu, Z_{2}=%1i, M_{2}=%1.1f amu" % (Z1, M1, Z2, M2)
    c1 = ROOT.TCanvas("c1", "c", 1200, 1000)


    ROOT.gStyle.SetOptStat(0)
    pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.5, 0.33, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0.33, 0.5, 0.66, 1.0)
    pad3 = ROOT.TPad("pad3", "pad3", 0.66, 0.5, 1.0, 1.0)
    pad4 = ROOT.TPad("pad4", "pad4", 0.0, 0.0, 0.33, 0.5)
    pad5 = ROOT.TPad("pad5", "pad5", 0.33, 0.0, 1.0, 0.5)
    pads = [pad1, pad2, pad3, pad4, pad5]
    map(lambda pad: pad.Draw(), pads)

    gr = range(5)
    couGr = range(5)
    impGr = range(5)

    for iE in range(5):
        pads[iE].cd()

        step = Bmax*a/n
        p = 0+0.001

        for j in range(n):
            p=step*j
            B=p/a+0.001
            X[j] = B
            Y[j] = sintheta(Z1, Z2, M1, M2, E[iE], p, method=1)
            impY[j] = sintheta(Z1, Z2, M1, M2, E[iE], p, method=2)
            couY[j] = sintheta(Z1, Z2, M1, M2, E[iE], p, method=0)

        text_A = "E=%1.0f keV" % E[iE]

        if Z1 == 1 and Z2 == 22: text_B = "D #RightarrowTi "
        if Z1 == 2 and Z2 == 6: text_B = "{}^{4}He #Rightarrow{}^{12}C "
        if Z1 == 3 and Z2 == 14: text_B = "{}^{7}Li #RightarrowSi"
        if Z1 == 1 and Z2 == 8: text_B = "p #Rightarrow{}^{16}O"
        if Z1 == 6 and Z2 == 18: text_B = "{}^{12}C #Rightarrow{}^{41}Ar"
        if Z1 == 6 and Z2 == 18: text_B = "{}^{12}C #Rightarrow{}^{41}Ar"
        if Z1 == 1 and Z2 == 3: text_B = "T #Rightarrow{}^{6}Li"

        gr[iE] = ROOT.TGraph(j-1, X, Y)
        gr[iE].SetTitle(text_A)
        gr[iE].GetXaxis().SetTitle("B")
        gr[iE].GetYaxis().SetTitle("sin^{2}(#theta/2)")
        gr[iE].SetMarkerStyle(8)
        gr[iE].SetMarkerSize(0.5)
        gr[iE].Draw("AP")

        impGr[iE] = ROOT.TGraph(j-1, X, impY)
        impGr[iE].SetMarkerColor(2)
        impGr[iE].SetMarkerColor(4)
        impGr[iE].SetMarkerStyle(8)
        impGr[iE].SetMarkerSize(0.3)
        impGr[iE].Draw("P")

        couGr[iE] = ROOT.TGraph(j-1, X, couY)
        couGr[iE].SetMarkerColor(3)
        couGr[iE].SetLineColor(3)
        couGr[iE].SetMarkerColor(3)
        couGr[iE].SetMarkerStyle(8)
        couGr[iE].SetMarkerSize(0.3)
        couGr[iE].Draw("P")

    leg = ROOT.TLegend(0.4, 0.4, 0.8, 0.85)
    leg.AddEntry(couGr[iE], text_B, "s")
    leg.AddEntry(couGr[iE], text_C, "s")
    leg.AddEntry(gr[iE], "Semi-empirical formula", "p")
    leg.AddEntry(impGr[iE], "Impulse approximation", "p")
    leg.AddEntry(couGr[iE], "Coulomb scattering", "p")
    leg.Draw()

    c1.Update()
