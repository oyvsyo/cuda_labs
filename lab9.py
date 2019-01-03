from __future__ import division
import numpy as np
import ROOT
import lab4
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


Na = 6.02e6

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

def length(x1, y1, z1, x2, y2, z2):
    return (x2-x1)**2+(y2-y1)**2+(z2-z1)**2

def calc_p(n0):
    return (np.random.uniform(0., 1., 1)[0] / n0 ** (-2. / 3) / np.pi) ** 0.5


def nap_cos(omega0, theta, phi):
    omega1 = list(range(3))
    mu = np.cos(theta)
    omega1[0] = omega0[0] * mu - (omega0[1] * np.sin(phi) + omega0[0] * omega0[2] * np.cos(phi)) * ((1 - mu ** 2) / (1 - omega0[2] ** 2)) ** .5
    omega1[1] = omega0[1] * mu + (omega0[0] * np.sin(phi) - omega0[1] * omega0[2] * np.cos(phi)) * ((1 - mu ** 2) / (1 - omega0[2] ** 2)) ** .5
    omega1[2] = omega0[2] * mu + (1 - omega0[2] ** 2) * np.cos(phi) * ((1 - mu ** 2) / (1 - omega0[2] ** 2)) ** .5
    return omega1


class Atom(object):
    """Class for Atom"""
    def __init__(self, Z=1., M=2., name="d"):
        self.Z = Z
        self.M = M
        self.name = name


class Particle(Atom):
    """Class for particle"""
    def __init__(self, X=[0, 0, 0], Omega=[1., 0., 0.], E=100., M=2., Z=1., name="-"):
        super(Particle, self).__init__(Z, M, name)
        self.X = X
        self.E = E
        self.Omega = Omega


class Molequla:
    """Class for molequle"""
    def __init__(self, atoms=[Atom()], weights=[]):
        self.atoms = atoms
        self.weights = weights
        self.probability = np.array(weights)/float(sum(weights))
        names = ""
        for i in range(len(atoms)):
            if weights[i] == 1:
                names += atoms[i].name
            else:
                names += "%s_{%d}" % (atoms[i].name, weights[i])
        self.name = names

    def __str__(self):
        return self.name


class Environ:
    """Class for enviroument"""
    def __init__(self, molequla, ro=1.):
        self.molequla = molequla
        self.n0s = list(map(lambda i: ro*Na/molequla.atoms[i].M*molequla.weights[i], range(len(molequla.atoms))))
        self.Ls = list(map(lambda n0: n0 ** (-1. / 3), self.n0s))
        self.ro = ro

    def inter_partner(self):
        return np.random.choice(self.molequla.atoms,
                                1,
                                p=self.molequla.probability)[0]

    def do_interact(self, particle):
        partner = self.inter_partner()
        index = self.molequla.atoms.index(partner)
        n0 = self.n0s[index]
        L = self.Ls[index]
        p = calc_p(n0)
        M1, Z1 = particle.M, particle.Z
        M2, Z2 = partner.M, partner.Z
        E = particle.E
        Emax = 4*M1*M2*E/(M1+M2)
        sin2 = lab4.sintheta(Z1, Z2, M1, M2, E, p, method=0)

        dE = Emax*sin2
        # new energy
        particle.E -= dE

        cost = 1 - 2*sin2
        sint = 1 - cost**2
        theta_lab = np.arctan(sint/(cost+M1/M2))
        phi = np.random.uniform(0, 1, 1)[0]*2*np.pi
        # new Omega
        particle.Omega = nap_cos(particle.Omega, theta_lab, phi)

        particle.X = [
            particle.X[0] + particle.Omega[0] * L,
            particle.X[1] + particle.Omega[1] * L,
            particle.X[2] + particle.Omega[2] * L
        ]
        return index


H = Atom(name="H")
He = Atom(2., 4., name="He")
C = Atom(6., 12.01, name="C")
F = Atom(9., 18., name="F")
O = Atom(8., 16., name="O")
Si = Atom(14., 28.086, name="Si")
Li = Atom(3., 6.94, name="Li")
Fe = Atom(26., 55.84, name="Fe")
Ti = Atom(22., 47.86, name="Ti")
T = Atom(1., 3., name="T")

d = Particle(M=2., Z=1., E=110, name="d")
alpha = Particle(Z=2, M=4., E=100, name="alpha")
pLi = Particle(Z=3., M=6.94, E=390, name="Li")
p = Particle(Z=1., M=1., E=20, name="p")
t = Particle(Z=1., M=3., E=200, name="t")

CH2 = Molequla([C, H], [1, 2])
SiO2 = Molequla([Si, O], [1, 2])
Fe2O3 = Molequla([Fe, O], [2, 3])
LiFe = Molequla([Li, Fe], [1, 1])
C12H18O7 = Molequla([C, H, O], [12, 18, 7])
TiT2 = Molequla([Ti, T], [1, 2])
LiF = Molequla([Li, F], [1, 1])


P = p
Env = Environ(C12H18O7, 1.)

E = P.E
n = 100
natoms = len(Env.molequla.atoms)
vlth = []
lth = []

# X, Y, Z = range(n), range(n), range(n)
# for i in range(n):
#     X[i] = []
#     Y[i] = []
#     Z[i] = []
#     l = 0
#     while P.E > 0:
#         x, y, z = P.X[0], P.X[1], P.X[0]
#         Env.do_interact(P)
#         x1, y1, z1 = P.X[0], P.X[1], P.X[0]
#         l += length(x, y, z, x1, y1, z1)
#         print l
#         X[i] += [P.X[0]]
#         Y[i] += [P.X[1]]
#         Z[i] += [P.X[2]]
#     lth += [l]
#     vlth += [length(X[i][0], Y[i][0], Z[i][0], X[i][-1], Y[i][-1], Z[i][-1])]
#     print l, "   ", vlth[i]
#     P.E = E
#
# h1 = ROOT.TH1F("g", "g", 50, 0., .5)
# map(lambda x: h1.Fill(x), lth)
# h1.Draw()
coordinates = []
for n in range(30):
    # print(f"{n}\r")
    X, Y, Z = list(range(natoms)), list(range(natoms)), list(range(natoms))
    graps = list(range(natoms))
    x, y, z = [], [] ,[]
    while P.E > 0:
        index = Env.do_interact(P)
        if type(X[index])==type(1):
            X[index] = []
        if type(Y[index])==type(1):
            Y[index] = []
        if type(Z[index])==type(1):
            Z[index] = []
        x += [P.X[0]]
        y += [P.X[1]]
        z += [P.X[2]]
        X[index] += [P.X[0]]
        Y[index] += [P.X[1]]
        Z[index] += [P.X[2]]
    p.E = 1+10*n
    p.X = [0, 0, 0]
    p.Omega = [1., 0., 0.]
    try:
        print(f"{max(x)}   {p.E}")

    except:
        pass
    coordinates.append([x, y, z])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


max_max = -999.
min_min = 999.
for x, y, z in coordinates:
    min_ = min([min(x), min(y), min(z)])
    max_ = max([min(x), max(y), max(x)])
    if max_ > max_max:max_max = max_
    if min_ < min_min:min_min = min_

ax.set_xlim(min_min, max_max)
ax.set_ylim(min_min, max_max/5)
ax.set_zlim(min_min, max_max/5)
for x, y, z in coordinates:
    ax.plot(x, y, z)

plt.show()
# gr = ROOT.TGraph2D(len(x), np.array(x), np.array(y), np.array(z))
# gr.SetName("line")
# gr.Draw("LINE")
# gr.GetXaxis().GetXmax()
# gr.SetTitle("GEANT4")
# leg = ROOT.TLegend(0.4, 0.6, 0.8, 0.85)
# leg.SetNColumns(natoms+1)
# leg.AddEntry(gr, P.name+"#Rightarrow"+Env.molequla.name, "s")
# for i in range(natoms):
#     if type(X[i])==type(1):
#         continue
#     graps[i] = ROOT.TGraph2D(len(X[i]), np.array(X[i]), np.array(Y[i]), np.array(Z[i]))
#     graps[i].SetMarkerColor(i+2)
#     graps[i].SetMarkerStyle(i+20)
#     graps[i].SetName(str(i))
#     graps[i].SetMarkerSize(1.0)
#     graps[i].Draw("PSAME")
#     leg.AddEntry(graps[i], Env.molequla.atoms[i].name, "p")
#
# leg.Draw()
# # gr = ROOT.TGraph2D(len(X), np.array(X), np.array(Y), np.array(Z))
# # gr.Draw("LINE")
# ROOT.gPad.Update()
