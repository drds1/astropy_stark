import numpy as np
import astropy_stark.mytemp0 as mt0
twopi = np.pi * 2
deg2rad = np.pi / 180
planck = 6.626307004e-34
c = 2.99792458e8
boltz = 1.38064852e-23


from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
from pycallgraph import Config


class pytfbclass:
    def __init__(self,tau,embh, emdot, wavang, deginc,thcent=0,thfwhm=0,
                 t0vin=-1, t0iin=-1, alpha_visc=-0.75, hxsch=3.0, alpha_irad=-0.75,
                 eta=0.1, rlosch=3.0, norm=1, quick=1, xstop=15, udlnr=0.01, oldsmooth=0,
                 newsmooth=1, diagnose=0
                 ):
        self.dtau = tau[1]-tau[0]
        self.taus = tau
        self.wavang = wavang
        self.emdot = emdot
        self.embh = embh
        self.thcent = thcent
        self.thfwhm = thfwhm
        self.Ntau = len(tau)
        self.t0vin = t0vin
        self.t0iin = t0iin
        self.alpha_visc = alpha_visc
        self.alpha_irad = alpha_irad
        self.hxsch = hxsch
        self.eta=eta
        self.rlosch=rlosch
        self.norm = norm
        self.quick = quick
        self.xstop = xstop
        self.udlnr = udlnr
        self.oldsmooth = oldsmooth
        self.newsmooth = newsmooth
        self.diagnose = diagnose
        self.deginc = deginc
        self.psis = np.zeros(self.Ntau)



    def prep(self):
        wavang=self.wavang
        taus = self.taus
        psis = self.psis
        thcent = self.thcent
        thfwhm = self.thfwhm
        embh = self.embh
        emdot = self.emdot
        t0vin = self.t0vin
        t0iin = self.t0iin
        eta = self.eta
        hxsch = self.hxsch
        deginc = self.deginc
        rlosch = self.rlosch
        alpha_visc = self.alpha_visc
        alpha_irad = self.alpha_irad
        xstop = self.xstop
        udlnr = self.udlnr

        # if you want a top hat response then its easy else do what you had before
        if (wavang < 0.0):
            idxinc = np.where((taus > thcent - thfwhm / 2) & (taus < thcent + thfwhm / 2))[0]
            psis[idxinc] = 1


        else:
            # either input desired reference temperature t0 at 1 light day, or calculate the
            # value appropriate for black body accretion disc given m and mdot
            if (t0vin < 0):
                t0v = mt0.tv0(embh, emdot)
            else:
                t0v = t0vin

            if (t0iin < 0):
                t0i = mt0.ti0(embh, emdot, eta=eta)
            else:
                t0i = t0iin

                # need to calculate 1/T^3 (r/r0)^alpha_irad x^2/wav^2 / (cosh(x) - 1) delta (tau - tau(r,theta)) dtau

            t0v2 = t0v * t0v
            t0v4 = t0v2 * t0v2

            t0i2 = t0i * t0i
            t0i4 = t0i2 * t0i2

            rsch = 1.15821e-10 * embh  # scwarzchild radius in light days
            hx = hxsch * rsch
            hx2 = hx * hx
            cosi = np.cos(deginc * deg2rad)
            sini = np.sin(deginc * deg2rad)
            hxci = hx * cosi
            hc_kwav = planck * c / wavang / 1.e-10 / boltz
            wavang2 = wavang * wavang
            rlold = rlosch * rsch

            # this estimate for the cutoff radius is based on the max radius of the
            # highest lag (re-arrange equation 3 in Starkey et al 2017)
            # should be ok but might exclude some significant low lags for a VERY hot, edge on disk

            av4 = alpha_visc * 4
            ai4 = alpha_irad * 4

        taus  = self.taus
        # use a cutoff x_stop to determine when to stop the radius grid
        dtau = taus[1] - taus[0]
        rhilog = dtau
        rtemp = np.array([rhilog, 10 * rhilog])
        rtl = np.log(rtemp)
        ttemp4 = t0v4 * (rtemp) ** av4 + t0i4 * (rtemp) ** ai4
        ttemp = np.sqrt(np.sqrt(ttemp4))
        ttl = np.log(ttemp)
        tstop = hc_kwav / xstop
        grad = (rtl[1] - rtl[0]) / (ttl[1] - ttl[0])
        rhil = rtl[0] + grad * (np.log(tstop) - ttl[0])
        rhild = np.exp(rhil)

        rgridlog = np.exp(
            np.arange(np.log(rlold), np.log(rhilog), udlnr)[:-1])  # np.logspace(np.log(rlold),np.log(rhilog))
        rgridlin = np.arange(rhilog, rhild, rhilog)


        self.rgrid = np.concatenate((rgridlog, rgridlin))
        self.av4 = av4
        self.ai4 = ai4
        self.t0v4 = t0v4
        self.t0i4 = t0i4
        self.drad = dra



    def get_radii(self):
        pass



    def get_azimuths(self):
        # now azimuth grid
        rgrid = self.rgrid
        radlo = rgrid[:-1]
        radhi = rgrid[1:]
        drad = radhi - radlo

        azwidth = drad / radlo
        azgrid = np.arange(0.0, twopi, azwidth)


        X, Y = np.mgrid[-5:5.1:0.5, -5:5.1:0.5]
        a1, a2 = np.mgrid[0:twopi:, -5:5:21j]
        naz = np.shape(azgrid)[0]
        nazsub1 = naz - 1

        raz = np.random.uniform(radlo, radhi, nazsub1)
        daz = np.sqrt(raz * raz + hx2)
        az = np.random.uniform(low=azgrid[:-1], high=azgrid[1:],
                               size=nazsub1)  # np.random.uniform(low=0,high=1,size=nazsub1)*azgrid_s + azgrid[:-1]#np.random.uniform(azgrid[:-1],azgrid[1:],1)
        caz = np.cos(az)
        tdl = hxci - raz * caz * sini + daz

        pass


    def get_temps(self):
        av4 = self.av4
        ai4 = self.ai4
        t0v4 = self.t0v4
        t0i4 = self.t0i4
        rgrid = self.rgrid
        # calculate temperature at each radius grid
        tv4 = t0v4 * (rgrid) ** av4
        ti4 = t0i4 * (rgrid) ** ai4
        rir0b = ti4 / t0i4  # this is just (r/r0)^alpha_irad
        ttot4 = tv4 + ti4
        ttot2 = np.sqrt(ttot4)
        ttot = np.sqrt(ttot2)
        ttot3 = ttot2 * ttot
        self.ttot = ttot
        self.ttot2 = ttot2
        self.ttot3 = ttot3
        self.ttot4 = ttot4







if __name__ == '__main__':

    taugrid = np.arange(0, 30.1, 0.1)
    embh = 1.e7
    emdot = 1.0
    wavnow = 5000
    deginc = 0.0


    config = Config(max_depth=20)
    graphviz = GraphvizOutput(output_file="profile_mytfbquick.png")
    with PyCallGraph(output=graphviz, config=config):
        x = pytfbclass(
            taugrid,
            embh,
            emdot,
            wavnow,
            deginc,
            t0vin=-1,
            t0iin=-1,
            alpha_visc=-0.75,
            hxsch=3.0,
            alpha_irad=-0.75,
            eta=0.1,
            rlosch=3.0,
            norm=1,
            quick=1,
            xstop=15,
            udlnr=0.01,
            thcent=1.0,
            thfwhm=0.2,
            oldsmooth=0,
            newsmooth=1,
            diagnose=0,
        )

        x.prep()

        x.get_temps()

    rgrid = x.rgrid


