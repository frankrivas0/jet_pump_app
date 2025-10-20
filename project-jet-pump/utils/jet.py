import math
import numpy as np
from pandas import DataFrame
import pandas as pd

class ipr:
    def __init__(self, q, p, pr, pb):
        self.q = q
        self.p = p
        self.pr = pr
        self.pb = pb
        if self.pr>self.pb:
            self.j = ipr.productivity_index(self)
            self.qob = ipr.bp_oil_rate(self)

    def ipr_voguel_sat(self, pwf):
        qmax = self.q/(1-0.2*(self.p/self.pr)-0.8*(self.p/self.pr)**2)
        y = lambda x: qmax*(1-0.2*(x/self.pr)-0.8*(x/self.pr)**2)
        qo = y(pwf)
        return qo

    def productivity_index(self):
        if self.p>=self.pb:
            j = self.q/(self.pr-self.p)
        else:
            j = self.q/((self.pr-self.pb)+self.pb*(1-0.2*(self.p/self.pb)-0.8*(self.p/self.pb)**2)/1.8)
        return j

    def bp_oil_rate(self):
        qob = self.j*(self.pr-self.pb)
        return qob

    def ipr_voguel_und(self, pwf):
        if pwf<=self.pb:
            y = lambda x: self.qob+self.j*self.pb*(1-0.2*(x/self.pb)-0.8*(x/self.pb)**2)/1.8
        else:
            y = lambda x: self.j*(self.pr-x)
        qo = y(pwf)
        return qo

    def voguel(self, pwf):
        if self.pr<=self.pb:
            qo = ipr.ipr_voguel_sat(self, pwf)
        else:
            qo = ipr.ipr_voguel_und(self, pwf)
        return qo
    
class pvt_correlations:

    def __init__(self, api, pb, yg):
        self.api = api
        self.pb = pb
        self.yg = yg

    '''Solution Gas-Oil Ratio: Vasquez and Beggs (1980) correlation
    returns in scf/STB'''
    def rso(self, p, t):
        ygc = self.yg * (1 + 5.912 * 10 ** -5 * self.api * 60 * math.log10(25 / 114.7) / math.log(10))
        if self.api<= 30:
            c1 = 0.0362
            c2 = 1.0937
            c3 = 25.7240
        else:
            c1 = 0.0178
            c2 = 1.1870
            c3 = 23.9310
        if p<self.pb:
            rs = c1*ygc*(p**c2)*math.e**(c3*(self.api)/(t+460))
        else:
            rs = c1*ygc*(self.pb**c2)*math.e**(c3*(self.api)/(t+460))
        return rs

    '''Oil Formation Volume Factor: Vasquez and Beggs (1980) correlation
    returns in bbl/STB'''
    def bod(self, p, t):
        rs = pvt_correlations.rso(self, p, t)
        if self.api<= 30:
            c1 = 4.677*10**-4
            c2 = 1.751*10**-5
            c3 = -1.811*10**-8
        else:
            c1 = 4.670*10**-4
            c2 = 1.100*10**-5
            c3 = 1.377*10**-9
        if p<self.pb:
            bo = 1+c1*rs+c2*(t-60)*(self.api/self.yg)+c3*rs*(t-60)*(self.api/self.yg)
        else:
            Bob = 1+c1*rs+c2*(t-60)*(self.api/self.yg)+c3*rs*(t-60)*(self.api/self.yg)
            c0 = (-1433+5*rs+17.2*t-1180*self.yg+12.61*self.api)/(p*10**5)
            bo = Bob*math.e**(c0*(self.pb-p))
        return bo

    ''' Oil Viscosity: Beggs and Robinson (1975) correlation
    returns in centipoise cP'''
    def oviscosity(self, p, t):
        rs = pvt_correlations.rso(self, p, t)
        z = 3.0324-(0.02023*self.api)
        y = 10**z
        x = y*(t**(-1.163))
        uod = 10**x -1
        a=10.715*(rs+100)**(-0.515)
        b=5.44*(rs+150)**-0.338
        ub=a*(uod**b)
        if p<self.pb:
            u = ub
        else:
            x = 2.6*(p**1.187)*10**(-5-p*3.9*10**-5)
            u = ub*(p/self.pb)**x
        return u

    ''' Water Viscosity: Brill and Beggs (1978)
    returns in centipoise cP'''
    def wviscosity(self, t):
        wv = math.e**(1.003-(1.479*10**(-2))*(t)+(1.982*10**-5)*(t)**(2))
        return wv

    ''' Gas Viscosity: Lee, Gonzales and Eakin
    returns in centipoise cP'''
    def gas_viscosity(self, t, pg):
        M = self.yg*28.97
        X = 3.5 + 986/(t+460) + 0.01*M
        Y = 2.4 - 0.2*X
        K = ((9.4+0.02*M)*(t+460)**1.5)/(209+19*M+t+460)
        ug = K*(math.e**(X*(pg/62.4)**Y))/10000
        return ug

    ''' Natural Gas Compressibility: '''
    def compresibility(self, p, t):
        Ppc=756.8-131.07*(self.yg)-3.6*(self.yg)**2
        Tpc=169.2+349.5*(self.yg)-74*(self.yg)**2
        Pr=p/Ppc
        Tr=(t+460)/Tpc
        A = -0.101-(0.36*Tr)+1.3868*(Tr-0.919)**0.5
        B = 0.021 + (0.04275/(Tr-0.65))
        D = 0.6222-0.224*Tr
        E = (0.0657/(Tr-0.86)) - 0.037
        F = 0.32*(math.e**(-19.53*(Tr-1)))
        H = 0.122*(math.e**(-11.3*(Tr-1)))
        C = Pr*(D+E*Pr+F*(Pr**4))
        z = A+B*Pr +(1-A)*(math.e**(-C)) -H*(Pr/10)**4
        return z

    ''' Gas/Oil interfacial tension dyn/cm '''
    def o_tens(self, p, t):
        sig_68 = 39 - 0.2571*self.api
        sig_100 = 37.5 - 0.2571*self.api
        if t<=68:
            m = sig_68
        elif t>=100:
            m = sig_100
        else:
            m = sig_68 + (t-68)*(sig_100-sig_68)/(100-68)
        sigma_o = m*(1-0.024*p**0.45)
        if sigma_o<1:
            sigma_o = 1
        return sigma_o

    ''' Gas/Water interfacial tension dyn/cm'''
    def w_tens(self, p, t):
        sig_74 = 75 - 1.108*p**0.349
        sig_280 = 53 - 0.1048*p**0.637
        if t<=74:
            sigma_w = sig_74
        elif t>=280:
            sigma_w = sig_280
        else:
            sigma_w = (sig_74 + (t-74)*(sig_280-sig_74))/(280-74)
        return sigma_w
    
class gph():
    def __init__(self):
        pass

    def cnl(nl):
        #from math import log10
        x = 0
        c = [-1.93565953153, -0.149078681496, -1.87907830349, -8.16189071362,
             -18.134934767, -21.0627105031, -14.1647815099, -5.76192950961,
             -1.40168276048, -0.187791683544, -0.0106535356774]
        for i in range(11):
            x = x + c[i]*(math.log10(nl))**(i)
        x = 10**x
        return x

    def hly(ngv, nlv, cnl, nd, p):
        y = (nlv/(ngv**0.575))*((p/14.7)**0.1)*(cnl/nd)
        x = ((0.0047+1123.32*y+729489.64*y**2)/(1+1097.1566*y+722153.97*y**2))**(1/2)
        return x

    def sigma(ngv, nl, nd):
        y = (ngv*(nl**0.38))/(nd**2.14)
        if y<=0.025:
            x = 27170*y**3 - 317.52*y**2 + 0.5472*y + 0.9999
        elif y>0.025 and y<=0.055:
            x = -533.33*y**2 + 58.524*y + 0.1171
        else:
            x = 2.5714*y + 1.5962
        if x<=1:
            x = 1
        return x

def pressure_drop(delta_p, p0, pb, yg, yw, t0, t1, api, qprod, qiny, bsw, GOR, dcsg, dtgn, delta_h):
    # Determine average pressure and temperature in the interval, psi, °F
    pavg = p0 + delta_p/2
    tavg = (t0 + t1)/2
    # Determine total flow rate, BFPD
    q = qprod + qiny
    # Determine corrected bsw to account for inyection fluid. Assuming bsw of inyected fluid is cero
    bsw2 = bsw*(qprod)/(qiny+qprod)
    # Determine oil specific gravity
    yo = 141.5/(131.5+api)
    pvt = pvt_correlations(api, pb, yg)
    # Determine gas in solution, scf/STB
    rsavg = pvt.rso(pavg, tavg)
    # Determine oil volumetric factor, bbl/STB
    boavg = pvt.bod(pavg, tavg)
    # Determine oil viscosity, cp
    uoavg = pvt.oviscosity(pavg, tavg)
    # Determine water viscosity, cp
    uwavg = pvt.wviscosity(tavg)
    # Determine liquid viscosity, cp
    ul = uoavg*(1-bsw2)+uwavg*(bsw2)
    # Determine liquid density, lb/cu-ft
    pl = ((yo*62.4+rsavg*yg*0.0764/5.614)/boavg)*(1-bsw2)+(yw*62.4*(bsw2))
    # Determine gas compressibility factor
    z = pvt.compresibility(pavg, tavg)
    # Determine gas density, lb/cu-ft
    pg = yg*0.0764*(pavg/14.7)*(520/(tavg+460))*(1/z)
    # Determine flow area; anular area between casing and production pipe (tubing), sq-ft
    area = math.pi*(dcsg**2-dtgn**2)/(144*4)
    # Determine superficial velocity of liquid, ft/s
    vsl = 5.61*q*(boavg*(1-bsw2)+(bsw2))/(86400*area)
    # Determine superficial velocity of gas, ft/s
    vsg = (qprod*(1-bsw)*(GOR - rsavg))*(14.7/pavg)*((tavg+460)/520)*(z)/(86400*area)
    # If superficial gas velocity is negative, set it to be zero
    if vsg<0:
        vsg = 0
    # Determine mixture velocity, ft/s
    vm = vsl + vsg
    # Determine adimensional numbers Lb and lambda_g
    Lb = 1.071 - (0.2218*(vm)**2)*12/(dcsg-dtgn)
    lambda_g = vsg/vm
    # Set Lb=0.13 if <0.13
    if Lb<0.13:
        Lb = 0.13
    # Determine hydraulic diameter, in
    dh = dcsg - dtgn
    # Determine if the flow is bubble flow
    if lambda_g < Lb:
        # Griffith correlation
        gradient = griffith_corr(vm, vsg, vsl, dh, ul, pl, pg)
    else:
        # Determine gas viscosity
        ug = pvt.gas_viscosity(tavg, pg)
        # Determine Gas/Oil interfacial tension
        surteno = pvt.o_tens(pavg, tavg)
        # Determine Gas/Water interfacial tension
        surtenw = pvt.w_tens(pavg, tavg)
        # Hedgedorn and brown correlation
        gradient = hagedorn_brown_corr(vsl, vsg, vm, pl, surteno, surtenw, dh, ul, ug, pg, bsw2, pavg)
    # Determine total pressure change due to gravity and friction over the interval
    delta_p_new = gradient*delta_h
    return delta_p_new

def griffith_corr(vm, vsg, vsl, dh, ul, pl, pg):
    # Determine liquid hold up with slipping between phases
    Hl = 1-0.5*(1+vm/0.8 -((1+vm/0.8)**2 - 4*vsg/0.8)**(1/2))
    # Determine real liquid velocity
    vl = vsl/Hl
    # Determine Reynolds number
    Nre = (124*dh*pl*vl)/(ul)
    # Determine Moody friction factor
    f = 0.0055*(1+(20000*0.0018/dh+1000000/Nre)**(1/3))
    # Determine mixture density
    pm = pl*Hl + pg*(1-Hl)
    # Determine pressure gradient due to gravity and fricction loss
    gradient = (pm + (12*f*pl*vl**2)/(64.4*dh))/144
    return gradient

def hagedorn_brown_corr(vsl, vsg, vm, pl, surteno, surtenw, dh, ul, ug, pg, bsw2, pavg):
    # Determine no-slip liquid hold up
    lambda_l = vsl/vm
    # Determine no-slip mixture density
    pn = pl*lambda_l+pg*(1-lambda_l)
    # Determine surface tension of the liquid
    surtenl = surteno*(1-bsw2)+surtenw*(bsw2) #dynes/cm
    # Determine a-dimensional numbers
    nl = 0.15726*ul*(1/(pl*surtenl**3))**(1/4)
    cnlc = gph.cnl(nl)
    nlv = 1.938*vsl*(pl/surtenl)**(1/4)
    ngv = 1.938*vsg*(pl/surtenl)**(1/4)
    nd = 120.872*((dh)/12)*(pl/surtenl)**(1/2)
    # Determine a-dimensional numbers form correlation
    hlyc = gph.hly(ngv, nlv, cnlc, nd, pavg)
    y = gph.sigma(ngv, nl, nd)
    # Determine slip liquid hold up
    Hl = hlyc*y
    # Choose highest liquid hold up
    if Hl<lambda_l:
        Hl = lambda_l
    # Determine mixture viscosity
    u_m = (ul**(Hl))*(ug**(1-Hl))
    # Determine slip mixture density
    pm = pl*Hl+pg*(1-Hl)
    # Determine Reynolds number
    Nre = (124*dh*pn*vm)/(u_m)
    # Determine Moody friction factor
    f = 0.0055*(1+(20000*0.0018/dh+1000000/Nre)**(1/3))
    # Determine pf
    pf = (pn**2)/pm
    # Determine pressure gradient due to gravity and fricction loss
    gradient = (pm + (12*f*pf*vm**2)/(64.4*dh))/144
    return gradient

def pressure_profile(t_bottom, t_wellhead, p_wellhead, total_depth, pb,
                     qprod, qiny, bsw, GOR, dcsg, dtgn, api, yg, yw):
    '''
    Function to determine vertical lift performance VLP curve using Hagedorn and Brown correlation and Griffith correlation
    when bubble flow occurs. Assumes vertical well and known wellhead pressure. It was written to develop a VLP curve for a well
    producing with hydraulic jet pumping lift system.
    t_bottom:           Bottomhole temperature (°F)
    t_wellhead:         Wellhead temperature (°F)
    p_wellhead:         Wellhead pressure (psi)
    total_depth:        Bottomhole depth (ft)
    pb:                 Oil bubble point pressure (psi)
    qprod:              Oil production rate (stb/d)
    qiny:               Oil inyection rate (stb/d)
    bsw:                Production water cut
    GOR:                Production gas oil ratio (scf/stb)
    dcsg:               Casing internal diameter (in)
    dtgn:               Tubing external diameter (in)
    api:                Production Oil api gravity
    yg:                 Gas specific gravity
    yw:                 Water specific gravity
    surteno:            Gas/Oil superficial tension (dyn/cm)
    Surtenw:            Gas/Water superficial tension (dyn/cm)
    '''
    # set number of steps
    steps = 10
    # determine depth intervals
    delta_h = total_depth/steps
    H = np.linspace(0, total_depth, steps+1)
    T = np.linspace(t_wellhead, t_bottom, steps+1)
    P = [p_wellhead]
    stop_rel_error = 0.01
    for i in range(steps):
        p0 = P[i]
        t0 = T[i]
        t1 = T[i+1]
        delta_p = 25
        error = 100
        while error>stop_rel_error:
            delta_p_new = pressure_drop(delta_p, p0, pb, yg, yw, t0, t1, api, qprod, qiny, bsw, GOR, dcsg, dtgn, delta_h)
            error = abs(delta_p - delta_p_new)/delta_p
            delta_p = delta_p_new
        P.append(p0+delta_p)

    table = DataFrame({'Depth':H, 'Temperature': T, 'Pressure':P})
    return table

def pressure_gradient(yo, q, u, d, e=0.0018, flow_direction=-1):
    '''
    api:               Oil gravity
    q:                 Oil flow rate bbl/d
    u:                 Oil viscosity cp
    e:                 Tubing roughness in
    d:                 Tubing inside diameter in
    flow_direction:    Whether flow is upward (+1) or downward (-1)
    '''
    # Determine fluid velocity
    v = 0.01191386*q/d**2 # ft/s
    # Determine fluid density
    pl = 62.4*yo # lbm/cu-ft
    # Determine Reynolds number
    Nre = (124*d*pl*v)/(u)
    # Determine friction factor according to flow regime
    if Nre<2000:
        f = 64/Nre
    else:
        f = 0.0055*(1+(20000*e/d+1000000/Nre)**(1/3))
    # Determine pressure gradient
    grad = (pl + flow_direction*(12*f*pl*v**2)/(64.4*d))/144   # psi/ft
    return grad

def nozzle_q_p(p_inj, pwf, mu_inj, yo_inj, d_tbgID, total_depth, aj, e=0.0018, flow_direction=-1):
    '''
    Function to calculate the pressure of the power fluid at the nozzle and its flow rate. They are
    calculated through an iterative process until convergence.
    p_inj                Power fluid injection pressure psi
    pwf                  Bottomhole flowing pressure psi
    mu_inj               Power fluid viscosity cp
    yo_inj               Power fluid specific gravity
    d_tbgID              Tubing ID in
    total_depth          Bottomhole depth ft
    aj                   Nozzle area in2
    e                    Tubing roughness in
    flow_direction       Whether flow is upward (+1) or downward (-1)
    '''
    # Calculate an initial pressure at the nozzle
    p_nozzle = p_inj + (62.4/144)*yo_inj*total_depth
    # Define a desired relative error
    stop_rel_error = 0.001
    error = 1
    while error>stop_rel_error:
        # Determine power fluid rate
        q_inj = 832*(aj)*(((p_nozzle-pwf)/(62.4*yo_inj/144))**(1/2))
        # Determine single phase pressure gradient
        gradient = pressure_gradient(yo_inj, q_inj, mu_inj, d_tbgID, e=0.0018, flow_direction=-1)
        # Calculate the pressure at the nozzle
        p_nozzle_new = p_inj + gradient*total_depth
        error = abs(p_nozzle - p_nozzle_new)/p_nozzle
        # Update pressure nozzle
        p_nozzle = p_nozzle_new
    # Update power fluid rate
    q_inj = 832*(aj)*(((p_nozzle-pwf)/(62.4*yo_inj/144))**(1/2))
    return q_inj, p_nozzle

def M_graph(r, H, kn=0.03, ktd=0.20):
    '''
    r      Nozzle and throat area ratio
    H      Pressure ratio
    kn     Nozzle loss coeficient
    ktd    Throat loss coeficient
    '''
    a2 = 2*r
    b2 = (1-2*r)*(r**2)/(1-r)**2
    c2 = (1+ktd)*r**2
    dd2 = 1+kn
    M = (2*c2-((-2*c2)**2-4*(b2-c2)*(a2-c2-H*dd2/(H+1)))**0.5)/(2*(b2-c2))
    return M

def jet_pump(pwf, q_prod, aj, at, p_inj, p_wellhead, pb, t_bottom, t_wellhead, total_depth,
            dcsg, dtbgID, dtbgOD, api, GOR, bsw, yg, yw, mu_inj):
    '''
    Function to determine the operational parameters of a jet pump. It assumes that the api gravity
    is the same for the production and injection fluid. Also assumes the bsw and GOR of injection
    fluid to be negligible.
    pwf                Bottomhole flowing pressure psi
    q_prod             Desired production fluid flow bbl/d
    aj                 Nozzle area sq-in
    at                 Throat area sq-in
    p_inj              Power fluid injection pressure psi
    p_wellhead         Wellhead pressure of the returned fluid psi
    pb                 Bubble point pressure of the production fluid psi
    t_bottom           Bottomhole temperature F
    t_wellhead         Wellhead temperature F
    total_depth        Bottomhole depth ft
    dcsg               Casing ID in
    dtbgID             Tubing ID in
    dtbgOD             Tubing OD in
    api                Production fluid api gravity
    GOR                Production fluid gas oil ratio scf/stb
    bsw                Production fluid bsw
    yg                 Gas gravity lb/cu-ft
    yw                 Water gravity lb/cu-ft
    mu_inj             Power fluid viscosity cp
    '''
    # Calculate nozzle throat area ratio
    r = aj/at
    # Calculate power fluid gravity
    yo_inj = 141.5/(api+131.5)
    # Calculate production fluid gravity
    yo = 141.5/(api+131.5)
    # Calculate production fluid pressure gradient
    suc_gradient = 62.4*(yo*(1-bsw) + yw*bsw)/144
    # Calculate power fluid rate and pressure at the nozzle
    q_inj, p_nozzle = nozzle_q_p(p_inj, pwf, mu_inj, yo_inj, dtbgID, total_depth, aj, e=0.0018, flow_direction=-1)
    stop_rel_error = 0.01
    error = 1
    while error>stop_rel_error:
        # Calculate pump discharge pressure using Hagedorn and Brown multiphase correlation
        VLP_data = pressure_profile(t_bottom, t_wellhead, p_wellhead, total_depth, pb,
                     q_prod, q_inj, bsw, GOR, dcsg, dtbgOD, api, yg, yw)
        p_discharge = VLP_data.iloc[-1]['Pressure']
        # Calculate pressure ratio
        H = (p_discharge-pwf)/(p_nozzle-p_discharge)
        # Calculate flow ratio using graph
        Mgraph = M_graph(r, H)
        # Calculate actual flow ratio
        Mequation = q_prod*((1+2.8*((GOR/pwf)**1.2))*(1-bsw)+bsw)*suc_gradient/(q_inj*0.433*yo_inj)
        # Calculate new production fluid rate
        q_prod_new = q_prod*Mgraph/Mequation
        error = abs(q_prod - q_prod_new)/q_prod
        # Update production fluid rate
        q_prod = q_prod_new
    # Calculate power requirement (does not include pump efficiency)
    power = q_inj*p_inj*0.000017
    # Calculate jet pump efficiency
    efficiency = Mgraph*H*100
    # Calculate minimun area to avoid cavitation
    acm = q_prod*((1/691)*(suc_gradient/pwf)**(1/2)+((1-bsw)*GOR)/(24650*pwf))
    # Calculate maximun production fluid rate to avoid cavitation
    qcav = q_prod*(at-aj)/acm
    return q_prod, pwf, q_inj, p_discharge, p_nozzle, power, efficiency, acm, qcav

