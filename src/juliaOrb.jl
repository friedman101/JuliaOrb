module juliaOrb

export true2meanAnom, mean2trueAnom, cart2orb, orb2cart, twoBodyProp
export posix2julian, julian2posix, ecef2eciSimple
export earthMu, eartMeanRadius

earthMu = 3.986004418e14
earthMeanRadius = 6371e3

function true2eccAnom(trueAnom, ecc)
    eccAnom = atan2(sqrt(1-ecc^2)*sin(trueAnom), ecc + cos(trueAnom))
    return eccAnom
end

function true2meanAnom(trueAnom, ecc)
    eccAnom = true2eccAnom(trueAnom, ecc)
    meanAnom = eccAnom - ecc*sin(eccAnom)
    return meanAnom
end

function meanAnomResidual(trueAnom, meanAnom, ecc)
    residual = meanAnom - true2meanAnom(trueAnom, ecc)
    return residual
end

function meanAnomResidualWrtTrue(trueAnom, ecc)
    # pardon the copy-pasta from Maxima
    grad = (ecc*sqrt(1-ecc^2)*cos(trueAnom))/sqrt((1-ecc^2)*sin(trueAnom)^2+(cos(trueAnom)+ecc)^2)-(sqrt(1-ecc^2)*sin(trueAnom)^2)/((1-ecc^2)*sin(trueAnom)^2+(cos(trueAnom)+ecc)^2)-
(sqrt(1-ecc^2)*cos(trueAnom)*(cos(trueAnom)+ecc))/((1-ecc^2)*sin(trueAnom)^2+(cos(trueAnom)+ecc)^2)-
(ecc*sqrt(1-ecc^2)*sin(trueAnom)*(2*(1-ecc^2)*cos(trueAnom)*sin(trueAnom)-2*(cos(trueAnom)+ecc)*sin(trueAnom)))/(2*((1-ecc^2)*sin(trueAnom)^2+(cos(trueAnom)+ecc)^2)^(3/2));
    return grad
end

function mean2trueAnom(meanAnom, ecc)
    trueAnom = meanAnom
    tol = 1e-7
    maxIter = 30

    for i = 1:maxIter
        residual = meanAnomResidual(trueAnom, meanAnom, ecc)
        
        if abs(residual) < tol
            return trueAnom
        end

        grad = meanAnomResidualWrtTrue(trueAnom, ecc)
        trueAnom = trueAnom - residual/grad
    end

    assert(false)

end

function twoBodyPropOrb(trueAnom, SMA, ecc, mu, delT)
   meanAnom = true2meanAnom(trueAnom, ecc) 
   n = sqrt(mu/SMA^3)
   meanAnom = meanAnom + n*delT
   trueAnom = mean2trueAnom(meanAnom, ecc)
   return trueAnom
end

function twoBodyProp(pos, vel, mu, delT)
   ecc, inc, RAAN, AOP, trueAnom, SMA = cart2orb(pos, vel, mu)
   trueAnom = twoBodyPropOrb(trueAnom, SMA, ecc, mu, delT)
   pos, vel = orb2cart(ecc, inc, RAAN, AOP, trueAnom, SMA, mu)
   return pos, vel
end

function cart2orb(pos, vel, mu)
    r = pos
    v = vel
    h = cross(r,v)
    n = cross([0, 0, 1], h)
    nMag = norm(n)
    vMag = norm(v)
    rMag = norm(r)
    hMag = norm(h)
    e = 1/mu*( (vMag^2 - mu/rMag)*r - dot(r,v)*v)
    eMag = norm(e)
    zeta = (vMag^2)/2 - mu/rMag

    a = -mu/(2*zeta)
    p = a*(1 - eMag^2)
    i = acos(h[3]/hMag)

    if n[1] != 0
        O = acos(n[1]/nMag)
    else
        O = 0
    end

    if dot(n,e) != 0
        o = acos(dot(n,e)/(nMag*eMag))
    else
        o = 0
    end

    if dot(e,r) != 0
        nu = acos(dot(e,r)/(eMag*rMag))
    else
        nu = 0
    end

    ecc = eMag
    inc = i
    RAAN = O
    AOP = o
    trueAnom = nu
    SMA = p/(1-eMag^2)

    return ecc, inc, RAAN, AOP, trueAnom, SMA
end

function orb2cart(ecc, inc, RAAN, AOP, trueAnom, SMA, mu)
    e = ecc
    i = inc
    O = RAAN
    o = AOP
    nu = trueAnom
    p = SMA*(1 - e^2)

    rPQW = [p*cos(nu)/(1 + e*cos(nu)), p*sin(nu)/(1 + e*cos(nu)), 0]
    vPQW = [-sqrt(mu/p)*sin(nu), sqrt(mu/p)*(e + cos(nu)), 0]

    PQW2IJK = zeros(3,3)
    cO = cos(O)
    sO = sin(O)
    co = cos(o)
    so = sin(o)
    ci = cos(i)
    si = sin(i)

    PQW2IJK[1,1] = cO*co-sO*so*ci
    PQW2IJK[1,2] = -cO*so-sO*co*ci
    PQW2IJK[1,3] = sO*si
    PQW2IJK[2,1] = sO*co+cO*so*ci
    PQW2IJK[2,2] = -sO*so+cO*co*ci
    PQW2IJK[2,3] = -cO*si
    PQW2IJK[3,1] = so*si
    PQW2IJK[3,2] = co*si
    PQW2IJK[3,3] = ci

    r = PQW2IJK*rPQW
    v = PQW2IJK*vPQW

    pos = r
    vel = v

    return pos, vel

end

function ecef2eciSimple(posixTime)
    # http://aa.usno.navy.mil/publications/docs/Circular_179.pdf
    jd = posix2julian(posixTime)
    # UT1 days from J2000
    DU = jd - 2451545.0
    # earth rotation angle
    theta = 0.7790572732640 + 1.00273781191135448*DU
    # julian centuries from J2000
    T = DU/36525
    # Greenwich mean sidereal time in seconds
    GMST = 86400*theta + (0.014506 + 4612.156534*T + 1.3915817*T^2 
    -0.00000044*T^3 - 0.000029956*T^4 - 0.0000000368*T^5)/15
    # Approximate GAST as GMST, ignore equation of equinoxes
    GAST = GMST
    GMSTrad = GMST*2*pi/86400 
    R_e2i = R3(-GMSTrad)
    return R_e2i
end

function posix2julian(posixTime)
    #http://stackoverflow.com/questions/466321/convert-unix-timestamp-to-julian
    jd = posixTime/86400.0 + 2440587.5
    return jd
end

function julian2posix(jd)
    #http://stackoverflow.com/questions/466321/convert-unix-timestamp-to-julian
    posixTime = 86400.0*(jd - 2440587.5)
    return posixTime
end

function R3(theta)
    y = [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta) 0;
        0 0 1]
    return y
end

end
