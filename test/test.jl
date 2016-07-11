include("../src/JuliaOrb.jl")

using JuliaOrb

ecc = 0.02
trueAnom = 1.5
@printf("starting ecc: %0.2f\n", ecc)
@printf("starting true anom: %0.2f rad\n", trueAnom)
meanAnom = true2meanAnom(trueAnom, ecc)
@printf("converting to mean anom: %0.2f rad\n", meanAnom)
trueAnom = mean2trueAnom(meanAnom, ecc)
@printf("converting back to true anom: %0.2f rad\n\n", trueAnom)

pos = [-4845, -4014, 2609]*1e3
vel = [1.26, -5.93, -4.83]*1e3
mu = earthMu
@printf("starting pos: [%0.2f %0.2f %0.2f] m\n", pos[1], pos[2], pos[3])
@printf("starting vel: [%0.2f %0.2f %0.2f] m/s\n", vel[1], vel[2], vel[3])
ecc, inc, RAAN, AOP, trueAnom, SMA = cart2orb(pos, vel, mu)
@printf("ecc: %0.2f inc: %0.2f rad RAAN: %0.2f rad\n", ecc, inc, RAAN)
@printf("AOP: %0.2f rad trueAnom: %0.2f rad SMA: %0.2f m\n", AOP, trueAnom, SMA)
pos, vel = orb2cart(ecc, inc, RAAN, AOP, trueAnom, SMA, mu)
@printf("converting back to pos: [%0.2f %0.2f %0.2f] m\n", pos[1], pos[2], pos[3])
@printf("converting back to vel: [%0.2f %0.2f %0.2f] m/s\n", vel[1], vel[2], vel[3])
delT = 1000
pos, vel = twoBodyProp(pos, vel, mu, delT)
@printf("pos after 1000 second 2-body propagation: [%0.2f %0.2f %0.2f] m\n", pos[1], pos[2], pos[3])
@printf("vel after 1000 second 2-body propagation: [%0.2f %0.2f %0.2f] m/s\n", vel[1], vel[2], vel[3])
pos, vel = twoBodyProp(pos, vel, mu, -delT)
@printf("pos after -1000 second 2-body propagation: [%0.2f %0.2f %0.2f] m\n", pos[1], pos[2], pos[3])
@printf("vel after -1000 second 2-body propagation: [%0.2f %0.2f %0.2f] m/s\n\n", vel[1], vel[2], vel[3])

posixTime = 1468202088
R_e2i = ecef2eciSimple(posixTime)
@printf("rotation from ECEF to ECI at posix time %i is:\n", posixTime)
println(R_e2i)
