import math
# import numpy as np

Infinity = 1e50
epsilon = 1e-12

def new_mass(mass, pos, vel):
    return {'mass': mass, 'pos': pos, 'vel': vel}


def show_masses(masses, flag=True):
    momentum = 0
    for m in masses:
        if flag:
            print(f"{m['pos']:.3f} @ {m['vel']:.3f}", end='\t')
        momentum += m['mass'] * m['vel']
    print(f"{momentum:.3f}")


def bounce(mass1, mass2):
    m1, m2 = mass1['mass'], mass2['mass']
    v1, v2 = mass1['vel'], mass2['vel']
    v1_new = (2*m2*v2 + (m1-m2)*v1) / (m1+m2)
    v2_new =  (2*m1*v1 + (m2-m1)*v2) / (m1+m2)
    return (v1_new, v2_new)


def wall_collide(mass):
    if mass['vel'] >= 0:
        return Infinity
    return -mass['pos'] / mass['vel']


def time_to_collision(m1, m2):
    x1, x2 = m1['pos'], m2['pos']
    v1, v2 = m1['vel'], m2['vel']
    delta_v = v1 - v2
    dist = x2 - x1
    if dist == 0 or delta_v * dist < 0 or abs(delta_v) < 1e-8:
        return Infinity # masses will not collide
    return dist / delta_v


def run_bounces(masses) -> dict:
    # masses is list of [mass, position, velocity]
    locations = [[m[1] for m in masses]]
    times = [0]
    masses = [new_mass(m[0], m[1], m[2]) for m in masses]
    N = len(masses)
    # this is the initial energy (doubled), which is constant
    energy = sum([m['vel']**2 * m['mass'] for m in masses])
    norm_coeffs = [math.sqrt(m['mass']/energy) for m in masses]
    collision_count = 0
    total_time = 0
    vel_list = [ [m['vel'] for m in masses] ]
    bounce_data = [ [m['pos'] for m in masses] ]
    while True:
        tmin = wall_collide(masses[0])
        indices = [(0, -1)]
        for i in range(N):
            for j in range(i+1,N):
                t = time_to_collision(masses[i], masses[j])
                if t == Infinity:
                    continue
                elif abs(t-tmin)<epsilon: # simultaneous collisions
                    indices.append((i, j))
                elif t < tmin: # new shortest time
                    tmin = t
                    indices = [(i, j)]
        if tmin == Infinity: # no collisions detected
            break
        collision_count += 1
        total_time += tmin
        bounce_data.append([tmin, [m['vel'] for m in masses]])
        if collision_count == 10000:
            print("no conclusion after 10000 steps")
            collision_count = -1
            break
        # update all positions
        for m in masses:
            m['pos'] += m['vel'] * tmin
        locations.append([ m['pos'] for m in masses] )
        times.append(total_time)
        # Abort if there are any 3-way collisions (including
        # masses 1 and 2 colliding at the wall)
        if len(set([j for k in indices for j in k])) < 2*len(indices):
            print("Oops - 3 way collision! Bailing out")
            collision_count = -1
            break
        # process all collisions
        for i, j in indices:
            m1, m2 = masses[i], masses[j]
            if j == -1: # first mass hits wall ...
                m1['vel'] *= -1
            else:
                m1['vel'], m2['vel'] = bounce(m1, m2)
                if m1['pos'] > m2['pos']:
                    # this can happen because of roundoff errors,
                    # and should be fixed if we simply
                    # reset the positions to match.
                    m1['pos'] = m2['pos']
        vel_list.append([m['vel'] for m in masses])
    # all collisions are done, just need to record the final details
    t = total_time * 0.1
    times.append(total_time + t)
    locations.append([m['pos'] + m['vel']*t for m in masses])
    # TODO: add enough time so that the two most recent objects to
    # collide are slightly separated ... (0.16 is a kludge)
    bounce_data.append([0.16, [m['vel'] for m in masses]])
    vel_list.append([m['vel']for m in masses])
    return {
        "masses": masses,
        "collisions": collision_count,
        "total_time": total_time,
        "vel": vel_list,
        "coeffs": norm_coeffs,
        "locations": locations,
        "times": times,
        "bounces": bounce_data
    }


def find_momentum(masses):
    return sum([m['mass'] * m['vel'] for m in masses])
    # momentum = 0
    # for m in masses:
    #     momentum += m['mass'] * m['vel']
    # return momentum


def get_final_momentum(masses):
    result = run_bounces(masses)
    return find_momentum(result['masses']), result


def mom3(x1=0.5, m2=1, m3=2):
    # For an initial configuration of
    #    mass 1 at x1 (velocity 0)
    #    mass m2 at 1 (vel 0)
    #    mass m3 at 1.1 (vel -1)
    # returns the final momentum as a fraction of the
    # maximum possible momentum for those masses.
    mom, result = get_final_momentum([(1,x1,0), (m2,1,0), (m3, 1.1, -1)])
    if result['collisions'] == -1: # bailed out
        print(x1, m2, m3)
    # Maximum momentum occurs when all three masses end up with equal
    # positive velocity Ve. Because the energy must be the same as the
    # initial configuration, we have Ve^2 ( 1 + m2 + m3 ) = m3,
    # so Ve = sqrt( m3/(1+m2+m3) ).
    return mom / math.sqrt(m3 * (1 + m2 + m3))


def time3(x1, m3):
    # returns the total time until all bounces are complete
    result = run_bounces([(1,x1,0), (1,1,0), (m3, 1.1, -1)])
    # mom, result = get_final_momentum([(1,x1,0), (1,1,0), (m3, 1.1, -1)])
    if result["collisions"] == -1: # bailed out
        print(x1, m3)
    return result["total_time"]


def find_breaks(m2=1, m3=2):
    # Starting with
    #    mass 1 at x1 (velocity 0)
    #    mass m2 at 1 (vel 0)
    #    mass m3 at 1.1 (vel -1)
    # this finds the different values of x1 where the final momentum
    # changes. At the exact value of these points (x=x1), a three-way
    # collision occurs, so the difference between x<x1 and x>x1 results
    # in the rest of the path playing out differently.
    pos_list = [0.00000001] \
        + [i/50+epsilon for i in range(1,50)] \
        + [0.99999999]
    mom_list = [ (x1, mom3(x1, m2, m3)) for x1 in pos_list ]
    breaks = []
    for i in range(len(mom_list)-1):
        x_left, mom_left = mom_list[i]
        x_right, mom_right = mom_list[i+1]
        if mom_left != mom_right:
            # bisection search for the break ...
            while x_right - x_left > 1e-8:
                x_mid = (x_right+x_left) / 2
                mom_mid = mom3(x_mid, m2, m3)
                if abs(mom_mid-mom_left) < 1e-6:
                    x_left = x_mid
                else:
                    x_right = x_mid
            breaks.append(x_left)
    return breaks
 
# miscellaneous testing code
# print(run_bounces([[1, 0.3, 0], [1, 1, 0], [9.1, 1.1, -1]])['bounces'])
# for x1 in [0.3, 0.35, 0.7, 0.8]:
#     print(x1, mom3(x1, 9.1))
#     run_bounces([[1, x1, 0], [1, 1, 0], [9.1, 1.1, -1]], True)
# print(find_breaks(1, 5.90001))
# exit(0)

# print(find_breaks(1, 8.01))
# print(run_bounces([[1,1,0], [1,2,-1]]))
# print([ (x1, mom3(x1, 2, 8.1)) for x1 in [0.1, 0.2, 0.3] ])
