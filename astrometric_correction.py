"""
solve for the astrometric translation between two images
given two matched source lists, find the relationship between coordinates

include the following effects (in order):
  - rotation about an arbitrary point
  - scaling in both x and y
  - zeropoint offsets in x and y

"""
import numpy as np
from scipy.optimize import minimize

def apply_translation(pos, rangle, rpos=(0,0)):
    '''
    rotate about some arbitrary position

    :param pos: input position, tuple (x,y)
    :param rangle: rotation angle (radians)
    :param rpos: rotation coords, tuple (x,y)
    '''

    R = [ [np.cos(rangle), -np.sin(rangle), 0],
          [np.sin(rangle),  np.cos(rangle), 0],
          [0,             0,                1]]

    # shift w.r.t. origin of rotation, rotate, shift back to original frame
    pos1 = apply_zeropoint(pos, -1. * np.array(rpos))
    pos2 = np.dot(R, list(pos1)+[1])[0:2]
    pos3 = apply_zeropoint(pos2, rpos)
    return pos3

def apply_scaling(pos, c):
    '''
    scale x, y coords by different values

    :param pos: input position, tuple (x,y)
    :param c: scale factors, tuple (c_x, c_y)
    ''' 
    C = np.diag(list(c)+[1])
    pos0 = list(pos)+[1]
    new_pos = np.dot(C, pos0)
    return new_pos[0:2]

def apply_zeropoint(pos, d):
    '''
    shift x, y coords by different values

    :param pos: input position, tuple (x,y)
    :param d: scalar offsets, tuple (d_x, d_y)
    '''
    T = np.diag([1.]*3)
    offset = list(d)+[0]
    T[:,2] += offset

    pos0 = list(pos)+[1]
    new_pos = np.dot(T, pos0)
    return new_pos[0:2]

def shift_coords(pos, rangle, rpos, c, d):
    '''
    apply translation:
      - rotate about arbitrary point
      - scale
      - zeropoint offset

    :param pos: input position, tuple (x,y)
    :param rangle: rotation angle (radians)
    :param rpos: rotation coords, tuple (x,y)
    :param c: scale factors, tuple (c_x, c_y)
    :param d: scalar offsets, tuple (d_x, d_y)
    '''

    pos1 = apply_translation(pos, rangle, rpos)
    pos2 = apply_scaling(pos1, c)
    pos3 = apply_zeropoint(pos2, d)

#    pos1 = apply_zeropoint(pos, d)
#    pos2 = apply_scaling(pos1, c)
#    pos3 = apply_translation(pos2, rangle, rpos)


    return pos3

def solve_offsets(pos1, pos2):
    '''
    solve for the rotation angle, rotation origin, scale factor, and zeropoint
    to translate from pos2 --> pos1
    metric: minimize the sum of square distances

    pos1, pos2: coordinate lists (index-matched)
        shape (n, 2) where n=number of sources

    >>> example usage:
        # plot offsets prior to translation:
        sns.jointplot(pos1-pos2)
        # solve for translation and plot the offsets:
        x = solve_offsets(p1, p2)
        p2_shifted = [ shift_coords(p_i, **x) for p_i in p2 ]
        sns.jointplot(pos1 - p2_shifted)

    '''

    def get_sum_sq_dist(theta, coords):
        rangle, rpos_x, rpos_y, c_x, c_y, d_x, d_y = theta
        rpos = (rpos_x, rpos_y)
        c = (c_x, c_y)
        d = (d_x, d_y)
        pos1, pos2 = coords
        sqdist = []
        for i, (p1, p2) in enumerate(zip(pos1, pos2)):
            p1_pred = shift_coords(p2, rangle, rpos, c, d)
            dist = ((np.array(p1) - np.array(p1_pred))**2.).sum()
            sqdist.append(dist)
        return np.array(sqdist).sum()

    theta0 = [0.]*7
    result = minimize(get_sum_sq_dist, theta0, args=([pos1, pos2]))
    if result.success:
        x=result.x
        solution = {'rangle':x[0], 'rpos':x[1:3], 'c':x[3:5], 'd':x[5:]}
        print(solution)
        return solution
    else:
        raise ValueError(result.message)

