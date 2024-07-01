# RadonUtils.py
# m.biester

import math
import sys
import numpy as np


def intersections(d, phi_deg, x_l, x_u, y_l, y_u):
    """
    assumption:

    a projection line is defined by parameters [t, s, phi_rad]

    x(t,s,phi_rad) := t * cos(phi_rad) - s * sin(phi_rad)
    y(t,s,phi_rad) := t * sin(phi_rad) + s * cos(phi_rad)

    -------
    t       : length of vector perpendicular to projection line
    s       : shift parameter along projection line
    phi_rad : angle of vector perpendicular vector
    
    x_l     : minimum of x coordinate
    x_u     : maximum of x coordinate
    y_l     : minimum of y coordinate
    y_u     : maximum of y coordinate

    returns:

    
    
    """
    phi_rad = math.radians(phi_deg)
    sin_phi = math.sin(phi_rad)
    cos_phi = math.cos(phi_rad)

    # handling special cases
    # case#1 (vertical projection line)
    if phi_deg == 0:
        if (t >= x_l) and (t <= x_u):
            pCount = 2
            return pCount, [[t, y_l], [t, y_u]]
        else:
            pCount = 0
            return pCount, []
            
    # case#2 (horizontal projection line)  
    if phi_deg == 90:
        if (t >= y_l) and (t <= y_u):
            pCount = 2
            return pCount, [[x_l, t], [x_u, t]]
        else:
            pCount = 0
            return pCount, []
            
    # standard cases 3, 4, 5, 6
    # case 3: intersection for x=x=x_l (left y boundary)
    sc3 = (t*cos_phi - x_l)/sin_phi
    y_sc3 = t*sin_phi + sc3*cos_phi

    # case 4: intersection for x=x=x_u (right y boundary)
    sc4 = (t*cos_phi - x_u)/sin_phi
    y_sc4 = t*sin_phi + sc4*cos_phi

    # case 5: intersection for y=y_u (upper x boundary: top)
    sc5 = (y_u - d*sin_phi)/cos_phi
    x_sc5 = d*cos_phi - sc5*sin_phi

    # case 6: intersection for y=
    sc6 = (y_l - d*sin_phi)/cos_phi
    x_sc6 = d*cos_phi - sc6*sin_phi

    pCount = 0        
    iPoints = []

    # case#3: x=x_l; y_l <= y <= y_u
    if (y_sc3 >= y_l) and (y_sc3 <= y_u):
        iPoints.append([x_l, y_sc3])
        pCount += 1
    # case#4: x=x_u; y_l <= y <= y_u        
    if (y_sc4 >= y_l) and (y_sc4 <= y_u):
        iPoints.append([x_u, y_sc4])
        pCount += 1
    # case#5: y=y_u; x_l <= x <= x_u        
    if (x_sc5 >= x_l) and (x_sc5 <= x_u):
        iPoints.append([x_sc5, y_u])
        pCount += 1 
    # case#6: y=y_l; x_l <= x <= x_u        
    if (x_sc6 >= x_l) and (x_sc6 <= x_u):
        iPoints.append([x_sc6, y_l])
        pCount += 1 
    
    if pCount == 0:
        return pCount, []
        
    # point with smaller x coordinate shall be first item in iPoints (convention)
    if iPoints[0][0] > iPoints[1][0]:
        iPoints.reverse()
    # pCount should be either 2 (intersections)
    # or 0 -> no intersections; other values indicate some error in the implementation.
    
    return pCount, iPoints


def digitalProjectionLine(iPoints, Nx, Ny, x_l, x_u, y_l, y_u, debug=False):

    if len(iPoints ) == 0:
        return None, None

    # get intersection points (x1, y1) , (x2, y2)
    # by convention x1 <= x2
    x1 = iPoints[0][0]
    y1 = iPoints[0][1]
    x2 = iPoints[1][0]
    y2 = iPoints[1][1]

    del_Ax = (x_u - x_l)/Nx
    del_Ay = (y_u - y_l)/Ny
    dx = x2 - x1
    dy = y2 - y1

    # discrete starting and end points
    # point (x1, y1) -> (dp_x1, dp_y1)
    n_temp = int(math.round((x1 - x_l)/del_Ax))
    n_x1 = dp_temp if dp_temp < Nx else Nx-1
    n_temp = int(math.round((x2 - x_l)/del_Ax))
    n_x2 = dp_temp if dp_temp < Nx else Nx-1
    del_nx = n_x2 - n_x1
    # point (x2, y2) -> (dp_x2, dp_y2)
    n_temp = int(math.round((y1 - y_l)/del_Ay))
    n_y1 = dp_temp if dp_temp < Ny else Ny-1
    n_temp = int(math.round((y2 - y_l)/del_Ay))
    n_y2 = dp_temp if dp_temp < Ny else Ny-1
    del_ny = n_y2 - n_y1
    
    # slope of discrete projection line
    
    slope_d = ((n_y2 - n_y1) * del_Ay)
    

    if del_nx == 0:
        # vertical projection line (Ny points)
        vec_nx = np.repeat(n_x1, Ny, dtype=np.int32)
        vec_ny = np.arange(Ny, dtype=np.int32)
        return pass
    elif del_ny == 0:
        vec_ny = np.repeat(n_y1, Ny, dtype=np.int32)
        vec_nx = np.arange(Nx, dtype=np.int32)
        return pass
    else:
        # projection line with defined slope
        slope_d = ((n_y2 - n_y1) * del_Ay) / ((n_x2 - n_x1) * del_Ax)
### ab hier weiter machen.
    # case: xStepSize >= yStepSize -> generate x array first
    if abs(slope) <= 1:
        xIndexStart = int(math.floor(x1))
        xIndexStop = int(math.floor(x2))
        yStepSize = slope
        # nx = xIndexStop - xIndexStart + 1
        nx = xIndexStop - xIndexStart
        if debug:
            print(f"nx : {nx}")
        xIndices = np.arange(nx, dtype=np.int32) + xIndexStart
        yIndices = np.int32(np.rint(np.arange(nx) * yStepSize + y1))
        
        
    # case: xStepSize < yStepSize -> generate y array first
    else:
        yIndexStart = int(math.floor(y1))
        yIndexStop = int(math.floor(y2))
        xStepSize = 1.0/abs(slope)
        # ny = abs(yIndexStop - yIndexStart) + 1
        ny = abs(yIndexStop - yIndexStart)
        if debug:
            print(f"ny : {ny}")
        if yIndexStop > yIndexStart:
            yIndices = np.arange(ny, dtype=np.int32) + yIndexStart
        else:
            yIndices = yIndexStart - np.arange(ny, dtype=np.int32) 
            
        xIndices = np.int32(np.rint(np.arange(ny) * xStepSize + x1))

    return xIndices, yIndices


def projectionSingleLine(dVal, phi_deg, imgMat, x_l, x_u, y_l, y_u, Nx, Ny, placeholder=math.nan, debug=True):
    """
    dVal    : a scalar d-value
    phi_deg : angle of d-Vector with respect to x-axis
    imgMat  : matrix of image; nr of rows: Ny ; nr of columns: Nx; imgMat[0, 0] is image point
              at upper left corner; element of imgMat must be scalar (eg. greyscale image)
    x_l     : minimum of x coordinate
    x_u     : maximum of x coordinate
    y_l     : minimum of y coordinate
    y_u     : maximum of y coordinate
    Nx      : nr of pixel in x direction (nr. of columns of image matrix)
    Ny      : nr of pixel in y direction (nr. of rows of image matrix)
    debug   : if True -> additional output is displayed
    
    returns:

    projectionVal : if the projection line crosses image, the corresponding points on the line are summed, the 
                    result is returned
                    if the projection line is outside the image a placeholder value is returned
    """
    # compute intersection
    pCount, iPoints = intersections(dVal, phi_deg, x_l, x_u, y_l, y_u)

    if pCount == 0:
        # projection line outside image
        projectionVal = placeholder    
    elif pCount == 2:
        # if intersection is found -> compute projection line
        xIndices, yIndices = projectionIndices(iPoints, Nx, Ny, debug=debug)
        projectionVal = np.sum(imgMat[Ny - 1 - yIndices, xIndices])
    else:
        # this should never happen (only 0 or 2 intersections are valid results)
        sys.exit(f"error -> len(iPoints) : {len(iPoints)} ; iPoints: {iPoints}")
    
    return projectionVal


def projectionMultiLine(dVec, phi_deg, imgMat, x_l, x_u, y_l, y_u, Nx, Ny, placeholder=math.nan, debug=True):
    """
    dVec    : d-vector
    phi_deg : angle of d-Vector with respect to x-axis
    imgMat  : matrix of image; nr of rows: Ny ; nr of columns: Nx; imgMat[0, 0] is image point
              at upper left corner; element of imgMat must be scalar (eg. greyscale image)
    x_l     : minimum of x coordinate
    x_u     : maximum of x coordinate
    y_l     : minimum of y coordinate
    y_u     : maximum of y coordinate
    Nx      : nr of pixel in x direction (nr. of columns of image matrix)
    Ny      : nr of pixel in y direction (nr. of rows of image matrix)
    debug   : if True -> additional output is displayed
    
    returns:

    projectionVal : if the projection line crosses image, the corresponding points on the line are summed, the 
                    result is returned
                    if the projection line is outside the image a placeholder value is returned
    """
    projectionVec = []

    # for a fixed angle phi_deg iterate over all values of dVec and compute a projection value for 
    # each element of dVec
    for dVal in dVec:
        projectionVal = projectionSingleLine(dVal, phi_deg, imgMat, x_l, x_u, y_l, y_u, Nx, Ny, placeholder=math.nan, debug=debug)
        projectionVec.append(projectionVal)
    return np.array(projectionVec)
        