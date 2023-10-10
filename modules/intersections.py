import math
import sys
import numpy as np


def intersections(d, phi_deg, x_l, x_u, y_l, y_u):
    """
    d       : length if d-Vector
    phi_deg : angle of d-Vector with respect to x-axis
    x_l     : minimum of x coordinate
    x_u     : maximum of x coordinate
    y_l     : minimum of y coordinate
    y_u     : maximum of y coordinate
    """
    # phi_deg to rad
    phi_rad = math.radians(phi_deg)
    sin_phi = math.sin(phi_rad)
    cos_phi = math.cos(phi_rad)

    # handling special cases
    # case#5
    if phi_deg == 0:
        if (d >= x_l) and (d <= x_u):
            pCount = 2
            return pCount, [[d, y_l], [d, y_u]]
        else:
            pCount = 0
            return pCount, []
    # case#6  
    if phi_deg == 90:
        if (d >= y_l) and (d <= y_u):
            pCount = 2
            return pCount, [[x_l, d], [x_u, d]]
        else:
            pCount = 0
            return pCount, []
    # standard cases 1 ... 4        
    tc1 = (d*cos_phi - x_l)/sin_phi
    y_tc1 = d*sin_phi + tc1*cos_phi
    tc2 = (d*cos_phi - x_u)/sin_phi
    y_tc2 = d*sin_phi + tc2*cos_phi
    
    tc3 = (y_u - d*sin_phi)/cos_phi
    x_tc3 = d*cos_phi - tc3*sin_phi
    tc4 = (y_l - d*sin_phi)/cos_phi
    x_tc4 = d*cos_phi - tc4*sin_phi

    pCount = 0        
    iPoints = []

    # case#1: x=x_l; y_l <= y <= y_u
    if (y_tc1 >= y_l) and (y_tc1 <= y_u):
        iPoints.append([x_l, y_tc1])
        pCount += 1
    # case#2: x=x_u; y_l <= y <= y_u        
    if (y_tc2 >= y_l) and (y_tc2 <= y_u):
        iPoints.append([x_u, y_tc2])
        pCount += 1
    # case#3: y=y_u; x_l <= x <= x_u        
    if (x_tc3 >= x_l) and (x_tc3 <= x_u):
        iPoints.append([x_tc3, y_u])
        pCount += 1 
    # case#4: y=y_l; x_l <= x <= x_u        
    if (x_tc4 >= x_l) and (x_tc4 <= x_u):
        iPoints.append([x_tc4, y_l])
        pCount += 1 
    
    if pCount == 0:
        return pCount, []
        
    # point with smaller x coordinate shall be first item in iPoints
    if iPoints[0][0] > iPoints[1][0]:
        iPoints.reverse()
    # pCount should be either 2 (intersections)
    # or 0 -> no intersections; other values indicate some error in the implementation.
    
    return pCount, iPoints


def projectionIndices(iPoints, Nx, Ny, debug=False):

    if len(iPoints ) == 0:
        return None, None

    # get intersection points (x1, y1) , (x2, y2)
    # by convention x1 <= x2
    x1 = iPoints[0][0]
    y1 = iPoints[0][1]
    x2 = iPoints[1][0]
    y2 = iPoints[1][1]
    dx = x2 - x1
    dy = y2 - y1

    if dx == 0:
        # vertical projection line (Ny points)
        return int(round(x1)) * np.ones(Ny, dtype=np.int32), np.arange(Ny, dtype=np.int32)
    elif dy == 0:
        # horizontal projection line (Nx points)
        return np.arange(Nx, dtype=np.int32), int(round(y1)) * np.ones(Nx, dtype=np.int32)
    else:
        # projection line with slope
        slope = dy / dx

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
        