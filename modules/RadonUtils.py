# RadonUtils.py
# m.biester

import math
import sys
import numpy as np


def intersections(t, phi_deg, x_l, x_u, y_l, y_u):
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

    pCount : number of intersections (expected values are: 0 -> no intersection
             2 -> intersection; any other value shall considered a failure of the computational procedure)
    
    
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
    sc5 = (y_u - t*sin_phi)/cos_phi
    x_sc5 = t*cos_phi - sc5*sin_phi

    # case 6: intersection for y=
    sc6 = (y_l - t*sin_phi)/cos_phi
    x_sc6 = t*cos_phi - sc6*sin_phi

    # initialisation
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
    """
    iPoints : list of intersection points 
    Nx      : number of pixels in x-direction
    Ny      : number of pixels in y-direction
    x_l     : minimum of x coordinate
    x_u     : maximum of x coordinate
    y_l     : minimum of y coordinate
    y_u     : maximum of y coordinate

    debug   : if True -> produce additional output

    returns:

    xIndices : ndarray of x-indices
    yIndices : ndarray of y-indices

    list_nxny : list of discretized start- and end points [n_x1, n_y1, n_x2, n_y2]

    note !!! xIndex[m], yIndices[m] indicate an image point
          (0,0) indicate the lower left corner point of an image
          => in some application transformation to image coordinates may be necessary !!
    """
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
    n_temp = int(round((x1 - x_l)/del_Ax))
    n_x1 = n_temp if n_temp < Nx else Nx-1
    n_temp = int(round((x2 - x_l)/del_Ax))
    n_x2 = n_temp if n_temp < Nx else Nx-1
    # number of pixels in x-direction
    del_nx = n_x2 - n_x1
    if debug:
        print(f"del_nx: {del_nx}")

    # point (x2, y2) -> (dp_x2, dp_y2)
    n_temp = int(round((y1 - y_l)/del_Ay))
    n_y1 = n_temp if n_temp < Ny else Ny-1
    n_temp = int(round((y2 - y_l)/del_Ay))
    n_y2 = n_temp if n_temp < Ny else Ny-1
    # number of pixels in y-direction
    del_ny = n_y2 - n_y1
    if debug:
        print(f"del_ny: {del_ny}")
    
    if del_nx == 0:
        # vertical projection line (Ny points)
        vec_nx = np.repeat(n_x1, Ny)
        vec_ny = np.arange(Ny, dtype=np.int32)
        return vec_nx, vec_ny, [n_x1, 0, n_x1, Ny-1]
    elif del_ny == 0:
        # horizontal projection line (Nx points)
        vec_ny = np.repeat(n_y1, Nx)
        vec_nx = np.arange(Nx, dtype=np.int32)
        return vec_nx, vec_ny, [0, n_y1, Nx-1, n_y1]
    else:
        # projection line with defined slope
        slope_d = (n_y2 - n_y1) / (n_x2 - n_x1)
        if debug:
            print(f"slope_d: {slope_d}")

    # case: -> generate x array first
    if del_nx >= abs(del_ny):
        yStepSize = slope_d
        xIndices = np.arange(n_x1, n_x2+1, dtype=np.int32)
        yIndices = n_y1 + np.int32(np.rint(np.arange(len(xIndices)) * yStepSize)) 
        if debug:
            print(f"yStepSize: {yStepSize}; yStepSize/del_Ay: {yStepSize/del_Ay}")
            
    # -> generate y array first
    else:
        xStepSize = 1.0/abs(slope_d)
        if debug:
            print(f"xStepSize: {xStepSize}; xStepSize/del_Ax: {xStepSize/del_Ax}")
    
        if n_y2 > n_y1:
            yIndices = np.arange(n_y1, n_y2+1, dtype=np.int32)
        else:
            yIndices = np.arange(n_y1, n_y2-1, -1, dtype=np.int32)
            
        xIndices = n_x1 + np.int32(np.rint(np.arange(len(yIndices)) * xStepSize))
    return xIndices, yIndices, [n_x1, n_y1, n_x2, n_y2]

def projectionSingleLine(tVal, phi_deg, imgMat, x_l, x_u, y_l, y_u, placeholder=math.nan, debug=True):
    """
    tVal    : a scalar t-value
    phi_deg : angle of d-Vector with respect to x-axis
    imgMat  : matrix of image; nr of rows: Ny ; nr of columns: Nx; imgMat[0, 0] is image point
              at upper left corner; element of imgMat must be scalar (eg. greyscale image)
    x_l     : minimum of x coordinate
    x_u     : maximum of x coordinate
    y_l     : minimum of y coordinate
    y_u     : maximum of y coordinate
    debug   : if True -> additional output is displayed
    
    returns:

    projectionVal : if the projection line crosses image, the corresponding points on the line are summed, the 
                    result is returned
                    if the projection line is outside the image a placeholder value is returned
    """
    Ny, Nx = imgMat.shape

    # compute intersection
    intersectionCount, intersectionPoints = intersections(tVal, phi_deg, x_l, x_u, y_l, y_u)

    if intersectionCount == 0:
        # projection line outside image
        projectionVal = placeholder    
    elif intersectionCount == 2:
        # if intersection is found -> compute projection line
        xIndices, yIndices, n_points = digitalProjectionLine(intersectionPoints, Nx, Ny, x_l, x_u, y_l, y_u, debug=False)
        projectionVal = np.sum(imgMat[Ny - 1 - yIndices, xIndices])
    else:
        # this should never happen (only 0 or 2 intersections are valid results)
        sys.exit(f"error -> len(intersectionPoints) : {len(intersectionPoints)} ; intersectionPoints: {intersectionPoints}")
    # the projection value
    return projectionVal

def projectionMultiLine(tVec, phi_deg, imgMat, x_l, x_u, y_l, y_u, placeholder=math.nan, debug=True):
    """
    tVec    : t-vector
    phi_deg : angle of d-Vector with respect to x-axis
    imgMat  : matrix of image; nr of rows: Ny ; nr of columns: Nx; imgMat[0, 0] is image point
              at upper left corner; element of imgMat must be scalar (eg. greyscale image)
    x_l     : minimum of x coordinate
    x_u     : maximum of x coordinate
    y_l     : minimum of y coordinate
    y_u     : maximum of y coordinate
    debug   : if True -> additional output is displayed
    
    returns:

    projectionVec : if the projection line crosses image, the corresponding points on the line are summed, the 
                    result is returned
                    if the projection line is outside the image a placeholder value is returned
    """
    projectionVec = []

    # for a fixed angle phi_deg iterate over all values of dVec and compute a projection value for 
    # each element of dVec
    for tVal in tVec:
        projectionVal = projectionSingleLine(tVal, phi_deg, imgMat, x_l, x_u, y_l, y_u, placeholder=placeholder, debug=debug)
        projectionVec.append(projectionVal)
    return np.array(projectionVec)

def sinoGram(tVec, phi_deg_vec, imgMat, x_l, x_u, y_l, y_u, placeholder=math.nan, debug=True):
    """
    tVec        : t-vector
    phi_deg_vec : vector of angles of t-Vector with respect to x-axis
    imgMat  : matrix of image; nr of rows: Ny ; nr of columns: Nx; imgMat[0, 0] is image point
              at upper left corner; element of imgMat must be scalar (eg. greyscale image)
    x_l     : minimum of x coordinate
    x_u     : maximum of x coordinate
    y_l     : minimum of y coordinate
    y_u     : maximum of y coordinate
    debug   : if True -> additional output is displayed
    
    returns:

    sinogram: rows -> angles, columns -> t values
    
    """
    # initialise
    sinogram = np.zeros( (len(phi_deg_vec), len(tVec)) )

    for nrow, phi_deg in enumerate(phi_deg_vec):
        sinogram[nrow, :] = projectionMultiLine(tVec, phi_deg, imgMat, x_l, x_u, y_l, y_u, placeholder=placeholder, debug=True)
    return sinogram


def backprojection(sinogram, x_min, x_max, y_min, y_max, Nx, Ny, t_vec, theta_vec_deg, filter=None):
    """
    sinogram      : matrix (rows: angle, columns: projection values
    x_min         : bounding box of image 
    x_max         : ""
    y_min         : ""
    y_max         : ""
    Nx            : nr of pixels in x-direction
    Ny            : nr of pixels in y-direction
    t_vec         : vector of t-values
    theta_vec_deg : vector of projection angles
    filter        : if defined filter is specified of np.fft.fftfreq(sinogram.shape) values in frequency domain
    """
    rows, cols = sinogram.shape
    if filter is not None:
        if len(filter) != cols:
            sys.exit('length of filter does not match length of projections in sinograms')
        # initialise memory for filtered sinogram
        sinogram_filtered = np.zeros_like(sinogram)
        # for each projection angle, filter in frequency domain and transform back to signal space
        for row in range(rows):
            sinogram_filtered[row, :] = np.real(np.fft.ifft(np.fft.fft(sinogram[row, :]) * filter))

    x_vec = np.linspace(x_min, x_max, Nx)
    y_vec = np.linspace(y_min, y_max, Ny)

    # initialise 
    img_backprojection = np.zeros( (Ny, Nx), dtype=np.float64)
    # change orientation -> needed for image coordinates
    y_vec_r = np.flip(y_vec)

    # iterate over angles
    for theta_index, theta_deg in enumerate(theta_vec_deg):
        theta_rad = math.radians(theta_deg)
        cos_v = math.cos(theta_rad)
        sin_v = math.sin(theta_rad)

        # iterate over image rows
        for nr in range(Ny):
            yval = y_vec_r[nr]
            t_values = x_vec * cos_v + yval * sin_v

            if filter is not None:
                sino_theta = sinogram_filtered[theta_index, :]
            else:
                sino_theta = sinogram[theta_index, :]

            # interpolate
            proj_i = np.interp(t_values, t_vec, sino_theta, left=0, right=0)
            # accumulate interpolated projections
            img_backprojection[nr, :] = img_backprojection[nr, :] + proj_i
        
    if filter is not None:
        return img_backprojection, sinogram_filtered
    else:
        return img_backprojection
    
def backprojection2(sinogram, x_min, x_max, y_min, y_max, Nx, Ny, t_vec, theta_vec_deg, filter=None):
    """
    sinogram      : matrix (rows: angle, columns: projection values
    x_min         : bounding box of image 
    x_max         : ""
    y_min         : ""
    y_max         : ""
    Nx            : nr of pixels in x-direction
    Ny            : nr of pixels in y-direction
    t_vec         : vector of t-values
    theta_vec_deg : vector of projection angles
    filter        : if defined filter is specified of np.fft.fftfreq(sinogram.shape) values in frequency domain
    """
    rows, cols = sinogram.shape
    if filter is not None:
        if len(filter) != cols:
            sys.exit('length of filter does not match length of projections in sinograms')
        # initialise memory for filtered sinogram
        sinogram_filtered = np.zeros_like(sinogram)
        # for each projection angle, filter in frequency domain and transform back to signal space

        for row in range(rows):
            sinogram_row = sinogram[row, :]
            sino_f_tmp = np.real(np.fft.ifft(np.fft.fft(np.fft.fftshift(sinogram_row)) * filter))
            sinogram_filtered[row, :] = np.fft.fftshift(sino_f_tmp)

    x_vec = np.linspace(x_min, x_max, Nx)
    y_vec = np.linspace(y_min, y_max, Ny)

    # initialise 
    img_backprojection = np.zeros( (Ny, Nx), dtype=np.float64)
    # change orientation -> needed for image coordinates
    y_vec_r = np.flip(y_vec)

    # iterate over angles
    for theta_index, theta_deg in enumerate(theta_vec_deg):
        theta_rad = math.radians(theta_deg)
        cos_v = math.cos(theta_rad)
        sin_v = math.sin(theta_rad)

        # iterate over image rows
        for nr in range(Ny):
            yval = y_vec_r[nr]
            t_values = x_vec * cos_v + yval * sin_v

            if filter is not None:
                sino_theta = sinogram_filtered[theta_index, :]
            else:
                sino_theta = sinogram[theta_index, :]

            # interpolate
            proj_i = np.interp(t_values, t_vec, sino_theta, left=0, right=0)
            # accumulate interpolated projections
            img_backprojection[nr, :] = img_backprojection[nr, :] + proj_i
        
    if filter is not None:
        return img_backprojection, sinogram_filtered
    else:
        return img_backprojection