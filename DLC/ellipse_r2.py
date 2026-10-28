import numpy as np

def get_residuals_ellipse_fit(ellipse_params,data_pts):
    r2=0;
    npts = data_pts.shape[0]
    #rotation of the ellipse
    ellipse_theta = ellipse_params[4]
    #major axis
    a=np.max(ellipse_params[2:4])
    #minor axis
    b=np.min(ellipse_params[2:4])
    #eccentricity of the ellipse
    e = np.sqrt(1-(b ** 2)/(a ** 2))
    d=0
    for i in range(npts):
        # shift data points so that the center of the ellipse is in the origin
        xp = data_pts[i, 0] - ellipse_params[0]
        yp = data_pts[i, 1] - ellipse_params[1]
        #transform data points to polar coordinates
        rp = np.sqrt(xp ** 2 + yp ** 2)
        thetap = np.arctan2(yp, xp)-ellipse_theta
        #radius of the ellipse at the same angle where the point is
        r_ellipse = b / np.sqrt(1 - (e*np.cos(thetap))**2)
        #distance is the difference in polar radius of the point and the ellipse radius
        d = d + (rp - r_ellipse)**2
    return d

def get_best_fit_ellipse(ellipse_params,data_pts):
    nellipse = ellipse_params.shape[0]
    #init R2 to the area of the ellipse
    r2min=np.median(ellipse_params[:,2])*np.median(ellipse_params[:,3])*np.pi
    minindex = -1
    for i in range(nellipse):
        r2i=get_residuals_ellipse_fit(ellipse_params[i],data_pts)
        if r2i<r2min:
            r2min=r2i
            minindex = i
    return (r2min, minindex)
