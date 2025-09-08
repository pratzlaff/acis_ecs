import argparse
import astropy.io.fits
import itertools
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy.polynomial import polynomial

# from https://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent

def polyfit2d_old(x, y, z, orders):
    degs = [ order+1 for order in orders ]
    ncols = math.prod(degs)
    G = np.zeros((x.size, ncols))
    ij = itertools.product(*(range(deg) for deg in degs))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z, rcond=-1)
    m = np.reshape(m, degs)
    return m

def polyval2d_old(x, y, m):
    ij = itertools.product(range(m.shape[0]), range(m.shape[1]))
    z = np.zeros_like(x)
    for (i,j) in ij:
        z += m[i,j] * x**i * y**j
    return z

def polyfit2d(x, y, f, deg):
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = polynomial.polyvander2d(x, y, deg)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f, rcond=-1)[0]
    return c.reshape(deg+1)

def read_img(img):
    hdulist = astropy.io.fits.open(img)
    data = hdulist[0].data
    hdulist.close()
    return data

def surface_fit(args):

    # shape = (16, 96)
    flux = np.ma.masked_invalid(read_img(args.img))
    xvals = np.arange(flux.shape[1])
    yvals = np.arange(flux.shape[0])

    # S1, S3 ignored (backside illuminated)
    # for i in 1, 3:
    #     xvals[16*i:16*i+16] = np.nan

    # remove S1,3
    mask = np.logical_and(
        np.logical_or(xvals<16*1,xvals>16*1+15),
        np.logical_or(xvals<16*3,xvals>16*3+15),
    )
    xvals = xvals[mask]
    flux = flux[:,mask]

    yvals = yvals*64. + 32.5
    xvals = xvals*64. + 32.5

    # insert 18 pixels between each S chip
    xvals += (xvals/1024).astype(int) *18
    print(xvals)

    # tile/repeat coordinates, flatten everything
    xvals = np.tile(xvals, flux.shape[0])
    yvals = yvals.repeat(flux.shape[1])
    flux = flux.flatten()

    nx, ny = 100, 100
    xx, yy = np.meshgrid(np.linspace(xvals.min(), xvals.max(), nx), 
                         np.linspace(yvals.min(), yvals.max(), ny))

    if False:
        m = polyfit2d_old(xvals, yvals, flux, [2, 2])
        zz = polyval2d_old(xx, yy, m)
    else:
        m = polyfit2d(xvals, yvals, flux, [2,2])
        zz = polynomial.polyval2d(xx, yy, m)
    print(m)

    # observed / modeled
    ratio = flux / polynomial.polyval2d(xvals, yvals, m)

    fig = plt.figure(figsize=plt.figaspect(0.25))
    ax = fig.add_subplot(1, 2, 1, projection='3d')

    ax.scatter(xvals, yvals, flux)
    surf = ax.plot_surface(xx, yy, zz, cmap=cm.viridis)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('flux')

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.scatter(xvals, yvals, ratio)
    ax.plot_surface(np.array([[0,6300], [0,6300]]), np.array([[0,0],[1024,1024]]), np.array([[1,1],[1,1]]))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('observed / modeled')
    plt.show()

    # Plot
    plt.scatter(xvals, yvals, c=flux)
    #plt.imshow(zz, origin='lower')#extent=(xvals.min(), xvals.max(), yvals.max(), yvals.min()), origin='lower')
    plt.imshow(zz, extent=(xvals.min(), xvals.max(), yvals.max(), yvals.min()))#, origin='lower')
    plt.show()

    plt.scatter(xvals, yvals, c=ratio)
    plt.show()

def main():
    parser = argparse.ArgumentParser(
        description='Perform surface fit of ECS flux.'
    )
    parser.add_argument('img', help='FITS image containing flux measurements for each 64x64 pixel tile.')
    args = parser.parse_args()
    surface_fit(args)

if __name__ == '__main__':
    main()
