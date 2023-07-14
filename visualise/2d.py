import os
import sys
import numpy as np
import matplotlib
from matplotlib import pyplot

def calc_wavenumbers(nx, ny):
    kny = int(ny // 2) + 1
    kx = np.zeros(nx)
    ky = np.zeros(kny)
    for i in range(nx):
        kx[i] = i if i < nx // 2 else i - nx
    for j in range(kny):
        ky[j] = j if j < ny // 2 else j - ny
    kx, ky = np.meshgrid(kx, ky)
    return kx, ky

def create_mask(kx, ky):
    ny, nx = kx.shape
    mask = np.logical_and(np.abs(kx) < nx // 3, np.abs(ky) < ny // 3)
    return mask

def transform(arr):
    arr = np.fft.ifft(arr, axis=1, norm="backward")
    arr = np.fft.irfft(arr, axis=0, norm="backward")
    return arr

def main(is_2_3, show_scalar):
    root = "output/save"
    dnames = sorted([f"{root}/{dname}" for dname in os.listdir(root) if dname.startswith("step")])
    if show_scalar:
        fig = pyplot.figure(figsize=(8., 6.), facecolor="#000000", edgecolor="#000000")
        axs = fig.add_subplot(121), fig.add_subplot(122)
    else:
        fig = pyplot.figure(figsize=(6., 6.), facecolor="#000000", edgecolor="#000000")
        axs = [fig.add_subplot(111)]
    for cnt, dname in enumerate(dnames):
        # load arrays
        glsizes = np.load(f"{dname}/glsizes.npy")
        nx, ny = glsizes
        ux = np.load(f"{dname}/ux.npy")
        uy = np.load(f"{dname}/uy.npy")
        if show_scalar:
            sc = np.load(f"{dname}/sc.npy")
        # initialisation
        if 0 == cnt:
            kx, ky = calc_wavenumbers(nx, ny)
            if is_2_3:
                mask = create_mask(kx, ky)
                assert mask.shape == ux.shape, f"{mask.shape} != {ux.shape}"
        # vorticity
        vz = 1j * kx * uy - 1j * ky * ux
        # masking
        if is_2_3:
            vz *= mask
            if show_scalar:
                sc *= mask
        # to physical domain (normalised)
        vz = transform(vz)
        if show_scalar:
            sc = transform(sc)
        # visualise
        fig.suptitle(f"{cnt+1} / {len(dnames)}")
        keywords = {
                "xticks": [],
                "yticks": [],
                "aspect": "equal",
        }
        axs[0].clear()
        axs[0].contourf(np.real(vz), cmap="seismic", levels=51)
        axs[0].set(**keywords)
        if show_scalar:
            axs[1].clear()
            axs[1].contourf(np.abs(sc), cmap="rainbow", levels=51)
            axs[1].set(**keywords)
        if len(dnames) - 1 == cnt:
            pyplot.show(block=True)
        else:
            pyplot.show(block=False)
            pyplot.pause(5.e-1)
    pyplot.close()

if __name__ == "__main__":
    matplotlib.rcParams["text.color"] = "#FFFFFF"
    is_2_3 = True
    show_scalar = True
    main(is_2_3, show_scalar)

