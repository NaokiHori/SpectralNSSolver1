import os
import sys
import numpy as np

def p2s(arr):
    # transform from physical to spectral
    # in z first
    arr = np.fft.rfft(arr, axis=0, norm="backward")
    # in y first
    arr = np.fft.fft (arr, axis=1, norm="backward")
    # followed by x
    arr = np.fft.fft (arr, axis=2, norm="backward")
    return arr

def s2p(arr):
    # transform from spectral to physical
    # in x first
    arr = np.fft.ifft (arr, axis=2, norm="backward")
    # followed by y
    arr = np.fft.ifft (arr, axis=1, norm="backward")
    # followed by z
    arr = np.fft.irfft(arr, axis=0, norm="backward")
    return arr

def compute_wavenumber(domain):
    def kernel(nitems, index):
        # e.g. nitems == 8
        #  -> 0 +1 +2 +3 -4 -3 -2 -1
        return index if index < nitems // 2 else index - nitems
    lx = domain["lx"]
    ly = domain["ly"]
    lz = domain["lz"]
    nx = domain["nx"]
    ny = domain["ny"]
    nz = domain["nz"]
    knz = nz // 2 + 1
    # x, y: normal
    kx = [2. * np.pi / lx * kernel(nx, i) for i in range( nx)]
    ky = [2. * np.pi / ly * kernel(ny, j) for j in range( ny)]
    # y: Hermite symmetry
    kz = [2. * np.pi / lz * kernel(nz, k) for k in range(knz)]
    kz, ky, kx = np.meshgrid(kz, ky, kx, indexing="ij")
    return kx, ky, kz

def gauss(std, x, mean):
    factor = 1. / std / np.sqrt(2. * np.pi)
    num = (x - mean) ** 2.
    den = 2. * std ** 2.
    return factor * np.exp(-1. * num / den)

def init_grid(domain):
    xs = np.linspace(0., domain["lx"], domain["nx"], endpoint=False)
    ys = np.linspace(0., domain["ly"], domain["ny"], endpoint=False)
    zs = np.linspace(0., domain["lz"], domain["nz"], endpoint=False)
    zs, ys, xs = np.meshgrid(zs, ys, xs, indexing="ij")
    return xs, ys, zs

def init_scalar(domain, xs, ys, zs):
    std = 5.e-1
    arr = 1. \
            * gauss(std, xs, 0.45 * domain["lx"]) \
            * gauss(std, ys, 0.45 * domain["ly"]) \
            * gauss(std, zs, 0.45 * domain["lz"])
    vmin = np.min(arr)
    vmax = np.max(arr)
    arr = (arr - vmin) / (vmax - vmin)
    return arr

def check_metrics(domain, pux, puy, puz, psc):
    kx, ky, kz = compute_wavenumber(domain)
    div = kx * p2s(pux) + ky * p2s(puy) + kz * p2s(puz)
    print(f"maximum divergence: {np.max(np.abs(div)): .1e}")
    print(f"max(ux): {np.max(pux): .1e}, min(ux): {np.min(pux): .1e}")
    print(f"max(uy): {np.max(puy): .1e}, min(uy): {np.min(puy): .1e}")
    print(f"max(uz): {np.max(puz): .1e}, min(uz): {np.min(puz): .1e}")
    print(f"max(sc): {np.max(psc): .1e}, min(sc): {np.min(psc): .1e}")

def initialiser(domain):
    def rand(vmin, vmax):
        val = np.random.random_sample()
        return (vmax - vmin) * val + vmin
    def compute_vel(kx, ky, kz, x, y, z):
        psix = kx * x + rand(-np.pi, +np.pi)
        psiy = ky * y + rand(-np.pi, +np.pi)
        psiz = kz * z + rand(-np.pi, +np.pi)
        t = rand(-np.pi, +np.pi)
        sx = np.sin(t + 0. * 2. * np.pi / 3.)
        sy = np.sin(t + 1. * 2. * np.pi / 3.)
        sz = np.sin(t + 2. * 2. * np.pi / 3.)
        ux = sx / kx * np.sin(psix) * np.cos(psiy) * np.cos(psiz)
        uy = sy / ky * np.cos(psix) * np.sin(psiy) * np.cos(psiz)
        uz = sz / kz * np.cos(psix) * np.cos(psiy) * np.sin(psiz)
        return ux, uy, uz
    print("Decaying turbulence")
    xs, ys, zs = init_grid(domain)
    ux = np.zeros(xs.shape, dtype=np.float64)
    uy = np.zeros(ys.shape, dtype=np.float64)
    uz = np.zeros(zs.shape, dtype=np.float64)
    for _ in range(10):
        kx = 2. * np.pi / domain["lx"] * np.random.randint(1, 4)
        ky = 2. * np.pi / domain["ly"] * np.random.randint(1, 4)
        kz = 2. * np.pi / domain["lz"] * np.random.randint(1, 4)
        ux_, uy_, uz_ = compute_vel(kx, ky, kz, xs, ys, zs)
        ux += ux_
        uy += uy_
        uz += uz_
    vmax = max(np.max(ux), np.max(uy), np.max(uz))
    vmin = min(np.min(ux), np.min(uy), np.min(uz))
    ux /= max(vmax, -vmin)
    uy /= max(vmax, -vmin)
    uz /= max(vmax, -vmin)
    sc = init_scalar(domain, xs, ys, zs)
    return ux, uy, uz, sc

def rectify(arr):
    nz, ny, nx = arr.shape
    vec = np.zeros(nx * ny * nz, dtype=np.complex128)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                vec[k * ny * nx + j * nx + i] = arr[k, j, i]
    return np.reshape(vec, (nz, ny, nx))

def main():
    domain = {
            "nx": 64,
            "ny": 64,
            "nz": 64,
            "lx": 2. * np.pi,
            "ly": 2. * np.pi,
            "lz": 2. * np.pi,
    }
    pux, puy, puz, psc = initialiser(domain)
    check_metrics(domain, pux, puy, puz, psc)
    root = "output"
    try:
        os.mkdir(root)
    except FileExistsError:
        pass
    except:
        msg = "mkdir failed"
        raise RuntimeError(msg)
        exit(1)
    np.save(f"{root}/step.npy", np.array(0 , dtype=np.uint64))
    np.save(f"{root}/time.npy", np.array(0., dtype=np.float64))
    np.save(f"{root}/glsizes.npy", np.array([domain["nx"], domain["ny"], domain["nz"]], dtype=np.uint64))
    np.save(f"{root}/lengths.npy", np.array([domain["lx"], domain["ly"], domain["lz"]], dtype=np.float64))
    np.save(f"{root}/ux.npy", rectify(p2s(pux)))
    np.save(f"{root}/uy.npy", rectify(p2s(puy)))
    np.save(f"{root}/uz.npy", rectify(p2s(puz)))
    np.save(f"{root}/sc.npy", rectify(p2s(psc)))

if __name__ == "__main__":
    main()

