import os
import sys
import numpy as np

def p2s(arr):
    # transform from physical to spectral
    # in y first
    arr = np.fft.rfft(arr, axis=0, norm="backward")
    # followed by x
    arr = np.fft.fft (arr, axis=1, norm="backward")
    return arr

def s2p(arr):
    # transform from spectral to physical
    # in x first
    arr = np.fft.ifft (arr, axis=1, norm="backward")
    # followed by y
    arr = np.fft.irfft(arr, axis=0, norm="backward")
    return arr

def compute_wavenumber(domain):
    def kernel(nitems, index):
        # e.g. nitems == 8
        #  -> 0 +1 +2 +3 -4 -3 -2 -1
        return index if index < nitems // 2 else index - nitems
    lx = domain["lx"]
    ly = domain["ly"]
    nx = domain["nx"]
    ny = domain["ny"]
    kny = ny // 2 + 1
    # x: normal
    kx = [2. * np.pi / lx * kernel(nx, i) for i in range( nx)]
    # y: Hermite symmetry
    ky = [2. * np.pi / ly * kernel(ny, j) for j in range(kny)]
    kx, ky = np.meshgrid(kx, ky)
    return kx, ky

def gauss(std, x, mean):
    # gaussian distribution
    factor = 1. / std / np.sqrt(2. * np.pi)
    num = (x - mean) ** 2.
    den = 2. * std ** 2.
    return factor * np.exp(-1. * num / den)

def init_grid(domain):
    xs = np.linspace(0., domain["lx"], domain["nx"], endpoint=False)
    ys = np.linspace(0., domain["ly"], domain["ny"], endpoint=False)
    xs, ys = np.meshgrid(xs, ys)
    return xs, ys

def init_scalar(domain, xs, ys):
    # scalar field, gaussian by default
    std = 5.e-1
    arr = 1. \
            * gauss(std, xs, 0.45 * domain["lx"]) \
            * gauss(std, ys, 0.45 * domain["ly"])
    vmin = np.min(arr)
    vmax = np.max(arr)
    arr = (arr - vmin) / (vmax - vmin)
    arr -= np.mean(arr)
    return arr

def sf2vel(domain, sf):
    # convert streamfunction to velocity
    sf = p2s(sf)
    ux = np.zeros((domain["ny"] // 2 + 1, domain["nx"]), dtype=np.complex128)
    uy = np.zeros((domain["ny"] // 2 + 1, domain["nx"]), dtype=np.complex128)
    # ensure divergence free
    kx, ky = compute_wavenumber(domain)
    k2inv = np.power(kx, 2.) + np.power(ky, 2.)
    k2inv[0, 0] = np.inf
    k2inv = 1. / k2inv
    ux = +1j * ky * sf * k2inv
    uy = -1j * kx * sf * k2inv
    ux = s2p(ux)
    uy = s2p(uy)
    return ux, uy

def check_metrics(domain, pux, puy, psc):
    def kernel(title, arr):
        print(f"{title}: (max: {np.max(arr): .1e}), (min: {np.min(arr): .1e}), (sum: {np.sum(arr): .1e})")
    kx, ky = compute_wavenumber(domain)
    div = kx * p2s(pux) + ky * p2s(puy)
    print(f"maximum divergence: {np.max(np.abs(div)): .1e}")
    kernel("ux", pux)
    kernel("uy", puy)
    kernel("sc", psc)

def initialiser0(domain):
    print("coalescing vortices")
    xs, ys = init_grid(domain)
    # initialise stream function
    std = 3.e-1
    sf = \
            + 1. \
            + 3. * gauss(std, xs, 0.4 * domain["lx"]) \
                 * gauss(std, ys, 0.4 * domain["ly"]) \
            + 3. * gauss(std, xs, 0.6 * domain["lx"]) \
                 * gauss(std, ys, 0.6 * domain["ly"])
    ux, uy = sf2vel(domain, sf)
    sc = init_scalar(domain, xs, ys)
    return ux, uy, sc

def initialiser1(domain):
    print("Kelvin-Helmholtz instability")
    xs, ys = init_grid(domain)
    # initialise stream function
    std = 1.e-1
    sf = \
            + 2. * (1. + 0.01 * np.sin(xs)) * gauss(std, ys, 0.25 * domain["ly"]) \
            - 2. * (1. + 0.01 * np.sin(xs)) * gauss(std, ys, 0.75 * domain["ly"])
    ux, uy = sf2vel(domain, sf)
    sc = init_scalar(domain, xs, ys)
    return ux, uy, sc

def initialiser2(domain):
    def compute_vel(kx, ky, x, y):
        def rand(vmin, vmax):
            val = np.random.random_sample()
            return (vmax - vmin) * val + vmin
        psix = kx * x + rand(-np.pi, +np.pi)
        psiy = ky * y + rand(-np.pi, +np.pi)
        norm = np.sqrt(kx * kx + ky * ky)
        ux = + ky / norm * np.cos(psix) * np.sin(psiy)
        uy = - kx / norm * np.sin(psix) * np.cos(psiy)
        return ux, uy
    print("Decaying turbulence")
    xs, ys = init_grid(domain)
    ux = np.zeros(xs.shape, dtype=np.float64)
    uy = np.zeros(ys.shape, dtype=np.float64)
    for _ in range(100):
        kx = 2. * np.pi / domain["lx"] * np.random.randint(4, 12)
        ky = 2. * np.pi / domain["ly"] * np.random.randint(4, 12)
        ux_, uy_ = compute_vel(kx, ky, xs, ys)
        ux += ux_
        uy += uy_
    vmax = max(np.max(ux), np.max(uy))
    vmin = min(np.min(ux), np.min(uy))
    ux /= max(vmax, -vmin)
    uy /= max(vmax, -vmin)
    sc = init_scalar(domain, xs, ys)
    return ux, uy, sc

def initialiser3(domain):
    def compute_vel(kx, ky, x, y):
        psix = kx * x
        psiy = ky * y
        norm = np.sqrt(kx * kx + ky * ky)
        ux = + ky / norm * np.cos(psix) * np.sin(psiy)
        uy = - kx / norm * np.sin(psix) * np.cos(psiy)
        return ux, uy
    print("Taylor-Green vortex")
    xs, ys = init_grid(domain)
    ux = np.zeros(xs.shape, dtype=np.float64)
    uy = np.zeros(ys.shape, dtype=np.float64)
    kx = 2. * np.pi / domain["lx"]
    ky = 2. * np.pi / domain["ly"]
    ux, uy = compute_vel(kx, ky, xs, ys)
    vmax = max(np.max(ux), np.max(uy))
    vmin = min(np.min(ux), np.min(uy))
    ux /= max(vmax, -vmin)
    uy /= max(vmax, -vmin)
    sc = init_scalar(domain, xs, ys)
    return ux, uy, sc

def visualise(pux, puy, psc):
    try:
        from matplotlib import pyplot
    except ModuleNotFoundError:
        print("could not import matplotlib, abort visualiser")
        return
    fig = pyplot.figure()
    ax0 = fig.add_subplot(131)
    ax1 = fig.add_subplot(132)
    ax2 = fig.add_subplot(133)
    ax0.contourf(np.real(pux), extend="both")
    ax1.contourf(np.real(puy), extend="both")
    ax2.contourf(np.real(psc), extend="both")
    ax0.set_title("ux")
    ax1.set_title("uy")
    ax2.set_title("scalar")
    keywords = {
            "xticks": [],
            "yticks": [],
            "aspect": "equal",
    }
    ax0.set(**keywords)
    ax1.set(**keywords)
    ax2.set(**keywords)
    pyplot.show()
    pyplot.close()

def rectify(arr):
    # sometimes array can be in Fortran order,
    #  which is converted to C order
    ny, nx = arr.shape
    vec = np.zeros(nx * ny, dtype=np.complex128)
    for j in range(ny):
        for i in range(nx):
            vec[j * nx + i] = arr[j, i]
    arr = np.reshape(vec, (ny, nx))
    return arr

def main(initialiser):
    domain = {
            "nx": 256,
            "ny": 256,
            "lx": 2. * np.pi,
            "ly": 2. * np.pi,
    }
    # init fields in physical domain
    pux, puy, psc = initialiser(domain)
    # check outcome
    check_metrics(domain, pux, puy, psc)
    visualise(pux, puy, psc)
    # save to files
    # flow fields are in the spectral domain
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
    np.save(f"{root}/glsizes.npy", np.array([domain["nx"], domain["ny"]], dtype=np.uint64))
    np.save(f"{root}/lengths.npy", np.array([domain["lx"], domain["ly"]], dtype=np.float64))
    np.save(f"{root}/ux.npy", rectify(p2s(pux)))
    np.save(f"{root}/uy.npy", rectify(p2s(puy)))
    np.save(f"{root}/sc.npy", rectify(p2s(psc)))

if __name__ == "__main__":
    msg = "give one of [0, 1, 2, 3]"
    argv = sys.argv
    if 2 != len(argv):
        print(msg)
        exit(1)
    try:
        case = int(argv[1])
    except ValueError:
        print(msg)
        exit(1)
    if not case in [0, 1, 2, 3]:
        print(msg)
        exit(1)
    initialisers = (initialiser0, initialiser1, initialiser2, initialiser3)
    main(initialisers[case])

