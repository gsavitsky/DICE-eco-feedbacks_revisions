"""C4MIP beta/gamma (1%/yr to 4xCO2, BGC-COU) for Lenton variants (iii) K_M=120 and (iv) K_M=145 with ocean warming."""
import os, sys
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def f(a):
    import lenton_c4mip_check as C
    km, ow, coupled = a
    return a, C.run(km, coupled, ocean_warming=ow)


if __name__ == "__main__":
    jobs = [(km, ow, c) for km, ow in ((120.0, False), (145.0, True)) for c in (False, True)]
    with Pool(4) as p:
        r = dict(p.map(f, jobs))
    for km, ow in ((120.0, False), (145.0, True)):
        lb, ob, dppm, dT = r[(km, ow, False)]; lc, oc, _, _ = r[(km, ow, True)]
        print(f"{'iv' if ow else 'iii'} K_M={km:g}: land beta {lb/dppm:.2f} ocean beta {ob/dppm:.2f} | land gamma {(lc-lb)/dT:.1f} ocean gamma {(oc-ob)/dT:.1f}")
