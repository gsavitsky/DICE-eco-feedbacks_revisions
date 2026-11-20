"""Land/ocean beta (PgC/ppm, C4MIP 1%/yr to 4xCO2, biogeochemically coupled) as a function of Lenton's K_M."""
import os, sys
import pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
KMS = [15, 25, 40, 60, 80, 120, 160, 220, 300, 450]


def f(km):
    import lenton_c4mip_check as C
    dl, do, dppm, _ = C.run(km, False)
    return dict(KM=km, land_beta=dl / dppm, ocean_beta=do / dppm)


if __name__ == "__main__":
    with Pool(8) as p:
        D = pd.DataFrame(p.map(f, KMS))
    D.to_csv(os.path.join(HERE, "km_beta_map.csv"), index=False)
    print(D.round(3).to_string(index=False))
