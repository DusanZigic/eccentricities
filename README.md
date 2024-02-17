collection of scripts that calculates reaction plane angles and eccentricities
and averages rotated energy densities and temperatures

run.py
    python script that extracts temperature evolutions, renames files, compiles
    c++ code and runs executables;
    use: 'python3 run.py --help' to see available parameters

cpsin
    directory that contains c++ source code that calculates reaction plane angles,
    Psi_n based on equation 3 from arxiv:1212.1008;

cespn
    directory that contains c++ source code that calculates anisotropies, eps_n
    based on equation 2 from arxiv:1212.1008;

cavgevols
    directory that contains c++ source code that rotates temperature and energy
    density evolutions for Psin_n and averages them;