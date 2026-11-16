# Default parameters for droplet position generation.
# Every value here can be overridden from the command line (see `python main.py -h`).

# boundary conditions
#   "triperiodic" : periodic in x, y, z
#   "channel"     : periodic in x, y; walls at z = 0 and z = lz
bc = "channel"

# domain and mesh (no hidden coupling: d = ppd*dx with dx = lx/nx)
lx = 16.
ly = 16.
lz = 10.
nx = 256
ny = 256
nz = 160
ppd = 32          # grid points per droplet diameter

# volume fraction (number of drops = floor(alpha*V/Vdrop))
alpha = 0.1

# method for drop position generation: "monte-carlo" or "boxed"
method = "boxed"
epsf = 0.5        # boxed only: jitter amplitude as a fraction of available slack, in [0, 1]

# minimum clearances, all surface-to-surface and in units of droplet radius R
# droplets never cross or touch any face; gap_* adds extra clearance per axis
gap_x = 0.12      # 0 would allow exact tangency with the face
gap_y = 0.12
gap_z = 0.12      # minimum droplet-surface to z-face (wall) gap / R
drop_gap = 0.4    # minimum droplet-droplet surface gap / R (0.4 -> centre spacing 1.2*d)

# monte-carlo controls
max_attempts = 100000   # attempts per droplet before declaring a restart
n_restarts = 20         # full restarts before giving up

# reproducibility (None -> random seed)
seed = None

# write/plot
outfile = "spheres.in"   # one line per sphere: x y z r
do_write = True
do_plot = False
