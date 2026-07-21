# Hexagonal sticker for the 'sdetorus' package.
#
# The badge shows what the package computes: a *toroidal diffusion*. On the flat
# torus [-pi, pi)^2 we draw the stationary density of a Wrapped-Normal
# Ornstein-Uhlenbeck process (jet colormap, the sdetorus identity), the drift
# vector field (white arrows, a central swirl), and a few diffusion sample paths
# spiralling toward the mode. All of it is produced with the package's own
# functions (dStatWn2D / driftWn2D / euler2D / linesTorus / matlab.like.colorRamps).
#
# Run from the package root:
#
#   Rscript logo/logo.R
#
# Output: logo/logo.png (master) and man/figures/logo.png (shipped mirror).
# Requires: sdetorus, hexSticker, magick.

library(sdetorus)
library(hexSticker)
library(magick)

## ---- Shared logo standard (identical across egarpor packages) -------------
# Aller_Rg is bundled with (and auto-registered by) hexSticker.
FONT    <- "Aller_Rg"   # typeface for the package name and the GitHub URL
P_SIZE  <- 31.2         # package-name size (shared across packages)
U_SIZE  <- 9.0          # GitHub URL size (large enough to read)
U_X     <- 1.00         # GitHub URL position: along the lower-right hex edge
U_Y     <- 0.08
U_ANGLE <- 30
H_SIZE  <- 1.5          # hexagon border thickness
DPI     <- 600

border  <- "#08519C"    # deep jet-blue hexagon border
fill    <- "#08306B"    # dark navy hexagon fill
ghurl   <- "github.com/egarpor/sdetorus"
out_logo <- "logo/logo.png"
out_man  <- "man/figures/logo.png"

## ---- 1. The toroidal-diffusion picture (drawn to a square PNG) --------------

# Wrapped-Normal OU with a rotational drift matrix A = a*I + w*J (complex
# eigenvalues -> a stable spiral sink). The rotation does not change the
# stationary law, so the density is still the centred blob dStatWn2D(alpha).
mu    <- c(0, 0)
a_att <- 0.7                                  # attraction toward the mode
w_rot <- 1.3                                  # rotation (swirl); w/a large -> visible spiral
Arot  <- rbind(c(a_att, -w_rot), c(w_rot, a_att))   # eig a +/- w i: stable spiral sink
s_den <- c(1.05, 1.05)                        # blob width (density + drift shape)
s_dif <- c(0.25, 0.25)                        # trajectory diffusion (smoother paths)

# Navy-anchored jet ramp: low density -> navy (blends with the hex fill).
jet_navy <- grDevices::colorRampPalette(
  c("#08306B", "#0B2E8F", "#1E6BFF", "#12C2C2",
    "#5FE84E", "#F5E13C", "#FF7A22", "#D51111"))(128)

torus_png <- tempfile(fileext = ".png")
W <- 1200
grDevices::png(torus_png, width = W, height = W, bg = fill)
op <- graphics::par(mar = rep(0, 4))

# (a) stationary-density heatmap
M    <- 340
xth  <- seq(-pi, pi, l = M + 1)[-(M + 1)]
grid <- as.matrix(expand.grid(xth, xth))
z    <- matrix(dStatWn2D(x = grid, alpha = c(a_att, a_att, 0), mu = mu, sigma = s_den),
               M, M)
graphics::image(xth, xth, z, col = jet_navy, axes = FALSE, xlab = "", ylab = "",
                xlim = c(-pi, pi), ylim = c(-pi, pi))

# (b) drift vector field: white arrows on a coarse grid (swirl); the drift
# vanishes at the mode, so arrows naturally fade at the bright centre.
gg <- seq(-pi, pi, l = 15)
xy <- as.matrix(expand.grid(gg, gg))
b  <- driftWn2D(x = xy, A = Arot, mu = mu, sigma = s_den)
nb <- sqrt(rowSums(b^2))
sc <- 0.30 / max(nb)                 # longest arrow ~ one grid cell
keep <- nb > 0.2
graphics::arrows(xy[keep, 1], xy[keep, 2],
                 xy[keep, 1] + sc * b[keep, 1], xy[keep, 2] + sc * b[keep, 2],
                 length = 0.05, col = "#FFFFFFBB", lwd = 3)

# (c) a few diffusion sample paths spiralling into the mode
set.seed(5)
th0 <- c(2.3, 3.5, 4.7)              # start on the left/lower-left, clear of name & URL
x0  <- cbind(2.8 * cos(th0), 2.8 * sin(th0))
NN  <- 250
tr  <- euler2D(x0 = x0, A = Arot, mu = mu, sigma = s_dif, N = NN, delta = 0.01,
               type = 1)
sub <- seq(1, NN + 1, by = 3)        # subsample -> smoother line
for (i in seq_len(nrow(x0))) {
  linesTorus(tr[i, 1, sub], tr[i, 2, sub], col = "#08183A", lwd = 9)   # dark casing
  linesTorus(tr[i, 1, sub], tr[i, 2, sub], col = "#FFFFFF", lwd = 4)   # white path
}

# (d) soft top vignette so the white wordmark reads over the field
yb <- seq(pi * 0.30, pi, l = 80)
for (k in seq_len(length(yb) - 1)) {
  a <- ((yb[k] - pi * 0.30) / (pi - pi * 0.30))^1.6 * 0.85
  graphics::rect(-pi, yb[k], pi, yb[k + 1],
                 col = grDevices::rgb(8, 26, 58, alpha = 255 * a, maxColorValue = 255),
                 border = NA)
}

graphics::par(op)
grDevices::dev.off()

## ---- 2. Hex sticker + alpha-mask crop to the hexagon -----------------------

# hexSticker does not clip an oversized subplot, so we render (a) the full sticker
# and (b) a hexagon-only silhouette, then use (b) as an alpha mask to crop (a).
render <- function(subplot, s_width, pkg, url, file) {
  sticker(subplot = subplot, s_x = 1, s_y = 1,
          s_width = s_width, s_height = s_width,
          package = pkg, p_x = 1, p_y = 1.52, p_size = P_SIZE,
          p_color = "#FFFFFF", p_family = FONT,
          h_fill = fill, h_color = border, h_size = H_SIZE,
          url = url, u_x = U_X, u_y = U_Y, u_angle = U_ANGLE,
          u_size = U_SIZE, u_color = "#FFFFFF", u_family = FONT,
          filename = file, dpi = DPI)
  invisible(file)
}

full_f <- tempfile(fileext = ".png")
render(torus_png, 1.0, "sdetorus", ghurl, full_f)

blank  <- image_blank(4, 4, color = "none")
btmp   <- tempfile(fileext = ".png")
image_write(blank, btmp, format = "png")
mask_f <- tempfile(fileext = ".png")
render(btmp, 0.1, "", "", mask_f)               # pure hex silhouette

sticker_img <- image_composite(image_read(full_f), image_read(mask_f),
                               operator = "CopyOpacity")
image_write(sticker_img, out_logo, format = "png")

## ---- 3. Mirror to the conventional man/figures/logo.png --------------------

dir.create(dirname(out_man), showWarnings = FALSE, recursive = TRUE)
image_write(sticker_img, out_man, format = "png")

message("Wrote ", out_logo, " and ", out_man,
        " (", image_info(sticker_img)$width, "x",
        image_info(sticker_img)$height, ")")
