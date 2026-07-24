
# Hex sticker for the 'sdetorus' package: a toroidal diffusion on the torus
# [-pi, pi)^2 -- the stationary density of a Wrapped-Normal Ornstein-Uhlenbeck
# process (jet colormap), its rotational drift field (white arrows), and a few
# diffusion paths spiralling into the mode. Built with the package's own
# dStatWn2D / driftWn2D / euler2D / linesTorus.

library(sdetorus)
library(hexSticker)
library(magick)

# Shared logo standards
font <- "Aller_Rg"
name_size <- 31.2
url_size <- 9.0
url_x <- 1.00
url_y <- 0.08
url_angle <- 30
hex_border <- 1.5
dpi <- 600

# Custom options
border <- "#08519C"
fill <- "#08306B"
ghurl <- "github.com/egarpor/sdetorus"
out_logo <- "logo/logo.png"
out_man <- "man/figures/logo.png"

# Toroidal-diffusion picture (drawn to a square PNG)

# Wrapped-Normal OU with rotational drift A = a*I + w*J (complex eigenvalues ->
# a stable spiral sink); the rotation leaves the stationary density unchanged.
mu <- c(0, 0)
a_att <- 0.7 # attraction toward the mode
w_rot <- 1.3 # rotation strength (swirl)
Arot <- rbind(c(a_att, -w_rot), c(w_rot, a_att))
s_den <- c(1.05, 1.05) # density / drift blob width
s_dif <- c(0.25, 0.25) # trajectory diffusion

# Navy-anchored jet ramp: low density -> navy (blends with the hex fill).
jet_navy <- grDevices::colorRampPalette(
  c("#08306B", "#0B2E8F", "#1E6BFF", "#12C2C2",
    "#5FE84E", "#F5E13C", "#FF7A22", "#D51111"))(128)

# Open a square PNG canvas for the torus picture
torus_png <- tempfile(fileext = ".png")
canvas_px <- 1200
grDevices::png(torus_png, width = canvas_px, height = canvas_px, bg = fill)
op <- graphics::par(mar = rep(0, 4))

# (a) stationary-density heatmap
grid_n <- 340
xth <- seq(-pi, pi, l = grid_n + 1)[-(grid_n + 1)]
grid <- as.matrix(expand.grid(xth, xth))
z <- matrix(dStatWn2D(x = grid, alpha = c(a_att, a_att, 0), mu = mu,
                      sigma = s_den), grid_n, grid_n)
graphics::image(xth, xth, z, col = jet_navy, axes = FALSE, xlab = "", ylab = "",
                xlim = c(-pi, pi), ylim = c(-pi, pi))

# (b) drift vector field: white arrows on a coarse grid (swirl); the drift
# vanishes at the mode, so arrows naturally fade at the bright centre.
gg <- seq(-pi, pi, l = 15)
xy <- as.matrix(expand.grid(gg, gg))
b <- driftWn2D(x = xy, A = Arot, mu = mu, sigma = s_den)
nb <- sqrt(rowSums(b^2))
sc <- 0.30 / max(nb) # longest arrow ~ one grid cell
keep <- nb > 0.2
graphics::arrows(xy[keep, 1], xy[keep, 2],
                 xy[keep, 1] + sc * b[keep, 1], xy[keep, 2] + sc * b[keep, 2],
                 length = 0.05, col = "#FFFFFFBB", lwd = 3)

# (c) a few diffusion sample paths spiralling into the mode
set.seed(5)
th0 <- c(2.3, 3.5, 4.7) # start clear of the name and URL
x0 <- cbind(2.8 * cos(th0), 2.8 * sin(th0))
n_steps <- 250
tr <- euler2D(x0 = x0, A = Arot, mu = mu, sigma = s_dif, N = n_steps,
              delta = 0.01, type = 1)
sub <- seq(1, n_steps + 1, by = 3) # subsample -> smoother line
for (i in seq_len(nrow(x0))) {
  linesTorus(tr[i, 1, sub], tr[i, 2, sub], col = "#08183A", lwd = 9) # casing
  linesTorus(tr[i, 1, sub], tr[i, 2, sub], col = "#FFFFFF", lwd = 4) # path
}

# (d) soft top vignette so the white wordmark reads over the field
yb <- seq(pi * 0.30, pi, l = 80)
for (k in seq_len(length(yb) - 1)) {
  a <- ((yb[k] - pi * 0.30) / (pi - pi * 0.30))^1.6 * 0.85
  graphics::rect(-pi, yb[k], pi, yb[k + 1], border = NA,
                 col = grDevices::rgb(8, 26, 58, alpha = 255 * a,
                                      maxColorValue = 255))
}

# Close the device
graphics::par(op)
grDevices::dev.off()

# Hex sticker + alpha-mask crop to the hexagon

# hexSticker does not clip an oversized subplot, so render (a) the full sticker
# and (b) a hexagon-only silhouette, then use (b) as an alpha mask on (a).
render <- function(subplot, s_width, pkg, url, file) {
  sticker(
    subplot = subplot, s_x = 1, s_y = 1, s_width = s_width, s_height = s_width,
    package = pkg, p_x = 1, p_y = 1.52, p_size = name_size,
    p_color = "#FFFFFF", p_family = font,
    h_fill = fill, h_color = border, h_size = hex_border,
    url = url,
    u_x = url_x, u_y = url_y, u_angle = url_angle, u_size = url_size,
    u_color = "#FFFFFF", u_family = font,
    dpi = dpi, filename = file
  )
  invisible(file)
}

# Render the full-bleed sticker
full_f <- tempfile(fileext = ".png")
render(torus_png, 1.0, "sdetorus", ghurl, full_f)

# Render a hexagon-only silhouette to use as the alpha mask
blank <- image_blank(4, 4, color = "none")
btmp <- tempfile(fileext = ".png")
image_write(blank, btmp, format = "png")
mask_f <- tempfile(fileext = ".png")
render(btmp, 0.1, "", "", mask_f)

# Crop to the hexagon using the silhouette as the alpha channel
sticker_img <- image_composite(image_read(full_f), image_read(mask_f),
                               operator = "CopyOpacity")
image_write(sticker_img, out_logo, format = "png")

# Mirror to man/figures/logo.png
dir.create(dirname(out_man), showWarnings = FALSE, recursive = TRUE)
image_write(sticker_img, out_man, format = "png")
