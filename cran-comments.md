## Test environments

* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel, release)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (on R-hub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on R-hub)
* Fedora Linux, R-devel, clang, gfortran (on R-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (on R-hub)

NOTEs double-checked and correctly spelled:

- "Possibly mis-spelled words in DESCRIPTION"
- "Found the following (possibly) invalid URLs"

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## CRAN changes

* Exported and documented the dBvM function to avoid :::.
* Replaced \dontrun{} with \donttest{} as much as possible in examples.
* Removed instances of "<<-"".
* Drop the plotVmdSurface3D() function.

