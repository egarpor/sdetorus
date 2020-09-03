## Test environments

* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (release, devel)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (on R-hub)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit (on R-hub)
* Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit (on R-hub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on R-hub)
* Fedora Linux, R-devel, clang, gfortran (on R-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (on R-hub)
* macOS 10.13.6 High Sierra, R-release, brew (on R-hub)
* macOS 10.13.6 High Sierra, R-release, CRAN's setup (on R-hub)
* Oracle Solaris 10, x86, 32 bit, R-release (on R-hub)
* Oracle Solaris 10, x86, 32 bit, R-release, Oracle Developer Studio 12.6 (on R-hub)

## R CMD check results

0 errors | 0 warnings | 1 note

* "Days since last update: 1" I am resubmitting to fix the additional issues in https://www.stats.ox.ac.uk/pub/bdr/LTO/sdetorus.out

## Resubmission

This is a resubmission. In this version I have:

* Fix the redefinition of arma::mat forwardSweepPeriodicTridiag() as arma::vec.
