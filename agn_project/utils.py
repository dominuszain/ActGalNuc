"""Cosmology and helper utilities for the AGN project."""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from .config import H0, OMEGA_M

cosmo = FlatLambdaCDM(H0=H0, Om0=OMEGA_M)
