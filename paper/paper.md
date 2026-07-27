---
title: 'ioccultcalc: open-source prediction of stellar occultations by asteroids, with local orbit determination'
tags:
  - C++
  - astronomy
  - celestial mechanics
  - asteroids
  - stellar occultations
  - orbit determination
authors:
  - name: Michele Bigi
    orcid: 0009-0006-8787-0499
    affiliation: "1, 2, 3"
affiliations:
  - name: ITALOccult Team, Italy
    index: 1
  - name: Gruppo Astrofili Massesi, Massa, Italy
    index: 2
  - name: Fellow of the Royal Astronomical Society
    index: 3
date: 27 July 2026
bibliography: paper.bib
---

# Summary

When an asteroid passes in front of a star, the shadow it casts sweeps a narrow
path across the Earth. Observers along that path see the star disappear for a
few seconds, and the timings they record measure the asteroid's size and shape
far more precisely than any other ground-based technique. The method depends
entirely on knowing where the shadow will fall, and on knowing how uncertain
that prediction is: a track drawn a hundred kilometres from its true position
sends observers to the wrong place.

`ioccultcalc` predicts these events. It propagates asteroid orbits with a full
force model, searches a Gaia DR3 catalogue along the resulting sky path,
refines each candidate to sub-millisecond closest approach, and projects the
shadow onto the Earth's surface. It also **refines the orbit on its own
astrometric observations**, so that both the trajectory and its uncertainty are
computed locally rather than imported from an external service.

The software is built on `AstDyn`, a C++23 celestial-mechanics library
developed alongside it, and exports predictions in the interchange format used
by Occult4 [@Herald2020] and OccultWatcher [@Pavlov2020], the tools through
which the observing community coordinates its campaigns.

# Statement of need

Occultation prediction has long relied on closed-source software. Occult4 and
OccultWatcher are excellent and freely available, but their internals cannot be
inspected, extended, or cited; predictions are distributed as results rather
than as reproducible computations. Orbital elements and their covariances are
obtained from AstDyS [@Knezevic2003] or JPL Horizons [@Giorgini1996], so the
uncertainty an observer sees is one someone else computed.

For a body whose orbit is well determined this is unproblematic. For faint or
recently discovered objects — precisely those where an occultation would be most
informative — the uncertainty dominates the prediction, and being unable to
examine how it was derived is a real limitation.

`ioccultcalc` addresses this by making the whole chain open and local. An
observer can fit the orbit on the published astrometry, obtain the formal
covariance from that fit, and see the resulting error ellipse — all without
depending on an external service. On asteroid (820987) 2015 BK290 the fit
reaches an RMS of 0.26 arcsec over 90 observations spanning sixteen years, and
the resulting ellipse agrees with the AstDyS one to 9 % in the major axis with
coincident orientation.

The software also predicts occultations by the **satellites of binary
asteroids**, deriving the satellite's position either from an SPK kernel or from
the parameters of its mutual orbit as published in the literature — which is
what exists for most known binaries.

# Implementation

Reference frames and physical units are encoded in the type system. A position
in GCRF and one in ECLIPJ2000 are distinct types and cannot be added; neither
can a distance in astronomical units and one in metres. Frame confusion is a
classic source of silent error in this field, and making it a compile-time
failure removes an entire class of defects.

Six numerical integrators are provided, characterised against each other and
against JPL Horizons. Over a decade of propagation the implicit Radau IIA method
[@Everhart1985] is the most accurate at 2.1e-8 AU; the embedded RKF78 pair is
twenty times faster at 7.9e-7 AU and serves as the default. The five working
methods agree with one another to about 1e-7 AU while differing from Horizons by
about 2e-6 — which locates the residual discrepancy in the force model rather
than in the integration.

Planetary positions are read from JPL SPK kernels through a native, stateless
reader that requires no external toolkit and can be queried concurrently, which
is what makes the parallel corridor scan possible.

# Validation

The library is checked against external references rather than only against
itself, on the principle that internal consistency is not correctness. Every
defect found during development was invisible to the unit tests and emerged only
from comparison with an outside source.

Residuals are compared observation by observation with those AstDyS records in
the same file: on (820987) 2015 BK290 the mean absolute residual is 0.199 arcsec
against 0.169 for AstDyS, with biases of −0.05 and +0.08 arcsec. Positions are
checked against JPL Horizons to better than 1e-5 AU over arcs up to ten years.
Earth rotation is verified from first principles.

Predictions themselves are validated against Occult4 and OccultWatcher. For the
occultation of J175336.90−214720.9 by (1216) Askania on 15 August 2026, the
Besselian elements agree with the reference prediction to the fourth significant
figure in the second-order terms, and the exported file is accepted and
displayed by OccultWatcher.

# Acknowledgements

This work builds on the freely available services of AstDyS at the University of
Pisa and of JPL Horizons, and on the Gaia DR3 catalogue. The interchange format
is documented by David Herald as part of Occult4. The author thanks the
occultation observing community, whose published predictions provided the
external references against which this software was validated.

# References
