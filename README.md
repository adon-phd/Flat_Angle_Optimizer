# Rotator Angle Optimiser for Flats

Modern astrophotography projects can generate thousands of light frames taken at slightly different rotator angles. Without careful planning, ensuring that every light has a matching flat becomes tedious and error-prone, often requiring redundant calibration work. This tool addresses that gap by applying optimisation methods to automatically determine the minimal set of calibration centres needed for complete coverage. The result is a reproducible and efficient workflow that reduces calibration overhead while preserving scientific integrity.

## Why this tool is useful

Calibration of astronomical imaging data requires flat frames that closely match the rotator angle of each light frame in order to correct detector and optical artefacts. Modern imaging campaigns often produce hundreds or thousands of light frames distributed across many slightly different rotator angles. Acquiring flats at every exact angle is infeasible, so a tolerance (for example ±2°) is typically adopted. The central problem is to determine a reduced set of “centre” angles such that every light frame is covered by at least one flat within the specified tolerance. At present this task is usually addressed informally, through manual grouping or ad hoc binning, which is inefficient and prone to error in large, multi-filter datasets.

## Technical challenges

- **Clustering of discrete values:** Light frames often accumulate at many distinct angles separated by small increments. Selecting representative centres requires collapsing these clusters while retaining complete coverage.  
- **Overlapping tolerance windows:** Each centre represents an angular interval, so coverage analysis must account for overlaps and avoid gaps.  
- **Efficiency trade-offs:** Too many centres increase calibration effort, while too few risk leaving frames without appropriate flats.  
- **Consistent integration:** Once centres are selected, metadata must be updated consistently so that downstream processing recognises the grouping.

## Approach

We cast the selection of calibration centres as a **coverage optimisation problem** closely related to the classical set cover problem [Karp, 1972]. The objective is to minimise the number of centres while ensuring that every light frame is assigned to at least one centre within a specified tolerance. This task is NP-hard in general, motivating the use of heuristic strategies [Vazirani, 2001]. Our implementation applies fuzzy merging to collapse near-identical clusters, snaps candidate centres to existing flats, and allows frame counts to weight optimisation decisions. The resulting set of centres provides efficient coverage with minimal redundancy. FITS headers are then tagged with group labels, ensuring consistent metadata for subsequent calibration.  

This work adapts concepts from combinatorial optimisation [Korte & Vygen, 2018] and clustering [Jain, 2010] into a domain-specific workflow for astronomical calibration. By formalising flat planning as a reproducible optimisation task, we provide a practical method that reduces calibration overhead while ensuring robust coverage across large datasets.

---

## References

```bibtex
@incollection{karp1972reducibility,
  title={Reducibility among combinatorial problems},
  author={Karp, Richard M.},
  booktitle={Complexity of Computer Computations},
  pages={85--103},
  year={1972},
  publisher={Springer}
}

@book{vazirani2001approximation,
  title={Approximation Algorithms},
  author={Vazirani, Vijay V.},
  year={2001},
  publisher={Springer}
}

@book{korte2018combinatorial,
  title={Combinatorial Optimization: Theory and Algorithms},
  author={Korte, Bernhard and Vygen, Jens},
  edition={6},
  year={2018},
  publisher={Springer}
}

@article{jain2010data,
  title={Data clustering: 50 years beyond K-means},
  author={Jain, Anil K.},
  journal={Pattern Recognition Letters},
  volume={31},
  number={8},
  pages={651--666},
  year={2010},
  publisher={Elsevier}
}
