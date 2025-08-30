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


## How to use

The tool is a command-line Python script that scans light and flat frame directories, determines optimal calibration centres, and optionally tags FITS headers with the group assignments.

### Installation

Clone the repository and install dependencies:

    git clone https://github.com/adon-phd/rotator.git
    cd rotator
    pip install -r requirements.txt

Dependencies:
- astropy (required) for FITS handling  
- PuLP (optional) for integer linear programming (exact optimisation; otherwise a greedy fallback is used)  
- openpyxl (optional) for detailed XLSX debug reports  

### Basic workflow

1. **Run the optimiser**

    python rotator.py --lights /path/to/LIGHTS --flats /path/to/FLATS --tolerances 2 5 10

   Produces a coverage report showing how many centres are needed at different tolerances.

2. **Choose a tolerance**  
   After the report, you will be prompted to select a tolerance for tagging.  
   Example: entering `2` selects centres that cover all lights within ±2°.

3. **Tag FITS headers (optional)**  

    python rotator.py --lights /path/to/LIGHTS --flats /path/to/FLATS --tolerances 2 --overwrite

   Add `--overwrite` to update FITS headers with the keyword ROTGRP. Without it, existing tags are preserved.

4. **GUI mode (file picker only)**  

    python rotator.py --gui

5. **Dry run**  

    python rotator.py --lights /path/to/LIGHTS --flats /path/to/FLATS --tolerances 2 --dry-run

   Add `--dry-run` to test without modifying FITS files.

### Example output

    --- Tolerance ±2.0° ---
    Chosen centres (optimised): 246.81, 301.08, 310.21, 314.62, 322.39
    Covered unique angles:      25
    Covered light frames:       1290
    UNCOVERED unique angles:    1
    UNCOVERED light frames:     84
    Residual light angles:
      308.20

This shows that with ±2° tolerance, 1374 light frames are reduced to 5 calibration centres, with only a small fraction uncovered.

### Common options

| Option              | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `--lights PATH`     | Root directory of LIGHT frames (required)                                   |
| `--flats PATH`      | Root directory of FLAT frames (optional)                                    |
| `--tolerances N...` | List of tolerances (deg) to test (default: 1 2 5 10 15)                     |
| `--gui`             | Open GUI folder pickers if paths are not provided                           |
| `--overwrite`       | Overwrite existing FITS keyword values when tagging                         |
| `--dry-run`         | Run analysis and report without modifying FITS headers                      |
| `--noninteractive`  | Skip prompt and use `--tolerance` value directly                            |
| `--tolerance N`     | Tolerance (deg) to use with `--noninteractive`                              |
| `--assume-lights`   | Treat files under `--lights` as LIGHT if header type is missing             |
| `--assume-flats`    | Treat files under `--flats` as FLAT if header type is missing               |
| `--existing-flats`  | Supply known flat angles (e.g. `300,315.26` or `300@HA 300@OIII`)           |
| `--debug-nearest`   | Print detailed nearest-neighbour matching and write XLSX debug output       |

### Recommended usage

- **Planning**: Run in dry-run mode to identify the minimal set of angles where flats are required.  
- **Calibration**: Re-run with `--overwrite` to tag light frames, ensuring WBPP or other pipelines group frames correctly.  
- **Large projects**: For datasets spanning many filters and sessions, this tool ensures reproducible centre selection without manual trial-and-error.




---

## License

This project is released under a **non-commercial license**.  
You may use, modify, and share the code for personal or research purposes,  
but **commercial use is prohibited without explicit permission**.


[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License: NC](https://img.shields.io/badge/License-Non--Commercial-red.svg)](./LICENSE)

---

## Cite

```bibtex
@misc{rotator_v1_0_0,
  author       = {P, Adon},
  title        = {Rotator Angle Optimiser for Flats},
  year         = {2025},
  howpublished = {\url{https://github.com/adon-phd/rotator/releases/tag/v1.0.0}},
  note         = {Version 1.0.0, GitHub release}
}
```

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
```
