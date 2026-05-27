# spatial-small-area-estimation-belgium

A small R sandbox I built to pressure-test the methods I plan to use at Ghent University if my Fulbright proposal is funded. The idea is simple: before I land in Belgium, I want to be sure the spatial-statistics pipeline I described in the proposal actually runs end to end on a Belgian statistical-sector geography, with a synthetic outcome that behaves the way real hypertension data tends to behave.

**Author:** Romario Joseph, MPH (BU SPH, Epidemiology & Biostatistics)
**Stack:** R 4.x, sf, spdep, dplyr, ggplot2
**Data:** Public Statbel statistical-sector shapefile + a synthetic hypertension prevalence variable I simulate inside the script. No real patient data.

## Epidemiological Objective

The clinical question I am rehearsing here is the one I expect to ask in Belgium: does hypertension prevalence cluster geographically at the level of statistical sectors (the finest official small-area unit Statbel publishes), or is it scattered randomly across the country? If it clusters, that is a signal that hypertension burden in Belgium is being shaped by neighborhood-level factors (housing, income, ageing of the local population, food environment) rather than purely individual factors, and that is exactly the equity claim my Fulbright project wants to test on real Belgian Health Interview Survey data.

For this sandbox I am not trying to recover a true effect. I am checking that the pipeline I will hand to a Belgian collaborator on day one is plumbed correctly, that the spatial weights matrix behaves, and that the diagnostics I plan to report are actually estimable on the Statbel geography.

## Methodological Framework

I kept the model list deliberately short so the sandbox is honest about what it can and cannot show.

The headline test is **Global Moran's I**, the standard product-moment statistic for spatial autocorrelation. I build a row-standardized queen-contiguity neighbors list from the sector polygons using `spdep::poly2nb` and `spdep::nb2listw`, and I assess significance with 999 Monte Carlo permutations rather than the asymptotic normal approximation, because the asymptotic null tends to misbehave when the weights matrix has many islands or thin polygons (and the Belgian sector layer has both).

Alongside the global test I run a **Moran scatterplot** so the relationship between each sector's value and its spatial lag is visible, not just summarized in a single statistic. I also fit a baseline **simultaneous autoregressive (SAR) model** with `spdep::spautolm` against a simulated covariate, mostly as a smoke test that the weights object is compatible with the regression machinery I will need later for the real analysis.

I am explicit in the script about what I am not doing yet: I am not running Getis-Ord Gi* local hot-spot detection in this sandbox, I am not doing small-area Bayesian smoothing (BYM2), and I am not adjusting for any covariate that I have not simulated myself. Those are the next steps once I have the real Belgian Health Interview Survey extract in hand at Ghent.

## Data Architecture

The pipeline lives in a single R script (`R/01_belgium_spatial_baseline.R`) so a reviewer can read it top to bottom in one sitting. I structured it as five short stages:

**Stage 1, geometry ingest.** I download the Statbel statistical-sector shapefile (the 2023 release is the one I am using here) directly from the open Statbel data portal, unzip it to a local cache, and read it with `sf::st_read`. I reproject everything to Belgian Lambert 2008 (EPSG:3812) so distances and contiguity are computed in a CRS that is appropriate for Belgium rather than in WGS84.

**Stage 2, geometry cleaning.** Real shapefiles are messy. I run `sf::st_make_valid` on every polygon to repair self-intersections, drop polygons with zero area (digitization artefacts), and check for islands. The number of polygons dropped is logged so the reviewer can see exactly how much was removed.

**Stage 3, synthetic outcome.** I simulate sector-level hypertension prevalence using a spatially-structured logistic mixture. I draw a smooth spatial trend from a Gaussian random field with a moderate range parameter, add a small amount of independent noise, and pass the result through a logit so the simulated proportions live in (0, 1) and have a national mean near 0.28, which is roughly the Belgian adult hypertension prevalence reported in the most recent published HIS round. The point of simulating with spatial structure is to make sure my pipeline can actually detect clustering when clustering is present; if Moran's I came out non-significant on a clearly clustered simulation, that would tell me something is wrong with the weights, not with the data.

**Stage 4, spatial dependence test.** I build the neighbors list, row-standardize it, and run `spdep::moran.test` and `spdep::moran.mc`. I print both the analytic p-value and the permutation p-value so the reader can see they agree. I also save a Moran scatterplot to `outputs/`.

**Stage 5, baseline SAR regression.** I fit a simultaneous autoregressive lag model with `spdep::spautolm` against a simulated covariate (a sector-level "deprivation index" I generate the same way I generated the outcome). This stage exists mainly to confirm the weights object plugs into the regression API cleanly. I am not interpreting the coefficient.

```
spatial-small-area-estimation-belgium/
├── README.md
├── R/
│   └── 01_belgium_spatial_baseline.R   # end-to-end pipeline
├── data/
│   └── README.md                       # Statbel URL, license, SHA of the zip
├── outputs/
│   ├── moran_scatterplot.png
│   └── moran_results.csv
└── LICENSE
```

A few cleaning conventions I want to be explicit about. Sector codes from Statbel are zero-padded character strings; I keep them as character throughout and never coerce to integer, because dropping a leading zero on a join is the kind of silent bug that ruins a small-area analysis. Missing values in the synthetic outcome are not allowed by construction, but I still write a check for `any(is.na(.))` before the Moran step, because in the real analysis I will almost certainly have missing sectors and I want the failure mode to be loud. The neighbors list is saved to disk so subsequent runs do not have to recompute it, but I check that the cached neighbors object matches the current geometry hash before reusing it.

## Reproduce

```r
install.packages(c("sf", "spdep", "dplyr", "ggplot2"))
source("R/01_belgium_spatial_baseline.R")
```

## Why this exists

If a Fulbright reviewer asks whether I can actually run the analysis I am proposing, this repository is the answer. The Belgian geography is real, the spatial-statistics calls are the same ones I will use against the real HIS data, and the script runs in under a minute on a laptop. When I get to Ghent the only thing that changes is the outcome variable.

## Disclaimer

The hypertension prevalence values in this sandbox are simulated. Nothing in this repository should be interpreted as a real estimate of Belgian hypertension burden.

## Contact

Romario Joseph, rjoseph3@bu.edu, [LinkedIn](https://www.linkedin.com/in/romariojosephpublichealth/)
# spatial-small-area-estimation-belgium
Spatial small-area estimation of hypertension prevalence across Belgian statistical sectors: synthetic-data sandbox for the Fulbright Belgium / Ghent University proposal. spdep + Moran's I baseline.
