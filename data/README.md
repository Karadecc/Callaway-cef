# Data

Digitised from Q. J. Li et al., Phys. Rev. B 87, 214408 (2013), arXiv:1305.6853,
by colour segmentation of screenshots (scripts in `tools/`). Axes were calibrated on
the detected tick marks (log axes: least-squares over all major ticks). Curves are
the column-wise median of each colour, cleaned against a running median and averaged
in log-spaced bins. Overlays of the extraction on the original figures were checked visually; they are not
included because they reproduce the published figures.

| file | source | content |
|---|---|---|
| `y2ti2o7_kappa_H0.csv` | Fig. 1(b) | Y2Ti2O7, H = 0, 0.3-87 K |
| `tbyti2o7_kappa_H0.csv` | Fig. 1(b) | TbYTi2O7, H = 0 |
| `tb2ti2o7_sampleB_kappa_H0.csv` | Fig. 1(b) | Tb2Ti2O7 sample B, H = 0 |
| `tb2ti2o7_kappa.csv` | Fig. 2(a) | Tb2Ti2O7 sample A, H // [111] = 0, 3, 8.5, 11.5, 14 T (heat current // [111]) |
| `tb2ti2o7_ratio_vs_H.csv` | Fig. 2(c) | kappa(H)/kappa(0), H // [111], at T = 0.36, 0.5, 0.7, 0.97, 1.4, 1.95 and 5 K (legend of Fig. 2(c)); columns T, H, ratio. Points with ratio ~ 1 at low field overlap and are partly missing. |

Consistency check: the ratio kappa(H)/kappa(0) at 0.36 K computed from the Fig. 2(a)
curves agrees with the independently digitised Fig. 2(c) within 2-6 % (3, 8.5 and
11.5 T). Samples A and B differ by ~15 % at 0.36 K, consistent with the authors'
statement that kappa does not depend on sample size.

Typical digitisation uncertainty: a few % in kappa, ~1 % in T. The field data in
Fig. 2(a) and the zero-field data of the two samples should not be mixed in one fit
without allowing for that sample-to-sample difference.
