# Changelog

## 1.1.0 - 2026-07-23

- Preserve multi-chromosome cross imports by carrying the active project stem and
  linkage-map filename between `Rmap` and `Rcross`.
- Warn before `Rcross` uses its legacy one-chromosome fallback when a selected
  input dataset has no available linkage map.
- Synchronize the GUI filename stem and Rmap/Rcross file defaults from
  `qtlcart.rc`, with compatibility fallback to the latest QTL Cartographer log.
- Add a five-chromosome regression workflow using the original Nuzhdin example.
- Add explicit Windows application and package version metadata.

## 1.0.0 - 2026-07-22

- Initial Windows GUI and x64 port of the complete QTL Cartographer 1.17j suite.
