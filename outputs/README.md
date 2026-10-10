# Runtime outputs

This directory contains generated search, Packmol, LAMMPS, and post-processing
artifacts. New runs are stored below `runs/<run_id>/` and the latest run is
announced by `latest_run.json`.

The files are intentionally ignored by Git because trajectories, molecular
coordinates, logs, and intermediate data can be large and are tied to a
particular machine, executable, force-field configuration, and run.

Use `data/README.md` for the per-evaluation feedback directory. The source
repository tracks the input generators, analyzers, tests, configuration, and
documentation; generated numerical results should be archived separately with
their run metadata when they are ready for publication.
