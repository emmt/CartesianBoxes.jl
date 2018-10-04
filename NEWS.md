* `boundingbox` method can take a predicate function.

* Instances of `CartesianBox{N}` are as fast as `CartesianIndices{N}` or as
  `CartesianRange{CardinalIndex{N}}` (depending on your version of Julia) for
  multi-dimensional loops.  They can thus be used as a *fast* and *portable*
  replacement of these different representations of rectangular
  multi-dimensional regions.

* `CartesianIndices(B)` and, for Juli versions ≤ 0.6, `CartesianIndices(B)` can
  be used to convert `B`, an instance of `CartesianBox{N}` to another
  representation.
