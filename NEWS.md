* `isnonemptypartof(R,S)` yielding `false` for any unsupported argument types
  has been removed as it is better to get an exception in this case.

* Methods, `view`, `getindex` and `setindex!` have been extended so that it
  makes sense to do: `V = view(A,B)`, `C = A[B]` or `A[B] = C` with `A` any
  array and `B` a `CartesianBox`.

* `boundingbox` method can take a predicate function.

* Instances of `CartesianBox{N}` are as fast as `CartesianIndices{N}` or as
  `CartesianRange{CardinalIndex{N}}` (depending on your version of Julia) for
  multi-dimensional loops.  They can thus be used as a *fast* and *portable*
  replacement of these different representations of rectangular
  multi-dimensional regions.

* `CartesianIndices(B)` and, for Julia versions ≤ 0.6, `CartesianIndices(B)` can
  be used to convert `B`, an instance of `CartesianBox{N}` to another
  representation.
