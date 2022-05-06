# User visible changes in `CartesianBoxes`

## Branch 0.3

* Method `intersection(A,B)` has been deprecated.  Call
  `intersect(CartesianBox,A,B)` instead.

* Method `ranges` has been deprecated.  Call `CartesianBoxes.indices` instead.


## Branch 0.2

* Type `CartesianBox` is now a simple decoration wrapped around an instance of
  `CartesianIndices`.  This simplifies a lot the code and lessen the risk of
  incompatibilities with future Julia versions. However, compatibility with
  Julia < 0.7 has to be abandoned.

* Cartesian boxes are now abstract arrays.  Type `CartesianBox{N}` is a subtype
  of `AbstractArray{CartesianIndex{N},N}`.  Method `axes(B)` which used to
  yield the list of index ranges in a Cartesian box `B` now yield the list of
  index ranges to index `B`.  Call `ranges(B)` to get the list of index ranges
  in the Cartesian box `B`.

* Expressions `A ∩ B` and `A ⊆ B` (or equivalently `intersect(A,B)` and
  `issubset(A,B)`) yield a Cartesian box when at least one of `A` or `B` is a
  Cartesian box while the other is a Cartesian box, a list of integer valued
  ranges, a list of dimensions, an instance of `CartesianIndices`, or an
  (abstract) array.


## Branch 0.1

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

* Constructors of `CartesianBox{N}` no longer accept a mixture of dimension
  lengths and index intervals.  This suppress the ambiguity that integer `n` is
  interpreted as a dimension length when all other arguments are integers,
  while it is interpreted as the range `n:n` when any other argument is an
  index interval.
