# Barry
Making splines in barycentric coordinates.

If one could write splines in barycentric coordinates, one could make splines
in higher dimensions that rely on simplices rather than quad patches (e.g.
NURBS). This would make pretty much everything about the underlying topology
better.

This requires some extra mathematical machinery. Its mostly just linear
algebra:

    - The polynomial coefficients have extra constraints so that polynomial in 
      lambda are also polynomial in x.
    - Must be able to convert between lambda sets (coordinates based on
      different simplices).
    - Link simplices across boundaries.  In 1D these are vertices, in 2D edges,
      in 3D faces, etc.  Boundaries always separate exactly two simplices and
      have one of the barycentric coordinates set to zero.
    - Express BCs across higher dimensional things than a point.
    - Define "basis splines" as a set of constraints rather than an
      interpolation scheme.


"Sounds good doesn't it?  Yes it does, other Barry."
    - Barry Dylan
