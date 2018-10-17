# special_orthogonal_group.R

# The special orthogonal group SO(n),
# i.e. the Lie group of rotations in n dimensions.

library(reticulate)
gs <- import_from_path("geomstats", path = ".")

ATOL <- 1e-5

TAYLOR.COEFFS.1.AT.0 <- c(1., 0.,
                        - 1. / 12., 0.,
                        - 1. / 720., 0.,
                        - 1. / 30240., 0.)
TAYLOR.COEFFS.2.AT.0 <- c(1. / 12., 0.,
                        1. / 720., 0.,
                        1. / 30240., 0.,
                        1. / 1209600., 0.)
TAYLOR.COEFFS.1.AT.PI <- c(0., - pi / 4.,
                         - 1. / 4., - pi / 48.,
                         - 1. / 48., - pi / 480.,
                         - 1. / 480.)



GetMaskIFloat <- function(i, n){
  range.n <- gs$backend$arange(n)
  i.float <- gs$backend$cast(array(i), gs$backend$int32)[1]
  mask.i <- gs$backend$equal(range.n, i.float)
  mask.i.float <- gs$backend$cast(mask.i, gs$backend$float32)
  return(mask.i.float)
}

SpecialOrthogonalGroup <- setRefClass("SpecialOrthogonalGroup",
                                      fields = c("n", "dimension"),
                                      methods = list(

                                        initialize = function(n){
                                          stopifnot(n %% 1 == 0)
                                          stopifnot(n > 0)
                                          .self$n <- n
                                          .self$dimension <- (n * (n - 1)) / 2
                                        },

                                        GetIdentity = function(){
                                          "Get the identity of the group,
                                          as a vector if point_type == 'vector',
                                          as a matrix if point_type == 'matrix'."

                                        },

                                        Belongs = function(point){
                                          "Evaluate if a point belongs to SO(n)."
                                          # for point type vector

                                          point <- ToNdarray(point, to.ndim = 2)
                                          vec.dim <- c(dim(point)[1], dim(point)[2])

                                          return(vec.dim == .self$dimension)
                                        },

                                        Regularize = function(point){
                                          "In 3D, regularize the norm of the rotation vector,
                                          to be between 0 and pi, following the axis-angle
                                          representation's convention.
                                          If the angle angle is between pi and 2pi,
                                          the function computes its complementary in 2pi and
                                          inverts the direction of the rotation axis."

                                          point <- ToNdarray(point, to.ndim = 2)
                                          stopifnot(.self$Belongs(point = point))
                                          n.points <- dim(point)[1]

                                          regularized.point <- point
                                          if (.self$n == 3){
                                            angle <- gs$backend$linalg$norm(regularized.point, axis = 1)
                                            mask.0 <- gs$backend$isclose(angle, pi)
                                            mask.not.0 <- (-mask.0 - 1)
                                            mask.pi <- gs$backend$isclose(angle, pi)

                                            mask.0.float <- gs$backend$cast(mask.0, gs$backend$float32)
                                            mask.not.0.float <- gs$backend$cast(mask.not.0, gs$backend$float32)
                                            mask.pi.float <- gs$backend$cast(mask.pi, gs$backend$float32)

                                            k <- floor(angle / (2 * pi) + .5)

                                            norms.ratio <- gs$backend$zeros_like(angle)

                                            # This avoids division by 0.
                                            angle <- angle + mask.0.float * 1.

                                            norms.ratio <- norms.ratio + mask.not.0.float * (1 - 2 * pi * k / angle)
                                            norms.ratio <- norms.ratio + mask.0.float * 1.
                                            norms.ratio <- norms.ratio + mask.pi.float * (pi / angle - (1 - 2 * pi * k / angle))

                                            regularized.point <- gs$backend$einsum(
                                              'n,ni->ni', norms.ratio, regularized.point)

                                          }
                                          stopifnot(length(dim(regularized.point)) == 2)
                                          return(regularized.point)
                                        },

                                        RegularizeTangentVecAtIdentity = function(tangent.vec, metric=None){
                                          "In 3D, regularize a tangent_vector by getting its norm at the identity,
                                          determined by the metric, to be less than pi."

                                          tangent.vec <- ToNdarray(tangent.vec, to.ndim = 2)

                                          if (.self$n == 3){
                                            tangent.vec.metric.norm <- MetricNorm(tangent.vec) # where is this func defined?
                                            tangent.vec.canonical.norm <- gs$backend$linalg$norm(tangent.vec, axis = 1)
                                            if (length(dim(tangent.vec.canonical.norm)) == 1){
                                              tangent.vec.canonical.norm <- gs$backend$expand_dims(tangent.vec.canonical.norm, axis = 1)
                                            }
                                            mask.norm.0 <- gs$backend$isclose(tangent.vec.metric.norm, 0)
                                            mask.canonical.norm.0 <- gs$backend$isclose(tangent.vec.canonical.norm, 0)

                                            mask.0 <- mask.norm.0 | mask.canonical.norm.0 # Don't know if r understands |
                                            mask.else <- (-mask.0 - 1)

                                            mask.0.float <- gs$backend$cast(mask.0, gs$backend$float32)
                                            mask.else.float <- gs$backend$cast(mask.else, gs$backend$float32)

                                            coef <- gs$backend$zeros_like(tangent.vec.metric.norm)

                                            regularized.vec <- gs$backend$zeros_like(tangent.vec)

                                            regularized.vec <- regularized.vec + mask.0.float * tangent.vec

                                            # This avoids dividing by 0.
                                            tangent.vec.canonical.norm <- tangent.vec.canonical.norm + mask.0.float * 1

                                            coef <- coef + mask.else.float * (tangent.vec.metric.norm / tangent.vec.canonical.norm)
                                            regularized.vec <- regularized.vec + mask.else.float * .self$regularize(coef * tangent.vec)
                                            # This avois dividing by 0.
                                            coef <- coef + mask.0.float * 1
                                            regularized.vec <- mask.else.float * (regularized.vec / coef)
                                          }
                                          return(regularized.vec)
                                        },


                                          RegularizeTangentVec = function(tangent.vec, base.point, metric=None){
                                            "In 3D, regularize a tangent_vector by getting the norm of its parallel
                                            transport to the identity, determined by the metric,
                                            to be less than pi."
                                            tangent.vec <- ToNdarray(tangent.vec, to.ndim = 2)
                                            if (.self$n == 3){
                                              if (metric == None){
                                                metric <- .self$left.canonical.metric # Check this
                                              }
                                              base.point <- .self$Regularize(base.point)

                                              jacobian <- .self$jacobian_translation(  # Where is this func defined?
                                                point = base.point,
                                                left.or.right = metric$left.or.right)
                                              inv.jacobian <- gs$backend$linalg.inv(jacobian)
                                              tangent.vec.at.id <- gs$backend$einsum(
                                                'ni,nij->nj',
                                                tangent.vec,
                                                gs$backend$transpose(inv.jacobian, axes = c(0, 2, 1)))

                                              tangent.vec.at.id <- .self$RegularizeTangentVecAtIdentity(
                                                tangent.vec.at.id,
                                                metric)

                                              regularized.tangent.vec <- gs$backend$einsum(
                                                'ni,nij->nj',
                                                tangent.vec.at.id,
                                                gs$backend$transpose(jacobian,
                                                             axes = c(0, 2, 1)))
                                            }
                                            return(regularized.tangent.vec)

                                          },

                                        SkewMatrixFromVector = function(vec){
                                          "In 3D, compute the skew-symmetric matrix,
                                          known as the cross-product of a vector,
                                          associated to the vector vec."
                                          vec <- ToNdarray(vec, to.ndim = 2)
                                          n.vecs <- dim(vec)[1]
                                          vec.dim <- dim(vec)[2]

                                          # TODO(nina): Change gs.cast function for elementary types
                                          vec.dim <- gs$backend$cast(array(dim(vecv)), gs$backend$float32)[1]
                                          mat.dim <- gs$backend$cast(((1 + sqrt(1 + 8 * vec.dim)) / 2), gs$backend$int32)

                                          skew.mat <- gs$backend$zeros((n.vecs) + (.self$n) * 2)
                                          if (.self$n == 3){
                                            levi.civita.symbol <- array(data = c(
                                               c(0, 0, 0),
                                               c(0, 0, 1),
                                               c(0, -1, 0),
                                               c(0, 0, -1),
                                               c(0, 0, 0),
                                               c(1, 0, 0),
                                               c(0, 1, 0),
                                               c(-1, 0, 0),
                                               c(0, 0, 0)
                                               ), dim = c(3, 3, 3)
                                              )

                                            basis.vec.1 <- c(1, 0, 0)
                                            basis.vec.2 <- c(0, 1, 0)
                                            basis.vec.3 <- c(0, 0, 1)
                                            cross.prod.1 <- gs$backend$einsum(
                                              'ijk,ni,nj->nk',
                                              levi.civita.symbol,
                                              basis.vec.1,
                                              vec)
                                            cross.prod.2 <- gs$backend$einsum(
                                              'ijk,ni,nj->nk',
                                              levi.civita.symbol,
                                              basis.vec.2,
                                              vec)
                                           cross.prod.3 <- gs$backend$einsum(
                                              'ijk,ni,nj->nk',
                                             levi.civita.symbol,
                                              basis.vec.3,
                                              vec)

                                            cross.prod.1 <- ToNdarray(cross.prod.1, to.ndim = 3, axis = 1)
                                            cross.prod.2 <- ToNdarray(cross.prod.2, to.ndim = 3, axis = 1)
                                            cross.prod.3 <- ToNdarray(cross.prod.3, to.ndim = 3, axis = 1)
                                            skew.mat <- gs$backend$concatenate(
                                              c(cross.prod.1, cross.prod.2, cross.prod.3), axis = 1)
                                          }
                                          stopifnot(length(dim(skew.mat)) == 3)
                                          return(skew.mat)
                                        },
                                        VectorFromSkewMatrix = function(skew.mat){
                                          "In 3D, compute the vector defining the cross product
                                          associated to the skew-symmetric matrix skew mat.

                                          In nD, fill a vector by reading the values
                                          of the upper triangle of skew_mat."

                                          skew.mat <- ToNdarray(skew.mat, to.ndim = 3)
                                          n.skew.mats <- dim(skew.mat)[1]
                                          mat.dim.1 <- dim(skew.mat)[2]
                                          mat.dim.2 <- dim(skew.mat)[3]

                                          stopifnot(mat.dim.1 == mat.dim.2)
                                          stopifnot(mat.dim.2 == .self$n)

                                          vec.dim <- .self$dimension
                                          vec <- gs$backend$zeros(c(n.skew.mats, vec.dim))

                                          if (.self$n == 3){
                                            vec.1 <- ToNdarray(skew.mat[ , 2, 1], to.ndim = 2, axis = 1)
                                            vec.2 <- ToNdarray(skew.mat[ , 0, 2], to.ndim = 2, axis = 1)
                                            vec.3 <- ToNdarray(skew.mat[ , 1, 0], to.ndim = 2, axis = 1)
                                            vec <- gs$backend$concatenate(c(vec.1, vec.2, vec.3), axis = 1)
                                          }
                                         stopifnot(length(dim(vec)) == 2)
                                         return(vec)
                                        },

                                        RotationVectorFromMatrix = function(rot.mat){
                                          "In 3D, convert rotation matrix to rotation vector
                                          (axis-angle representation).

                                          Get the angle through the trace of the rotation matrix:
                                          The eigenvalues are:
                                          1, cos(angle) + i sin(angle), cos(angle) - i sin(angle)
                                          so that: trace = 1 + 2 cos(angle), -1 <= trace <= 3

                                          Get the rotation vector through the formula:
                                          S_r = angle / ( 2 * sin(angle) ) (R - R^T)

                                          For the edge case where the angle is close to pi,
                                          the formulation is derived by going from rotation matrix to unit
                                          quaternion to axis-angle:
                                          r = angle * v / |v|, where (w, v) is a unit quaternion.

                                          In nD, the rotation vector stores the n(n-1)/2 values of the
                                          skew-symmetric matrix representing the rotation."

                                          rot.mat <- ToNdarray(rot.mat, to.ndim = 3)
                                          n.rot.mats <- dim(rot.mat)[1]
                                          mat.dim.1 <- dim(rot.mat)[2]
                                          mat.dim.2 <- dim(rot.mat)[3]
                                          stopifnot(mat.dim.1 == mat.dim.2)
                                          stopifnot(mat.dim.1 == .self$n)

                                          if (.self$n == 3){
                                            trace <- gs$backend$trace(rot.mat, axis1 = 1, axis2 = 2)
                                            trace <- ToNdarray(trace, to.ndim = 2, axis = 1)
                                            stopifnot(dim(shape) == c(n.rot.mats, 1, dim(trace)))

                                            cos.angle <- .5 * (trace - 1)
                                            cos.angle <- gs$backend$clip(cos.angle, -1, 1)
                                            angle <- gs$backend$arccos(cos.angle)

                                            rot.mat.transpose <- gs$backend$transpose(rot.mat, axes = c(0, 2, 1))
                                            rot.vec <- .self$VectorFromSkewMatrix(rot.mat - rot.mat.transpose)

                                            mask.0 <- gs$backend$isclose(angle, 0)
                                            mask.0.float <- gs$backend$cast(mask.0, gs$backend$float32)

                                            rot.vec <- rot.vec * (1 + mask.0.float * (.5 - (trace - 3) / 12 - 1))

                                            mask.pi <- gs$backend$isclose(angle, pi)
                                            mask.pi.float <- gs$backend$cast(mask.pi, gs$backend$float32)

                                            mark.else <- ~mask_0 & ~mask_pi # What is this python notation?
                                            mark.else.float <- gs$backend$cast(mark.else, gs$backend$float32)

                                            mask.pi <- gs$backend$squeeze(mask.pi, axis = 1)

                                            # choose the largest diagonal element
                                            # to avoid a square root of a negative number
                                            rot.mat.pi <- gs$backend$einsum(
                                              'ni,njk->njk',
                                              mask.pi.float,
                                              rot.mat)

                                            a <- array(0)
                                            rot.mat.pi.00 <- ToNdarray(
                                              rot.mat.pi[ , 0, 0], to.ndim = 2, axis = 1)
                                            rot.mat.pi.11 <- ToNdarray(
                                              rot.mat.pi[ , 1, 1], to.ndim = 2, axis = 1)
                                            rot.mat.pi.22 <- ToNdarray(
                                              rot.mat.pi[ , 2, 2], to.ndim = 2, axis = 1)
                                            rot.mat.pi.diagonal <- gs$backend$hstack(
                                              c(rot.mat.pi.00,
                                                rot.mat.pi.11,
                                                rot.mat.pi.22))
                                            a <- gs$backend$argmax(rot.mat.pi.diagonal, axis = 1)[1]
                                            b <- (a + 1) %% 3
                                            c <- (a + 2) %% 3

                                            # compute the axis vector
                                            sq.root <- gs$backend$zeros(c(n.rot.mats, 1))

                                            aux <- sqrt(
                                              mask.pi.float * (
                                                rot.mat[ , a, a]
                                                - rot.mat[ , b, b]
                                                - rot.mat[ , c, c]) + 1)
                                            sq.root.pi <- gs$backend$einsum(
                                              'ni,nk->ni',
                                              mask.pi.float,
                                              aux)

                                            sq.root <- sq.root + sq.root.pi

                                            rot.vec.pi <- gs$backend$zeros(c(n.rot.mats, .self$dimension))
                                            mask.a.float <- GetMaskIFloat(a, 3)
                                            mask.b.float <- GetMaskIFloat(b, 3)
                                            mask.c.float <- GetMaskIFloat(c, 3)

                                            mask.a.float <- ToNdarray(mask.a.float, to.ndim = 2, axis = 1)
                                            mask.b.float <- ToNdarray(mask.b.float, to.ndim = 2, axis = 1)
                                            mask.c.float <- ToNdarray(mask.c.float, to.ndim = 2, axis = 1)


                                            mask.a.float <- gs$backend$transpose(mask.a.float)
                                            mask.b.float <- gs$backend$transpose(mask.b.float)
                                            mask.c.float <- gs$backend$transpose(mask.c.float)

                                            mask.a.float <- gs$backend$tile(mask.a.float, c(n.rot.mats, 1))
                                            mask.b.float <- gs$backend$tile(mask.b.float, c(n.rot.mats, 1))
                                            mask.c.float <- gs$backend$tile(mask.c.float, c(n.rot.mats, 1))

                                            rot.vec.pi <- rot.vec.pi + mask.pi.float * mask.a.float * sq.root / 2

                                            # This avoids division by 0.
                                            sq.root <- sq.root + mask.0.float * 1
                                            sq.root <- sq.root + mask.else.float * 1

                                            rot.vec.pi.b <- gs$backend$zeros_like(rot.vec.pi)
                                            rot.vec.pi.c <- gs$backend$zeros_like(rot.vec.pi)

                                            rot.vec.pi.b <- rot.vec.pi.b + gs$backend$einsum(
                                              'nk,ni->nk',
                                              mask.b.float,
                                              ((rot.mat[ , b, a]
                                                + rot.mat[ , a, b])
                                               / (2 * sq.root)))

                                            rot.vec.pi <- rot.vec.pi + mask.pi.float * gs$backend$einsum(
                                              'ni,nk->nk',
                                              mask.pi.float,
                                              rot.vec.pi.b)

                                            rot.vec.pi.c <- rot.vec.pi.c + gs$backend$einsum(
                                              'nk,ni->nk',
                                              mask.c.float,
                                              ((rot.mat[ , c, a]
                                                + rot.mat[ , a, c])
                                               / (2 * sq.root)))

                                            rot.vec.pi <- rot.vec.pi + mask.pi.float * gs$backend$einsum(
                                              'ni,nk->nk',
                                              mask.pi.float,
                                              rot.vec.pi.c)

                                           norm.rot.vec.pi <- norm.rot.vec.pi + gs$backend$linalg$norm(rot.vec.pi, axis = 1)

                                           # This avoids division by 0.
                                           norm.rot.vec.pi <- norm.rot.vec.pi + gs$backend$squeeze(mask.0.float, axis = 1)
                                           norm.rot.vec.pi <- norm.rot.vec.pi + gs$backend$squeeze(mask.else.float, axis = 1)

                                           rot.vec <- rot.vec + mask.pi.float * (
                                             gs$backend$einsum(
                                                 'nk,n->nk',
                                                 angle * rot.vec.pi,
                                                 1 / norm.rot.vec.pi)
                                           )

                                           # This avoid dividing by zero
                                           angle <- angle + mask.0.float

                                           angle <- ToNdarray(angle, to.ndim = 2, axis = 1)
                                           fact <- gs$backend$einsum(
                                             'ni,ni->ni',
                                             mask.else.float,
                                             (angle / (2 * sin(angle)) - 1))

                                           rot.vec <- rot.vec * (1 + fact)

                                          }



                                          return(.self$Regularize(rot.vec))

                                        },

                                        MatrixFromRotationVector = function(rot.vec){
                                          "Convert rotation vector to rotation matrix."

                                          stopifnot(.self$Belongs(rot.vec))
                                          rot.vec <- .self$Regularize(rot.vec)
                                          n.rot.vecs <- dim(rot.vec)[1]

                                          if (.self$n == 3){
                                           angle <- gs$backend$linalg$norm(rot.vec, axis = 1)
                                           angle <- ToNdarray(angle, to.ndim = 2, axis = 1)

                                           skew.rot.vec <- .self$SkewMatrixFromVector(rot.vec)

                                           coef.1 <- gs$backend$zeros_like(angle)
                                           coef.2 <- gs$backend$zeros_like(angle)

                                           mask.0 <- gs$backend$isclose(angle, 0)
                                           mask.0.float <- gs$backend$cast(mask.0, gs$backend$float32)

                                           coef.1 <- coef.1 + mask.0.float * (1 - (angle ** 2) / 6)
                                           coef.2 <- coef.2 + mask.0.float * (1 / 2 - angle ** 2)

                                           mask.else <- -mask.0 - 1
                                           mask.else.float <- gs$backend$cast(mask.else, gs$backend$float32)

                                           # This avoids division by 0.
                                           angle <- angle + mask.0.float * 1

                                           coef.1 <- coef.1 + mask.else.float * (sin(angle) / angle)
                                           coef.2 <- mask.else.float * (
                                             (1 - cos(angle)) / (angle ** 2))

                                           term.1 <- gs$backend$zeros((n.rot.vecs) + (.self$n) * 2)
                                           term.2 <- gs$backend$zeros_like(term.1)

                                           coef.1 <- gs$backend$squeeze(coef.1, axis = 1)
                                           coef.2 <- gs$backend$squeeze(coef.2, axis = 1)
                                           term.1 <- (gs$backend$eye(.self$dimension)
                                                     + gs$backend$einsum('n,njk->njk', coef.1, skew.rot.vec))

                                           squared.skew.rot.vec <- gs$backend$einsum(
                                             'nij,njk->nik', skew.rot.vec, skew.rot.vec)

                                           term.2 <- gs$backend$einsum('n,njk->njk', coef.2, squared.skew.rot.vec)

                                           rot.mat <- term.1 + term.2

                                          }
                                          return(rot.mat)
                                        },

                                        Compose = function(point.1, point.2){
                                          "Compose two elements of SO(n)."

                                          point.1 <- .self$MatrixFromRotationVector(point.1)
                                          point.2 <- .self$MatrixFromRotationVector(point.2)
                                          point.prod <- gs$backend$einsum('ijk,ikl->ijl', point.1, point.2)
                                          point.prod <- .self$RotationVectorFromMatrix(point.prod)
                                          point.prod <- .self$Regularize(point.prod)
                                          return(point.prod)

                                        },

                                        Inverse = function(point){
                                          "Compute the group inverse in SO(n)."
                                          if (.self$n == 3){
                                            inv.point <- -.self$Regularize(point)
                                            return(inv.point)
                                          }
                                          inv.point <- gs$backend$linalg$inv(point)
                                          inv.point <- .self$RotationVectorFromMatrix(inv.point)
                                          return(inv.point)

                                        },

                                        JacobianTranslation = function(point, left.or.right='left'){
                                          "Compute the jacobian matrix of the differential
                                          of the left/right translations from the identity to point in SO(n)."
                                          stopifnot(left.or.right == ("left" || "right"))

                                          if (.self$n == 3){
                                            point <- .self$Regularize(point)
                                            n.points <- dim(point)[1]
                                            angle <- gs$backend$linalg$norm(point, axis = 1)
                                            angle <- gs$backend$expand_dims(angle, axis = 1)

                                            coef.1 <- gs$backend$zeros(c(n.points, 1))
                                            coef.2 <- gs$backend$zeros(c(n.points, 1))

                                            mask.0 <- gs$backend$isclose(angle, 0)
                                            mask.0.float <- gs$backend$cast(mask.0, gs$backend$float32)

                                            coef.1 <- coef.1 + mask.0.float * (
                                              TAYLOR.COEFFS.1.AT.0[1]
                                              + TAYLOR.COEFFS.1.AT.0[3] * angle ** 2
                                              + TAYLOR.COEFFS.1.AT.0[5] * angle ** 4
                                              + TAYLOR.COEFFS.1.AT.0[7] * angle ** 6)


                                            coef.2 <- coef.2 + mask.0.float * (
                                              TAYLOR.COEFFS.2.AT.0[1]
                                              + TAYLOR.COEFFS.2.AT.0[3] * angle ** 2
                                              + TAYLOR.COEFFS.2.AT.0[5] * angle ** 4
                                              + TAYLOR.COEFFS.2.AT.0[7] * angle ** 6)

                                            mask.pi <- gs$backend$isclose(angle, pi)
                                            mask.pi.float <- gs$backend$cast(mask.pi, gs$backend$float32)

                                            delta.angle <- angle - pi

                                            coef.1 <- coef.1 + mask.pi.float * (
                                              TAYLOR.COEFFS.1.AT.PI[2]
                                              + TAYLOR.COEFFS.1.AT.PI[3] * angle ** 2
                                              + TAYLOR.COEFFS.1.AT.PI[4] * angle ** 3
                                              + TAYLOR.COEFFS.1.AT.PI[5] * angle ** 4
                                              + TAYLOR.COEFFS.1.AT.PI[6] * angle ** 5
                                              + TAYLOR.COEFFS.1.AT.PI[7] * angle ** 6)

                                            angle <- angle + mask.0.float
                                            coef.2 <- coef.2 + mask.pi.float * (
                                              (1 - coef.1) / angle ** 2)

                                            mask.else <- ~mask_0 & ~mask_pi # I don't know what
                                            mark.else.float <- gs$backend$cast(mark.else, gs$backend$float32)

                                            # This avoids division by 0.
                                            angle <- angle + mask.pi.float
                                            coef.1 <- coef.1 + mask.else.float * (
                                              (angle / 2) / tan(angle / 2))
                                            coef.2 <- coef.2 + mark.else.float * (
                                              (1 - coef.1) / angle ** 2)
                                            jacobian <- gs$backend$zeros(c(n.points, .self$dimension, .self$dimension))
                                            n.points.tensor <- array(n.points)

                                            for(i in range(n.points)[1]:range(n.points)[2]){
                                              mask.i.float <- GetMaskIFloat(i, n.points.tensor)
                                            }

                                            sign <- -1
                                            if (left.or.right == 'left'){
                                              sign <- 1
                                            }

                                            jacobian.i <- (
                                              coef.1[i] * gs$backend$eye(.self$dimension)
                                              + coef.2[i] * gs$backend$outer(point[i], point[i])
                                              + sign * .self$SkewMatrixFromVector(point[i]) / 2)
                                            jacobian.i <- gs$backend$squeeze(jacobian.i, axis = 0)

                                            jacobian <- jacobian + gs$backend$einsum(
                                              'n,ij->nij',
                                              mask.i.float,
                                              jacobian.i)
                                          }
                                            return(jacobian)
                                          },
                                          RandomUniform = function(n.samples = 1){
                                            "Sample in SO(n) with the uniform distribution."
                                            random.point <- runif(n.samples, .self$dimension) * 2 - 1
                                            random.point <- .self$Regularize(random.point)
                                            return(random.point)

                                          },
                                          GroupExpFromIdentity = function(tangent.vec){
                                            "Compute the group exponential of the tangent vector at the identity."
                                            point <- ToNdarray(tangent.vec, to.ndim = 2)
                                            return(point)
                                          },

                                          GroupLogFromIdentity = function(point){
                                            "Compute the group logarithm of the point at the identity."
                                            tangent.vec <- .self$Regularize(point)
                                            return(tangent.vec)
                                          }
                                      )
  )
