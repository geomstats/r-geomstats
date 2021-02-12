kATOL <- 1e-5

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

#' The special orthogonal group SO(n),
#' i.e. the Lie group of rotations in n dimensions.
SpecialOrthogonalGroup <- setRefClass("SpecialOrthogonalGroup",
  fields = c("n", "dimension"),
  methods = list(
    initialize = function(n){
      stopifnot(n %% 1 == 0)
      stopifnot(n > 0)
      .self$n <- n
      .self$dimension <- (n * (n - 1)) / 2
    },

    Belongs = function(point){
      "Evaluate if a point belongs to SO(n)."
      # for point type vector
      if (length(dim(point)) == 1) {
        point <- ToNdarray(point, to.ndim = 2)
      }
      vec.dim <-  dim(point)[length(dim(point))]
      return(vec.dim == .self$dimension)
    },

    Regularize = function(point){
      "In 3D, regularize the norm of the rotation vector,
      to be between 0 and pi, following the axis-angle
      representation's convention.
      If the angle angle is between pi and 2pi,
      the function computes its complementary in 2pi and
      inverts the direction of the rotation axis."
      if (length(dim(point)) == 1){
      point <- ToNdarray(point, to.ndim = 2)
      }
      n.points <- dim(point)[1]
      regularized.point <- point

      if (.self$n == 3) {
        angle <- apply(regularized.point, 1, function(x){sqrt(sum(x ^ 2))})
        mask.0 <- (abs(angle) < .Machine$double.eps ^ 0.5)
        mask.not.0 <- !mask.0
        mask.pi <- (abs(angle - pi) < .Machine$double.eps ^ 0.5)

        k <- floor(angle / (2 * pi) + .5)

        norms.ratio <- array(0, c(n.points, 1))

        # This avoids division by 0.
        angle <- angle + mask.0 * 1.

        norms.ratio[mask.not.0] <- 1 - 2 * pi * k / angle
        norms.ratio[mask.0] <- 1
        norms.ratio[mask.pi] <-  (pi / angle[mask.pi])

        regularized.point <- regularized.point * c(norms.ratio)
        if (is.null(dim(regularized.point))) {
          regularized.point <- ToNdarray(array(regularized.point), to.ndim = 2)
        }
      }
      regularized.point <- array(regularized.point, dim = c(n.points, 3))
      stopifnot(length(dim(regularized.point)) == 2)
      return(regularized.point)
    },

    SkewMatrixFromVector = function(vec){
      "In 3D, compute the skew-symmetric matrix,
      known as the cross-product of a vector,
      associated to the vector vec."
      if (length(dim(vec)) == 1) {
        vec <- ToNdarray(vec, to.ndim = 2)
      }
      n.vecs <- dim(vec)[1]
      vec.dim <- dim(vec)[2]
      if (.self$n == 3) {
        levi.civita.symbol <- array(data = c(
          c(0, 0, 0),
          c(0, 0, -1),
          c(0, 1, 0),
          c(0, 0, 1),
          c(0, 0, 0),
          c(-1, 0, 0),
          c(0, -1, 0),
          c(1, 0, 0),
          c(0, 0, 0)
        ), dim = c(3, 3, 3)
        )

        ApplyLeviCivitaSymbol <- function(vec) {
          array(cbind(
          vec %*% levi.civita.symbol[1, , ],
          vec %*% levi.civita.symbol[2, , ],
          vec %*% levi.civita.symbol[3, , ]),
          dim=c(3,3))
        }

        cross.prod.1 <- vec %*% ApplyLeviCivitaSymbol(array(c(1,0,0)))
        cross.prod.2 <- vec %*% ApplyLeviCivitaSymbol(array(c(0,1,0)))
        cross.prod.3 <- vec %*% ApplyLeviCivitaSymbol(array(c(0,0,1)))

        skew.mat <- -array(c(cross.prod.1,
                            cross.prod.2,
                            cross.prod.3), dim = c(n.vecs, .self$n, .self$n))
      }
      stopifnot(length(dim(skew.mat)) == 3)
      return(skew.mat)
    },

    JacobianTranslation = function(point, left.or.right="left"){
      "Compute the jacobian matrix of the differential
      of the left/right translations from the identity to point in SO(n)."

      stopifnot(left.or.right == "left" || "right")

      if (.self$n == 3) {
        point <- .self$Regularize(point)
        n.points <- dim(point)[1]
        angle <- apply(point, 1, function(x){sqrt(sum(x ^ 2))})
        angle <- ToNdarray((array(angle)), to.ndim = 2, axis = 1)

        coef.1 <- array(0, dim = c(n.points, 1))
        coef.2 <- array(0, dim = c(n.points, 1))

        mask.0 <- (angle < .Machine$double.eps ^ 0.5)

        coef.1[mask.0] <- (
          TAYLOR.COEFFS.1.AT.0[1]
          + TAYLOR.COEFFS.1.AT.0[3] * angle[mask.0] ^ 2
          + TAYLOR.COEFFS.1.AT.0[5] * angle[mask.0] ^ 4
          + TAYLOR.COEFFS.1.AT.0[7] * angle[mask.0] ^ 6)

        coef.2[mask.0] <- (
          TAYLOR.COEFFS.2.AT.0[1]
          + TAYLOR.COEFFS.2.AT.0[3] * angle ^ 2
          + TAYLOR.COEFFS.2.AT.0[5] * angle ^ 4
          + TAYLOR.COEFFS.2.AT.0[7] * angle ^ 6)

        mask.pi <- (abs(angle - pi) < .Machine$double.eps ^ 0.5)

        delta.angle <- angle - pi

        coef.1[mask.pi] <- (
          TAYLOR.COEFFS.1.AT.PI[2] * delta.angle
          + TAYLOR.COEFFS.1.AT.PI[3] * delta.angle ^ 2
          + TAYLOR.COEFFS.1.AT.PI[4] * delta.angle ^ 3
          + TAYLOR.COEFFS.1.AT.PI[5] * delta.angle ^ 4
          + TAYLOR.COEFFS.1.AT.PI[6] * delta.angle ^ 5
          + TAYLOR.COEFFS.1.AT.PI[7] * delta.angle ^ 6)
        coef.2[mask.pi] <- (1 - coef.1[mask.pi]) / angle[mask.pi] ^ 2

        mask.else <- !mask.0 & !mask.pi
        angle <- angle + mask.pi
        coef.1[mask.else] <- ((angle[mask.else] / 2)
                              / tan(angle[mask.else] / 2))
        coef.2[mask.else] <- ((1 - coef.1[mask.else])
                              / angle[mask.else] ^ 2)
        jacobian <- array(0, dim = c(n.points, .self$dimension, .self$dimension))

        for (i in 1:n.points) {
          sign <- -1
          if (left.or.right == "left") {
            sign <- 1
          }

          jacobian[i, , ] <- (
            coef.1[i] * diag(1, .self$dimension, .self$dimension)
            + coef.2[i] * outer(point[i, ], point[i, ])
            + sign * (.self$SkewMatrixFromVector(ToNdarray(array(point[1, ]), to.ndim = 2)) / 2)[1,,])
        }
      }
      stopifnot(length(dim(jacobian)) == 3)
      return(jacobian)
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
      vec <- array(0, dim = c(n.skew.mats, vec.dim))

      if (.self$n == 3) {
        vec.1 <- ToNdarray(array(skew.mat[ , 3, 2]), to.ndim = 2, axis = 1)
        vec.2 <- ToNdarray(array(skew.mat[ , 1, 3]), to.ndim = 2, axis = 1)
        vec.3 <- ToNdarray(array(skew.mat[ , 2, 1]), to.ndim = 2, axis = 1)
        vec <- cbind(c(vec.1, vec.2, vec.3))
      }
      vec <- array(vec, dim = c(n.skew.mats, 3))
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
      if (length(dim(rot.mat)) == 2) {
        rot.mat <- ToNdarray(rot.mat, to.ndim = 3)
      }
      n.rot.mats <- dim(rot.mat)[1]
      mat.dim.1 <- dim(rot.mat)[2]
      mat.dim.2 <- dim(rot.mat)[3]
      stopifnot(mat.dim.1 == mat.dim.2)
      stopifnot(mat.dim.1 == .self$n)

      if (.self$n == 3) {
        trace <- colSums(apply(rot.mat, 1, diag))

        stopifnot(dim(trace) == c(n.rot.mats, 1))

        clip <- function(x, x.min, x.max){
          x[x < x.min] <- x.min
          x[x > x.max] <- x.max
          return(x)
        }

        cos.angle <- .5 * (trace - 1)
        cos.angle <- clip(cos.angle, -1, 1)
        angle <- acos(cos.angle)

        rot.mat.transpose <- aperm(rot.mat, c(1, 3, 2))

        rot.vec <- .self$VectorFromSkewMatrix(rot.mat - rot.mat.transpose)

        mask.0 <- (abs(angle) < .Machine$double.eps ^ 0.5)
        rot.vec <- rot.vec * (1 + mask.0 * (.5 - (trace - 3) / 12 - 1))
        mask.pi <- (abs(angle - pi) < .Machine$double.eps ^ 0.5)

        mask.else <- !mask.0 & !mask.pi
        # choose the largest diagonal element
        # to avoid a square root of a negative number

        rot.mat.pi <- apply(rot.mat, c(2, 3), function(x){x * mask.pi})
        rot.mat.pi <- ToNdarray(rot.mat.pi, to.ndim = 3)

        a <- array(0)
        rot.mat.pi.00 <- ToNdarray(
          array(rot.mat.pi[ , 1, 1]), to.ndim = 2, axis = 1)
        rot.mat.pi.11 <- ToNdarray(
          array(rot.mat.pi[ , 2, 2]), to.ndim = 2, axis = 1)
        rot.mat.pi.22 <- ToNdarray(
          array(rot.mat.pi[ , 3, 3]), to.ndim = 2, axis = 1)
        rot.mat.pi.diagonal <- cbind(
            rot.mat.pi.00,
            rot.mat.pi.11,
            rot.mat.pi.22)
        a <- apply(rot.mat.pi.diagonal, 1, which.max)[1]
        b <- (a) %% 3 + 1
        c <- (a + 1) %% 3 + 1

        # compute the axis vector
        sq.root <- array(0, dim = c(n.rot.mats, 1))

        aux <- sqrt(
          mask.pi * (
            rot.mat[ , a, a]
            - rot.mat[ , b, b]
            - rot.mat[ , c, c]) + 1)

        sq.root.pi <- array(aux * mask.pi, dim = c(n.rot.mats, 1))

        sq.root <- sq.root + sq.root.pi

        rot.vec.pi <- array(0, dim = c(n.rot.mats, .self$dimension))

        mask.a <- cbind(
        rep(a == 1, n.rot.mats),
        rep(a == 2, n.rot.mats),
        rep(a == 3, n.rot.mats)
        )

        mask.b <- cbind(
          rep(b == 1, n.rot.mats),
          rep(b == 2, n.rot.mats),
          rep(b == 3, n.rot.mats)
        )

        mask.c <- cbind(
          rep(c == 1, n.rot.mats),
          rep(c == 2, n.rot.mats),
          rep(c == 3, n.rot.mats)
        )

        rot.vec.pi <- rot.vec.pi + mask.pi * mask.a * c(sq.root) / 2

        # This avoids division by 0.
        sq.root <- sq.root + mask.0
        sq.root <- sq.root + mask.else

        rot.vec.pi.b <- array(0, c(dim(rot.vec.pi)))
        rot.vec.pi.c <- array(0, c(dim(rot.vec.pi)))

        rot.vec.pi.b <- rot.vec.pi.b + mask.b * c((rot.mat[ , b, a]
                                                   + rot.mat[ , a, b]) / (2 * sq.root))

        rot.vec.pi <- rot.vec.pi + mask.b * mask.pi

        rot.vec.pi.c <- rot.vec.pi.c + mask.c * c((rot.mat[ , c, a]
                                                   + rot.mat[ , a, c]) / (2 * sq.root))

        rot.vec.pi <- rot.vec.pi + mask.pi * mask.pi * rot.vec.pi.c

        norm.rot.vec.pi <- sqrt(sum(rot.vec.pi ^ 2))

        # This avoids division by 0.
        norm.rot.vec.pi <- norm.rot.vec.pi + mask.0
        norm.rot.vec.pi <- norm.rot.vec.pi + mask.else

        rot.vec <- rot.vec + mask.pi * apply(angle * rot.vec.pi, 2,
                                             function(x){x * 1 / norm.rot.vec.pi})

        # This avoid dividing by zero
        angle <- angle + mask.0

        fact <- mask.else * (c(angle) / (2 * sin(c(angle))) - 1)

        rot.vec <- rot.vec * (1 + fact)

      }
      return(.self$Regularize(rot.vec))
    },

    MatrixFromRotationVector = function(rot.vec){
      "Convert rotation vector to rotation matrix."

      stopifnot(.self$Belongs(rot.vec))
      rot.vec <- .self$Regularize(rot.vec)
      n.rot.vecs <- dim(rot.vec)[1]

      if (.self$n == 3) {
        angle <- apply(rot.vec, 1, function(x){sqrt(sum(x ^ 2))})
        angle <- ToNdarray(array(angle), to.ndim = 2, axis = 1)

        skew.rot.vec <- .self$SkewMatrixFromVector(rot.vec)

        coef.1 <- array(0, dim = c(dim(angle)))
        coef.2 <- array(0, dim = c(dim(angle)))

        mask.0 <- (abs(angle) < .Machine$double.eps ^ 0.5)

        coef.1 <- coef.1 + mask.0 * (1 - (angle ^ 2) / 6)
        coef.2 <- coef.2 + mask.0 * (1 / 2 - angle ^ 2)

        mask.else <- !mask.0

        # This avoids division by 0.
        angle <- angle + mask.0

        coef.1 <- coef.1 + mask.else * (sin(angle) / angle)
        coef.2 <- mask.else * ((1 - cos(angle)) / (angle ^ 2))

        term.1 <- array(0, dim = c((n.rot.vecs + .self$n * 2),
                                   (n.rot.vecs + .self$n * 2)))
        term.2 <- term.1

        coef.1 <- c(coef.1)
        coef.2 <- c(coef.2)
        term.1 <- (aperm(replicate(n.rot.vecs, diag(3)), c(3, 2, 1)) + coef.1 * skew.rot.vec)

        squared.skew.rot.vec <- apply(skew.rot.vec, 1, function(x){x %*% x})
        squared.skew.rot.vec <- array(t(squared.skew.rot.vec), dim = c(n.rot.vecs, 3, 3))

        term.2 <- coef.2 * squared.skew.rot.vec

        rot.mat <- term.1 + term.2
      }
      return(rot.mat)
    },

    Compose = function(point.1, point.2){
      "Compose two elements of SO(n)."

      point.1 <- .self$Regularize(point.1)
      point.2 <- .self$Regularize(point.2)

      point.1 <- .self$MatrixFromRotationVector(point.1)
      point.2 <- .self$MatrixFromRotationVector(point.2)

      point.prod <- array(0, dim = dim(point.1))
      for (i in 1:dim(point.1)[1]) {
        point.prod[i, , ] <- point.1[i, , ] %*% point.2[i, , ]
      }
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
    },

    RandomUniform = function(n.samples = 1){
      "Sample in SO(n) with the uniform distribution."
      random.point <- array(runif(n.samples * .self$dimension) * 2 - 1,
                           dim = c(n.samples, .self$dimension))
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
    },

    GroupLog = function(point, base.point=NULL){
      "Compute the group logarithm at point base_point
      of the point point."
      point <- ToNdarray(point, to.ndim = 2)
      base.point <- ToNdarray(base.point, to.ndim = 2)
      point <- .self$Regularize(point)
      base.point <- .self$Regularize(base.point)
      n.points <- dim(point)[1]
      n.base.points <- dim(base.point)[1]
      jacobian <- .self$JacobianTranslation(base.point)
      point.near.id <- .self$Compose(.self$Inverse(base.point), point)
      group.log.from.id <- .self$GroupLogFromIdentity(point.near.id)
      group.log <- array(0, dim = c(n.points, 3))
      for (n in 1:n.points) {
        group.log[n,] <- aperm(jacobian, c(1, 3, 2))[n,,] %*% group.log.from.id[n,]
      }
      return(group.log)
    }
  )
)
