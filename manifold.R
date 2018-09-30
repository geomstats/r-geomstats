#Manifold, i.e. a topological space that locally resembles
#Euclidean space near each point.

library(reticulate)
gs <- import_from_path("geomstats", path = ".")

Manifold <- setRefClass("Manifold",
                        methods = list(

                          init <- function(self, dimension){
                            stopifnot(dimension %% 1 == 0, dimension > 0)
                            self$dimension = dimension
                          },

                          belongs <- function(self, point){
                            # Evaluate if a point belongs to the manifold.

                            # raise NotImplementedError('belongs is not implemented.')
                          },

                          regularize <- function(self, point){
                            # Regularize a point to the canonical representation
                            # chosen for the manifold.

                            # raise NotImplementedError('regularize is not implemented.')
                          }


                        )
                        )

# class Manifold(object):
#   """
# Class for manifolds.
# """
#
# def __init__(self, dimension):
#   assert isinstance(dimension, int) and dimension > 0
# self.dimension = dimension
#
# def belongs(self, point):
#   """
# Evaluate if a point belongs to the manifold.
# """
# raise NotImplementedError('belongs is not implemented.')
#
# def regularize(self, point):
#   """
# Regularize a point to the canonical representation
# chosen for the manifold.
# """
# raise NotImplementedError('regularize is not implemented.')
