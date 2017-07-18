# coding=utf-8

import math
from decimal import Decimal, getcontext

getcontext().prec = 15


class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MEG = 'No unique parallel component to the zero vector'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def plus(self, v):
        new_coordinates = [x + y for x,
                           y in zip(self.coordinates, v.coordinates)]
        # print zip(self.coordinates, v.coordinates)
        #[(8.218, -1.129), (-9.341, 2.111)]
        return Vector(new_coordinates)

    def minus(self, v):
        new_coordinates = [
            x - y for (x, y) in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def scalar_multiply(self, c):
        new_coordinates = [Decimal(c) * x for x in self.coordinates]
        return Vector(new_coordinates)

    def magnitude(self):
        coordinates_squared = [x**Decimal(2) for x in self.coordinates]
        return Decimal(math.sqrt(sum(coordinates_squared)))
        # math.sqrt的结果是float啊

    def normalization(self):  # make a unit vector
        try:
            magnitude = self.magnitude()
            return self.scalar_multiply(Decimal(1.0) / magnitude)
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot_product(self, v):
        return sum([x * y for (x, y) in zip(self.coordinates, v.coordinates)])

    def angel_of_two_vectors(self, v, in_degrees=False):
        try:
            u1 = self.normalization()  # unit vector
            u2 = v.normalization()
            angel_in_radians = math.acos(u1.dot_product(u2))

            if in_degrees:
                degrees_per_radian = 180. / math.pi
                return angel_in_radians * degrees_per_radian
            else:
                return angel_in_radians

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')
            else:
                raise e

    def parallelism(self, v, tolerance=1e-10):
        return (self.magnitude() < tolerance or
                v.magnitude() < tolerance or
                self.angel_of_two_vectors(v) == 0 or
                self.angel_of_two_vectors(v) == math.pi)

    def orthogonality(self, v, tolerance=1e-10):
        return abs(self.dot_product(v)) < tolerance

    def projection_vector(self, b):
        try:
            ub = b.normalization()
            weight = self.dot_product(ub)
            return ub.scalar_multiply(weight)
        except Exception as e: #when basis vector is 0
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MEG)
            else:
                raise e

    def orthogonal_vector(self, b):
        try:
            projection = self.projection_vector(b)
            return self.minus(projection)
        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MEG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MEG)
            else:
                e


v = Vector([3.039, 1.879])
b = Vector([0.825, 2.036])
print v.projection_vector(b)

v = Vector([-9.88, -3.264, -8.159])
b = Vector([-2.155, -9.353, -9.473])
print v.orthogonal_vector(b)

v = Vector([3.009, -6.172, 3.692, -2.51])
b = Vector([6.404, -9.144, 2.759, 8.718])
print v.projection_vector(b)
print v.orthogonal_vector(b)
